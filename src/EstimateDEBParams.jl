# EstimateDEBParams.jl
# High-level functionality to estimate DEB parameters from life-table data.
# Simon Hansul
# 2021-10-22

using Base: Float64
import Base.vcat
using DataFrames, DataFramesMeta
using Plots, StatsPlots

include("LifeTableUtils.jl")

"""
Combine vector of LifeTableDatasets into a single LifeTableDataset. Missing values will be skipped.
$(TYPEDSIGNATURES)
"""
function vcat(simulations::Vector{LifeTableDataset})
    combined_sims = LifeTableDataset(DataFrame(), DataFrame(), DataFrame())
    for sim in simulations
        append!(combined_sims.repro, sim.repro)
        append!(combined_sims.growth, sim.growth)
        append!(combined_sims.survival, sim.survival)
    end
    return combined_sims
end


"""
Load reference data from seperate files. <br>
Currently makes a lot of assumptions about organization of data and column names. Cf. source code in case of doubt.
$(TYPEDSIGNATURES)
"""
function load_calibration_data(
        growth_path::String,
        repro_path::String,  
        survival_path::String,
        survival_replicate_signatures::DataFrame;
        food_levels=[0.4, 0.1],
        treatment_names=["Co", "F"],
    )
    @warn("Assuming column names Length, metal, rep, t_day, repro, cum_repro, num_surviving, num_valid_reps.")
    growth_obs = CSV.read(growth_path, DataFrame)
    growth_obs.Length = growth_obs.Length ./ 10; # conversion from mm to cm
    rename!(growth_obs, Dict("replicate"=>"rep"))
    "tday" in names(growth_obs) ? rename!(growth_obs, Dict("tday"=>"t_day")) : nothing
    select!(growth_obs, :metal, :rep, :t_day, :Length)
    
    repro_obs = CSV.read(repro_path, DataFrame) 
    #repro_obs  = filter(df -> (df.food.=="d"), repro_obs)
    "tday" in names(repro_obs) ? rename!(repro_obs, Dict("tday"=>"t_day")) : nothing
    select!(repro_obs, :metal, :t_day, :rep, :repro, :cum_repro)
    
    timepoints = unique(repro_obs.t_day)
    pushfirst!(timepoints, 0)
    
    survival = @pipe CSV.read(survival_path, DataFrame) |>
    leftjoin(_, survival_replicate_signatures, on=:metal) |>
    @transform(_, :survival=:num_surviving ./ :num_valid_reps) 
    "tday" in names(survival) ? rename!(survival, Dict("tday"=>"t_day")) : nothing
    select!(survival, Not(:food))
    
    calibration_data = [repro_obs, growth_obs, survival]

    for df in calibration_data
        df[!,:food_level] .= 99.
        for (trt,food_level) in zip(treatment_names, food_levels)
            df[df.metal.==trt,:food_level] .= food_level
            df[df.metal.==trt,:food_level] .= food_level
            sort!(df, :food_level)
        end
    end

    return LifeTableDataset(calibration_data...)
end

function transfunct(x)
    return log10(x+1)
end

function scalefunct(x)
    return std(skipmissing(x))
end

function distfunct(a, b)
    return Distances.rmsd(a,b)
end

"""
error function to estimate DEB parameters from growth, reproduction and survival at two food densities.
$(TYPEDSIGNATURES)
"""
function error_function_DEB(
    predicted::LifeTableDataset, 
    calibration_data::LifeTableDataset; 
    return_partial_error=false
    )
    predicted.repro.t_day = round.(predicted.repro.t_day, sigdigits=2)
    predicted.growth.t_day = round.(predicted.growth.t_day, sigdigits=2)
    predicted.survival.t_day = round.(predicted.survival.t_day, sigdigits=2)

    replace!(predicted.growth.carapace_length, NaN=>0)

    scale_repro = scalefunct(calibration_data.repro.cum_repro)
    scale_growth = scalefunct(calibration_data.growth.Length)
    scale_survival = scalefunct(calibration_data.survival.survival)

    if nrow(drop_na(predicted.repro[:,[:cum_repro]]))>0
        error_repro = @pipe outerjoin(predicted.repro, calibration_data.repro, on=[:t_day, :food_level], makeunique=true) |>
        @transform(_, 
        :cum_repro_pred_scaled = :cum_repro ./ scale_repro, 
        :cum_repro_obs_scaled = :cum_repro_1 ./ scale_repro
        ) |>
        drop_na |>
        combine(groupby(_, :food_level)) do df
            if nrow(df)>0
                DataFrame(
                    MSD = distfunct(transfunct.(df.cum_repro_obs_scaled),transfunct.(df.cum_repro_pred_scaled))
                )
            else
                DataFrame(
                    MSD = NaN
                )
            end
        end |> x-> length_corrected_sum(x.MSD)

        error_growth = @pipe outerjoin(predicted.growth, calibration_data.growth, on=[:t_day, :food_level], makeunique=true) |>
        @transform(_, 
        :length_pred_scaled = :carapace_length ./ scale_growth, 
        :length_obs_scaled = :Length ./ scale_growth
        ) |>
        drop_na |>
        combine(groupby(_, :food_level)) do df
            DataFrame(
                MSD = distfunct(transfunct.(df.length_obs_scaled),transfunct.(df.length_pred_scaled))
            )
        end |> x-> length_corrected_sum(x.MSD)

        error_survival = @pipe leftjoin(predicted.survival, calibration_data.survival, on=[:t_day, :food_level], makeunique=true) |>
        @transform(_,
        :survival_pred_scale = :survival ./ scale_survival,
        :survival_obs_scaled = :survival_1 ./ scale_survival
        ) |>
        drop_na |>
        combine(groupby(_, :food_level)) do df
            DataFrame(
                MSD = distfunct(transfunct.(df.survival),transfunct.(df.survival_1))
            )
        end |> x-> length_corrected_sum(x.MSD)
        
        # zeror errors might occur if there is no prediction 
        # -> throw these parmaters out
        if ismissing(error_repro)
            error_repro=99.
        end
        if ismissing(error_growth)
            error_growth=99.
        end
        if ismissing(error_survival)
            error_survival=99.
        end
        if error_repro==0
            error_repro=99.
        end
        if error_growth==0
            error_growth=99.
        end
        if error_survival==0.
            error_survival=99.
        end

        total_error = error_repro + error_growth + error_survival 
        if return_partial_error
            return total_error, (error_repro, error_growth, error_survival) 
        else
            return total_error
        end
    else
        error_survival = @pipe leftjoin(predicted[2], calibration_data[3], on=[:t_day, :food_level], makeunique=true) |>
        drop_na |>
        combine(groupby(_, :food_level)) do df
            DataFrame(
                MSD = distfunct(df.survival, df.survival_1)
            )
        end |> x-> length_corrected_sum(x.MSD)
         # survival error has to be scaled up to account for missing values of three other error values
        total_error = error_survival*3
        if return_partial_error
            return (total_error, (Inf, Inf, error_survival, Inf), false)
        else
            return (total_error, false)
        end
    end
end

function plot_ts(prediction::LifeTableDataset, calibration_data::LifeTableDataset)
    pg = @df calibration_data.growth groupedlineplot(
        :t_day, :Length, :food_level, color=1, layout=(1,2), ylim=(0.025, maximum(:Length)*1.5), title=hcat([LaTeXString("$x mg/d") for x in unique(:food_level)]...), xticks=0:7:21
        )
    @df calibration_data.growth scatter!(:t_day, :Length, group=:food_level, color=1)
    @df prediction.growth groupedlineplot!(:t_day, :carapace_length, :food_level, color=2)
    @df prediction.growth groupedlineplot!(:t_day, :carapace_length, :food_level, (0.25, 0.75), lw=0, color=:red, fillalpha=0.15)

    pr = @df calibration_data.repro groupedlineplot(:t_day, :cum_repro, :food_level, color=1, layout=(1,2), xticks=0:7:21, ylim=(-1,maximum(:cum_repro)*1.5))
    @df calibration_data.repro scatter!(:t_day, :cum_repro, group=:food_level, color=1)
    @df prediction.repro groupedlineplot!(:t_day, :cum_repro, :food_level, color=2)
    @df prediction.repro groupedlineplot!(:t_day, :cum_repro, :food_level, (0.25, 0.75), lw=0, color=:red, fillalpha=0.15)


    ps = @df calibration_data.survival scatter(:t_day, :survival, group=:food_level, color=1, layout=(1,2), xticks=0:7:21, ylim=(-0.05, 1.05), xlabel="Time (d)")
    @df prediction.survival groupedlineplot!(:t_day, :survival, :food_level)
    @df prediction.survival groupedlineplot!(:t_day, :survival, :food_level, (0.25, 0.75), lw=0, color=:red, fillalpha=0.15)

    tsplot = plot(pg, pr, ps, layout=(3, 1), size=(800,600), ylabel=["Length (cm)" "" "Cumulative \n reproduction" "" "Survival" ""]) |>
    pad |>
    x-> plot(x, size=(1000,600))
    return tsplot
end
