using DataFrames, DataFramesMeta
using DocStringExtensions
include("Utils.jl")
import Base:copy


"""
A ``LifeTableDataset` consists of a `DataFrame` for reproduction data, a `DataFrame` for gorwth data, 
and a `DataFrame` for survival data.
$(TYPEDSIGNATURES)
"""
mutable struct LifeTableDataset
    repro::DataFrame
    growth::DataFrame
    survival::DataFrame
    LifeTableDataset() = new()
    function LifeTableDataset(rdf::DataFrame, gdf::DataFrame, sdf::DataFrame)
        ltd = LifeTableDataset()
        ltd.repro = rdf
        ltd.growth = gdf
        ltd.survival = sdf
        return ltd
    end
end

"""
Concatenate vector of LifeTableDatasets into a single LifeTableDataset. Missing values will be skipped.
$(TYPEDSIGNATURES)
"""
function concat(simulations::Vector{Union{Missing,LifeTableDataset}})
    combined_sims = LifeTableDataset(DataFrame(), DataFrame(), DataFrame())
    for sim in simulations[ismissing.(simulations).==false]
        append!(combined_sims.repro, sim.repro)
        append!(combined_sims.growth, sim.growth)
        append!(combined_sims.survival, sim.survival)
    end
    return combined_sims
end

"""
Concatenate vector of LifeTableDatasets into a single LifeTableDataset.
$(TYPEDSIGNATURES)
"""
function concat(simulations::Vector{LifeTableDataset})
    combined_sims = LifeTableDataset(DataFrame(), DataFrame(), DataFrame())
    for sim in simulations
        append!(combined_sims.repro, sim.repro)
        append!(combined_sims.growth, sim.growth)
        append!(combined_sims.survival, sim.survival)
    end
    return combined_sims
end

function subset(ltd::LifeTableDataset, condition::Function)
    subltd = copy(ltd)
    subltd.repro = subltd.repro[condition(subltd.repro),:]
    subltd.growth = subltd.growth[condition(subltd.growth),:]
    subltd.survival = subltd.survival[condition(subltd.survival),:]
    return subltd
end
copy(ltd::LifeTableDataset) = LifeTableDataset(copy(ltd.repro), copy(ltd.growth), copy(ltd.survival))

function control(df::DataFrame, treatment_label::Union{Symbol,String})
    treatment_label = Symbol(treatment_label)
    stressor = df[:,treatment_label]
    mask = stressor.==minimum(stressor)
    return df[mask,:]
end

"""
Add column with control-normalized values.
$(TYPEDSIGNATURES)
"""
function control_normalize(df::D, response::Symbol, treatment_label::Union{String,Symbol}) where D<:AbstractDataFrame
    treatment_label = Symbol(treatment_label)
    #df = df[:,[treatment_label, response, :t_day]]
    ref = combine(groupby(control(df, treatment_label), [:t_day]), response=>mean) # control averages
    ref_col = Symbol(String(response)*"_mean") # column name of the control average
    response_col = Symbol(string(response)*"_norm") # column name of the relative response
    df = leftjoin(df, ref, on=:t_day) # add column with control average
    # iterate over rows and calculate control-normalized values
    y = Float64[]
    for r in eachrow(df)
        r_ref = r[ref_col]
        # if the control average is 0, the response is 1
        if ismissing(r_ref)
            if ismissing(r[response])
                push!(y, missing)
            else
                push!(y, 1.)
            end
        elseif r_ref==0.
            push!(y, 1.)
        # otherwise, calculate response as usual
        else
            push!(y, r[response] / r[ref_col])
        end
    end
    df[!,response_col] = y
    return df
end

"""
Add control normalized values to LifeTableDataset.
$(TYPEDSIGNATURES)
"""
function control_normalize!(ltd::LifeTableDataset, responses::N, treatment_label::Union{Symbol,String}) where N<:NamedTuple
    treatment_label = Symbol(treatment_label)
    ltd.repro = control_normalize(ltd.repro, responses[:repro], treatment_label)
    ltd.growth = control_normalize(ltd.growth, responses[:growth], treatment_label)
end


"""
Apply diffvec to cumulative reproduction, yielding raw brood sizes.
$(TYPEDSIGNATURES)
"""
function add_repro_diffvec!(ltd::LifeTableDataset, treatment_label::Union{String,Symbol})
    # TODO: Can this replace calc_brood_size?
    rpdf = combine(groupby(ltd.repro, [:rep, :metal, Symbol(treatment_label)])) do df
        sort!(df, :t_day)
        DataFrame(
            t_day = df.t_day,
            repro = diffvec(df.cum_repro)
        )
    end 
    ltd.repro = leftjoin(ltd.repro, rpdf, on=[:t_day, :rep, :metal, Symbol(treatment_label)])
end

"""
Brood sizes back-calculated from cumulative reproduction for an individiual experimental unit.
$(TYPEDSIGNATURES)
"""
function calc_brood_size(repro::D) where D<:AbstractDataFrame
    repro = DataFrame(repro)
    repro[!,:brood_size] = vcat([0], diff(repro[:,:cum_repro]))
    repro = @subset(repro, :brood_size.>0)
    instar = 1:nrow(repro)
    brood_size = repro.brood_size
    DataFrame(
        instar=instar,
        brood_size=brood_size
    )
end

"""
Normalize brood numbers.
$(TYPEDSIGNATURES)
"""
function normalize_broods(broods::D, treatment_label::Union{Symbol,String}) where D<:AbstractDataFrame
    treatment_label = Symbol(treatment_label)
    ref_df = broods[broods[:,treatment_label].==minimum(broods[:,treatment_label]),:]
    ref_mean = combine(groupby(ref_df, :instar), :brood_size=>mean)
    broods = leftjoin(broods, ref_mean, on=:instar)
    brood_size_norm = []
    for r in eachrow(broods)
        if ismissing(r[:brood_size_mean])
            push!(brood_size_norm, missing)
        else
            if r[:brood_size_mean]==0
                push!(brood_size_norm, 1)
            else
                push!(brood_size_norm, r[:brood_size]/r[:brood_size_mean])
            end
        end
    end
    broods[!,:brood_size_norm] = brood_size_norm
    return broods
end

"""
Calculate age at first (observed) reproduction.
$(TYPEDSIGNATURES)
"""
function age_at_first_repro(ltd::LifeTableDataset, treatment_label::Union{String,Symbol})
    treatment_label = Symbol(treatment_label)
    afr = combine(groupby(ltd.repro, [:rep, treatment_label])) do df
        let afr
            rdf = df[df.cum_repro .> 0,:]
            if nrow(rdf)>0
                afr = minimum(rdf.t_day)
            else
                afr = Inf
            end
            DataFrame(
                afr = afr
            )
        end
    end
    return afr
end

function age_at_first_repro(repro::DataFrame, treatment_label::Union{String,Symbol})
    treatment_label = Symbol(treatment_label)
    afr = combine(groupby(repro, [:rep, treatment_label])) do df
        let afr
            rdf = df[df.cum_repro .> 0,:]
            if nrow(rdf)>0
                afr = minimum(rdf.t_day)
            else
                afr = Inf
            end
            DataFrame(
                afr = afr
            )
        end
    end
    return afr
end

"""
To DataFrame containing age at first reproduction over treatments, 
add column with control-normalized values.
$(TYPEDSIGNATURES)
"""
function normalize_afr!(afr::D, treatment_label::Union{Symbol,Any}) where D<:AbstractDataFrame
    treatment_label = Symbol(treatment_label)
    ref = mean(skipmissing(afr[afr[:,treatment_label].==minimum(afr[:,treatment_label]),:afr]))
    afr[!,:afr_mean] .= ref
    afr[!,:afr_norm] = afr.afr ./ afr.afr_mean
    select!(afr, Not(:afr_mean))
end

"""
Save a collection of life history statistics.
$(TYPEDSIGNATURES)
"""
function save_life_hist_stats(outdir, life_hist_stats::T) where T<: NamedTuple
    CSV.write(outdir*"broods.csv", life_hist_stats.broods)
    CSV.write(outdir*"repro_21.csv", life_hist_stats.repro21)
    CSV.write(outdir*"afr.csv", life_hist_stats.afr)
    CSV.write(outdir*"length21.csv", life_hist_stats.length21)
end

"""
Save contents of LifeTableDataset to csv. Results in three separate files for reproduction, growth and survival.
$(TYPEDSIGNATURES)
"""
function write_ltd(ltd::LifeTableDataset, path::String)
    CSV.write(path*"_01_repro.csv", ltd.repro)
    CSV.write(path*"_02_growth.csv", ltd.growth)
    CSV.write(path*"_03_survival.csv", ltd.survival)
end

"""
Load DataFrames as LifeTableDataset. Assuming files path*"_repro_csv." etc. to exist.
$(TYPEDSIGNATURES)
"""
function read_ltd(path::String)
    repro = path*"_01_repro.csv"
    growth = path*"_02_growth.csv"
    survival = path*"_03_survival.csv"
    lte = LifeTableDataset(
        CSV.read(repro, DataFrame),
        CSV.read(growth, DataFrame),
        CSV.read(survival, DataFrame)
    )
    #lte.repro.h_z = vectify.(lte.repro.h_z)
    return lte
end

"""
Aggregate life-table dataset by calculating mean and percentiles for each `t_day` and `treatment_label`.
In any case, returns aggregates for absolute values of cumulative reproduction, carapace_length and survival. 
If present, returns aggregates for control-normalized values.
$(TYPEDSIGNATURES)
"""
function aggregate(ltd::LifeTableDataset, treatment_label::Union{Symbol,String}; qntl::Tuple{Float64, Float64}=(0.05, 0.95))
    let repro_agg, growth_agg, survival_agg
        if nrow(ltd.repro)>0
            repro_agg = combine(groupby(drop_na(ltd.repro), [:t_day, Symbol(treatment_label)])) do df
                if ("repro_norm") in names(df)
                    DataFrame(
                        repro_mean = mean(df.repro),
                        repro_median = median(df.repro),
                        repro_lo = quantile(df.repro, qntl[1]),
                        repro_hi = quantile(df.repro, qntl[2]),
                        repro_norm_mean = mean(df.repro_norm),
                        repro_norm_median = median(df.repro_norm),
                        repro_norm_lo = quantile(df.repro_norm, qntl[1]),
                        repro_norm_hi = quantile(df.repro_norm, qntl[2]),
                        cum_repro_mean = mean(df.cum_repro),
                        cum_repro_median = median(df.cum_repro),
                        cum_repro_lo = quantile(df.cum_repro, qntl[1]),
                        cum_repro_hi = quantile(df.cum_repro, qntl[2])
                    )
                elseif ("cum_repro_norm") in names(df)
                    DataFrame(
                        cum_repro_mean = mean(df.cum_repro),
                        cum_repro_median = median(df.cum_repro),
                        cum_repro_lo = quantile(df.cum_repro, qntl[1]),
                        cum_repro_hi = quantile(df.cum_repro, qntl[2]),
                        cum_repro_norm_mean = mean(df.cum_repro_norm),
                        cum_repro_norm_median = median(df.cum_repro_norm),
                        cum_repro_norm_lo = quantile(df.cum_repro_norm, qntl[1]),
                        cum_repro_norm_hi = quantile(df.cum_repro_norm, qntl[2])
                    )
                else
                    DataFrame(
                        cum_repro_mean = mean(df.cum_repro),
                        cum_repro_median = median(df.cum_repro),
                        cum_repro_lo = quantile(df.cum_repro, qntl[1]),
                        cum_repro_hi = quantile(df.cum_repro, qntl[2])
                    )
                end
            end
        end
        if nrow(ltd.growth)>0
            growth_agg = combine(groupby(@subset(drop_na(ltd.growth), isnan.(:carapace_length).==false), [:t_day, Symbol(treatment_label)])) do df
                if ("carapace_length_norm") in names(df)
                    DataFrame(
                        carapace_length_mean = mean(df.carapace_length),
                        carapace_length_median = median(df.carapace_length),
                        carapace_length_lo = quantile(df.carapace_length, qntl[1]),
                        carapace_length_hi = quantile(df.carapace_length, qntl[2]),
                        carapace_length_norm_mean = mean(df.carapace_length_norm),
                        carapace_length_norm_median = median(df.carapace_length_norm),
                        carapace_length_norm_lo = quantile(df.carapace_length_norm, qntl[1]),
                        carapace_length_norm_hi = quantile(df.carapace_length_norm, qntl[2])
                    )
                else
                    DataFrame(
                        carapace_length_mean = mean(df.carapace_length),
                        carapace_length_median = median(df.carapace_length),
                        carapace_length_lo = quantile(df.carapace_length, qntl[1]),
                        carapace_length_hi = quantile(df.carapace_length, qntl[2])
                    )
                end
            end
        end
        if nrow(ltd.survival)>0
            survival_agg = combine(groupby(ltd.survival, [:t_day, Symbol(treatment_label)])) do df
                DataFrame(
                    survival_mean = mean(df.survival),
                    survival_median = median(df.survival),
                    survival_lo = quantile(df.survival, qntl[1]),
                    survival_hi = quantile(df.survival, qntl[2])
                )
            end
        end
        return LifeTableDataset(repro_agg, growth_agg, survival_agg)
    end
end


function add_id_col!(ltd::LifeTableDataset, col::Symbol, val::Any)
    ltd.repro[!,col] .= val
    ltd.growth[!,col] .= val
    ltd.survival[!,col] .= val
end

"""
Aggregate a probabilistic prediction (applies to repro and growth.).
$(TYPEDSIGNATURES)
"""
function agg_prob_predict(ltd::LifeTableDataset, treatment_label::Union{Symbol,String})
    treatment_label = Symbol(treatment_label)

    rpr_agg = combine(groupby(ltd.repro, [:t_day, :n_sample, treatment_label])) do df
        DataFrame(
            cum_repro = mean(df.cum_repro),
            cum_repro_norm = mean(df.cum_repro_norm)
        )
    end

    grg_agg = combine(groupby(ltd.growth, [:t_day, :n_sample, treatment_label])) do df
        DataFrame(
            carapace_length = mean(df.carapace_length),
            carapace_length_norm = mean(df.carapace_length_norm)
        )
    end
    return LifeTableDataset(
        rpr_agg,
        grg_agg,
        ltd.survival
    )
end