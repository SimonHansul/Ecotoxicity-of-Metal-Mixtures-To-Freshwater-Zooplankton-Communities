import StatsBase.mean
using Distributions
using DocStringExtensions, Pipe
using DataFrames, DataFramesMeta, DataStructures
using Plots, StatsPlots
using ProgressMeter


include("DTW.jl")
include("LifeTableUtils.jl")
include("Figures.jl")
include("Utils.jl")

if isdefined(Main, :FITTED_DEB_PARAMS)==false
    const FITTED_DEB_PARAMS =  vcat([
        :h_b
        :J_X_A_m_0
        :kappa_EX_min
        :kappa_0
        :E_G_0
        :p_M_V_0
        :v
        :F_m_0
        :E_H_b_0
        :E_H_p
        :beta_e]...
    )
end

function maketag()
    global tag = species[1:2]*"_"*treatment_label
end

offdist = Truncated(Normal(0, 1e-8), 0, 1e-8) # used to turn off effects

"""
Check if PMoA is engaged in prior distributions.
$(TYPEDSIGNATURES)
"""
function is_engaged(θ::OrderedDict{Symbol,Any}, PMoA::String)
    if abs(mean(θ[Symbol("S_max_"*PMoA)]))>1e-8
        return true
    else
        return false
    end
end

"""
Write prior distributions to csv file.
$(TYPEDSIGNATURES)
"""
function priors_to_file(
    priors::OrderedDict{Symbol,Any}, 
    file::String
    )
    priors_df = DataFrame(
        param=priors.keys,
        dist=[typeof(v) for v in priors.vals],
        params=[[d.untruncated.µ, d.untruncated.σ] for d in priors.vals],
        lower=[d.lower for d in priors.vals],
        upper=[d.upper for d in priors.vals]
        )
    CSV.write(file, priors_df)
end

function exposure_levels_df(df::DataFrame, treatment_label::Union{Symbol,String})
    exposure_levels = DataFrame(
        conc=unique(df[:,treatment_label]),
        level=0:length(unique(df[:,treatment_label]))-1
        )
    rename!(exposure_levels, [treatment_label, "exposure_level"])
    return exposure_levels
end 

add_exposure_levels(df) = leftjoin(df, exposure_levels_df(df, treatment_label), on=treatment_label)

function transfunct(x)
    return log10(x+1)
end

function distfunct(a,b)
    return sum((transfunct.(a)-transfunct.(b))^2)/length(a)
end

function distance(a, b, w)
    #w = w ./ sum(w)
    #return @. (sum(w*(transfunct.(a)-transfunct.(b))^2))/length(a)
    return dtw_cost(Vector{Float64}(a), Vector{Float64}(b), distfunct, 1)
end

"""
If the prediction is Missing or Nothing, return missing.
$(TYPEDSIGNATURES)
"""
function error_function_TKTD(miss::Union{Missing,Nothing}, obs::LifeTableDataset)
    return missing
end

"""
Distance function for estimation of TKTD parameters. May use dynamic time warping (see distance()).
$(TYPEDSIGNATURES)
"""
function error_function_TKTD(
    predicted::LifeTableDataset, 
    observed::LifeTableDataset; 
    return_partial_error=false
    )

    # repro values need to be adjusted for observed time-points
    [select!(predicted.repro, Not(x)) for x in [:repro, :repro_norm, :repro_mean]]
    [select!(predicted.growth, Not(x)) for x in [:carapace_length_norm, :carapace_length_mean]]
    predicted.repro = @subset(predicted.repro, [round(t, sigdigits=2) in unique(observed.repro.t_day) for t in :t_day])

    add_repro_diffvec!(predicted, treatment_label)
    control_normalize!(predicted, (repro=:repro, growth=:carapace_length), treatment_label)

    predicted.growth[(isnan.(predicted.growth.carapace_length)).|(ismissing.(predicted.growth.carapace_length)),:carapace_length] .= 0.
    predicted_agg = aggregate(predicted, :metal)
    
    error_repro = @pipe rightjoin(predicted_agg.repro, observed.repro,
    on=[:t_day, :metal], makeunique=true
    ) |> 
    select(_, :metal, :repro_norm_mean, :repro_norm_mean_1) |>
    drop_na |>
    leftjoin(_, weightsdf, on=:metal) |>
    groupby(_, :metal) |>
    combine(_) do df
        DataFrame(
            dist=distance(df.repro_norm_mean, df.repro_norm_mean_1, df.weight)
        )
    end |> x->length_corrected_sum(x.dist)

    # rightjoin => drop time column
    error_growth = @pipe rightjoin(predicted_agg.growth, observed.growth,
    on=[:metal], makeunique=true
    ) |>
    select(_, :metal, :carapace_length_norm_mean, :carapace_length_norm_mean_1) |>
    drop_na |>
    leftjoin(_, weightsdf, on=:metal) |>
    groupby(_, :metal) |>
    combine(_) do df
        DataFrame(
            dist=distance(df.carapace_length_norm_mean, df.carapace_length_norm_mean_1, df.weight)
        )
    end |> x->length_corrected_sum(x.dist)

    error_survival = @pipe rightjoin(predicted_agg.survival, observed.survival, 
    on=[:t_day, :metal], makeunique=true) |>
    drop_na |>
    leftjoin(_, weightsdf, on=:metal) |>
    combine(groupby(_, :metal)) do df
        DataFrame(
            dist=distance(df.survival_mean, df.survival_mean_1, df.weight)
        )
    end |> x->length_corrected_sum(x.dist)

    errors = [error_repro, error_growth, error_survival]
    errors = errors[(isnan.(errors).==false).&(ismissing.(errors).==false)]

    # sum errors 
    # and correct for different number of error terms if one is missing
    total_error = robust_sum(errors) * (3/length(errors))

    if return_partial_error
        return total_error, [error_repro, error_growth, error_survival]
    else
        return total_error
    end
end

function error_function_TKTD_life_hist_patterns(predicted::LifeTableDataset, observed::NT; return_partial_error=false) where NT<:NamedTuple
    predicted_life_hist_patterns = get_predicted_life_hist_patterns(predicted)

    error_broods = @pipe leftjoin(observed.broods, predicted_life_hist_patterns.broods, on=[:instar, Symbol(treatment_label)], makeunique=true) |>
        select(_, treatment_label, :brood_size_norm, :brood_size_norm_1) |>
        drop_na |>
        combine(groupby(_, treatment_label), x->DataFrame(
            RMSD = Distances.rmsd(Vector{Float64}(log10.(x.brood_size_norm .+1)), Vector{Float64}(log10.(x.brood_size_norm_1 .+1)))
            )) |> x-> length_corrected_sum(x.RMSD)

    error_repro21 = @pipe leftjoin(observed.repro21, predicted_life_hist_patterns.repro21, on=Symbol(treatment_label), makeunique=true) |>
        select(_, treatment_label, :cum_repro_norm, :cum_repro_norm_1) |> 
        drop_na |>
        combine(groupby(_, treatment_label), x->DataFrame(
            RMSD = Distances.rmsd(
                Vector{Float64}(log10.(x.cum_repro_norm .+1)), 
                Vector{Float64}(log10.(x.cum_repro_norm_1 .+1))
                )
            )) |> x-> length_corrected_sum(x.RMSD)

    error_afr = @pipe leftjoin(observed.afr, predicted_life_hist_patterns.afr, on=Symbol(treatment_label), makeunique=true) |>
        select(_, treatment_label, :afr_norm, :afr_norm_1) |>
        drop_na |>
        combine(groupby(_, treatment_label), x->DataFrame(
            RMSD = Distances.rmsd(
                Vector{Float64}(log10.(x.afr_norm .+1)), 
                Vector{Float64}(log10.(x.afr_norm_1 .+1))
            )
            )) |> x-> length_corrected_sum(x.RMSD)
    error_length21 = @pipe leftjoin(observed.length21, predicted_life_hist_patterns.length21, on=Symbol(treatment_label), makeunique=true) |>
        select(_, treatment_label, :carapace_length_norm, :carapace_length_norm_1) |>
        drop_na |> 
        combine(groupby(_, treatment_label), x->DataFrame(
            RMSD = Distances.rmsd(
                Vector{Float64}(log10.(x.carapace_length_norm .+1)),
                Vector{Float64}(log10.(x.carapace_length_norm_1 .+1))
                )
            )) |> x-> length_corrected_sum(x.RMSD)
    error_survival = @pipe leftjoin(observed.surv, predicted_life_hist_patterns.surv, on=Symbol(treatment_label), makeunique=true) |>
        select(_, treatment_label, :survival, :survival_1) |>
        drop_na |>
        combine(groupby(_, treatment_label), x->DataFrame(
            RMSD = Distances.rmsd(
                Vector{Float64}(log10.(x.survival .+1)),
                Vector{Float64}(log10.(x.survival_1 .+1))
            )
        )) |> x-> length_corrected_sum(x.RMSD)

    errors = (error_broods, error_repro21, error_afr, error_length21, error_survival)
    # the errors are weighted so that reproduction, growth and survival have in total the same weight
    total_error = ((error_broods + error_repro21 + error_afr)/3) + error_length21 + error_survival
    if return_partial_error
        return total_error, errors
    else
        return total_error
    end
end


"""
Return quantile of x, or `missing` if no non-missing elements are in x.
$(TYPEDSIGNATURES)
"""
function robust_quantile(x::AbstractArray{Union{R,Missing},1}, q::R) where R<:Real
    x = collect(skipmissing(x))
    if length(x)==0
        return missing
    else
        return quantile(x, q)
    end
end

"""
Read TKTD data, returning LifeTableDataset.
$(TYPEDSIGNATURES)
"""
function read_TKTD_data(
    exposure_path::String, 
    repro_path::String, 
    growth_path::String, 
    survival_path::String, 
    treatment_label::String; 
    replicate_signatures::Union{DataFrame,Nothing}=nothing
    )
    exposure_df = CSV.read(exposure_path, DataFrame)[:,[:metal, Symbol(treatment_label)]]
    exposure_df[:,treatment_label] = round.(exposure_df[:,treatment_label], sigdigits=3)
    exposure = @subset(exposure_df, (occursin.(treatment_label[1:2], :metal)).|(:metal.=="Co"))[:,treatment_label];
    
    growth = CSV.read(growth_path, DataFrame)
    growth = @subset(growth, :metal.!="F", :food.=="D")
    growth = leftjoin(growth, exposure_df, on=:metal)
    growth.Length = growth.Length ./ 10; # conversion from mm to cm
    "Length" in names(growth) ? rename!(growth, Dict(:Length=>:carapace_length)) : nothing
    "tday" in names(growth) ? rename!(growth, Dict(:tday=>:t_day)) : nothing
    growth = control_normalize(growth, :carapace_length, treatment_label)

    repro = CSV.read(repro_path, DataFrame)
    "tday" in names(repro) ? rename!(repro, Dict("tday"=>"t_day")) : nothing
    repro = @subset(repro, :metal.!="F", :food.=="d")
    repro = leftjoin(repro, exposure_df, on=:metal)
    repro = control_normalize(repro, :cum_repro, treatment_label)
    repro = control_normalize(repro, :repro, treatment_label)

    survival = CSV.read(survival_path, DataFrame)
    "tday" in names(survival) ? rename!(survival, Dict("tday"=>"t_day")) : nothing
    survival = @subset(survival, :metal.!="F", :food.=="d")
    survival = leftjoin(survival, exposure_df, on=:metal)
    if "pct_survival" in names(survival)
        survival[!,:survival] = survival.pct_survival ./ 100
    else
        if replicate_signatures==nothing
            error("Need replicate signatures.")
        end
        survival = leftjoin(survival, replicate_signatures, on=:metal)
        survival[!,:survival] = survival.num_surviving ./ survival.num_valid_reps
    end

    sort!(repro, Symbol(treatment_label))
    sort!(growth, Symbol(treatment_label))
    sort!(survival, Symbol(treatment_label))

    return LifeTableDataset(repro, growth, survival), exposure
end

function plot_ts(pred::LifeTableDataset, observed::LifeTableDataset, exposure::Vector; q=(0.05, 0.95))

    pred.repro = @subset(pred.repro, [round(t, sigdigits=2) in unique(observed.repro.t_day) for t in :t_day])

    sort!(pred.repro, :metal)
    sort!(pred.growth, :metal)
    sort!(pred.survival, :metal)
    sort!(observed.repro, :metal)
    sort!(observed.growth, :metal)
    sort!(observed.survival, :metal)

    let lmarg=[4mm 0mm 0mm 0mm 0mm 0mm], yticks=[0, 0.5, 1, 1.5, 2.], plt_r, plt_g, plt_s
        plt_r = @df observed.repro groupedlineplot(
            :t_day, :repro_norm, :metal, q, linetype=:path, estimator=mean,
            layout=(1,6), ylim=(0,2), size=(800,400), color=1, alpha=.5
            )
        @df pred.repro groupedlineplot!(
            :t_day, :repro_norm, :metal, q, estimator=mean,
            color=2, linetype=:path,
            title=hcat(fround.(exposure, sigdigits=2)...),
            ylabel=["Cumulative \n reproduction" "" "" "" "" ""],
            leftmargin=lmarg, yticks=yticks
            )

        if nrow(@subset(observed.growth, :metal.!="Co"))>0
            plt_g = @df observed.growth scatter(
                :t_day, :carapace_length_norm, group=:metal, estimator=mean,
                layout=(1,6), ylim=(0,1.5), size=(800,400)
            )
        else
            plt_g = plot(layout=(1,6))
        end
        @df pred.growth groupedlineplot!(
            :t_day, :carapace_length_norm, :metal, q, estimator=mean,
            layout=(1,5), ylim=(0,1.5), size=(800,400), color=2,
            xlabel="", title="", ylabel=["Carapace \n length" "" "" "" "" ""],
            leftmargin = lmarg,  yticks=yticks
        )

        plt_s = @df observed.survival groupedlineplot(
            :t_day, :survival, :metal, q, linetype=:path, estimator=mean,
            layout=(1,6), ylim=(0,1.01), size=(800,400), title=""
        )
        @df pred.survival groupedlineplot!(
            :t_day, :survival, :metal, q, estimator=mean,
            color=2, linetype=:path,
            title="", bottommargin=4mm, 
            xlabel=["" "" "      Time (d)" "" "" ""], ylabel=["Survival" "" "" "" "" ""],
            lefmargin=lmarg, yticks=yticks
        )

        p = plot(plt_r, plt_g, plt_s, layout=(3,1), xticks=0:7:21, lw=2) |>
        pad |> x->plot(x, size=(1000,500))
        return p
    end
end

function model(TKTD_params::OrderedDict{Symbol,Any})
    m = ""
    for (smax, mchar) in zip([:S_max_G, :S_max_M, :S_max_A, :S_max_R], ['G', 'M', 'A', 'R'])
        if TKTD_params[smax]>0
            m *= mchar
        end
    end 
    return m
end

"""
Generate posterior predictions for TKTD models. 
Takes 1,000 samples by default.
$(TYPEDSIGNATURES)
"""
function TKTD_posterior_prediction(
    accepted_df::DataFrame, 
    parents_df::DataFrame, 
    simulator_function::Function; 
    n_eval=1_000
    )
    error("Function needs to be updated to accomodate ApproxBayes.")
    ppred = LifeTableDataset[]
    ppred_TKTD_samples = DataFrame()
    ppred_DEB_samples = DataFrame()
    for i in 1:n_eval
        try
            TKTD_sample, DEB_sample = posterior_sample(accepted_df, parents_df)
            TKTD_sample_df_row, DEB_sample_df_row = frameify([(TKTD_sample,DEB_sample)])
            TKTD_sample_df_row[!,:n_sample] .= DEB_sample_df_row[!,:n_sample] .= i
            sampled_model = model(TKTD_sample)
            append!(ppred_TKTD_samples, TKTD_sample_df_row)
            append!(ppred_DEB_samples, DEB_sample_df_row)

            lte = simulator_function(TKTD_sample, DEB_sample)
            lte.repro[!,:model] .= sampled_model
            lte.growth[!,:model] .= sampled_model
            lte.survival[!,:model] .= sampled_model
            lte.repro[!,:n_sample] .= i
            lte.growth[!,:n_sample] .= i
            lte.survival[!,:n_sample] .= i
            push!(ppred, lte)
        catch
            @info("Failed to generate predictions.")
        end
    end
    ppred = concat(ppred)
    return ppred, ppred_TKTD_samples, ppred_DEB_samples
end

"""
Calculate benchmarks based on both bayesian predictions and point estimates, using NSE, MAE and RMSD.
$(TYPEDSIGNATURES)
"""
function benchmark(output::T) where T <: NamedTuple
       
    ppred_agg = aggregate(output.ppred, treatment_label)
    ppred_agg.repro.t_day = round.(ppred_agg.repro.t_day, sigdigits=2)
    ppred_agg.growth.t_day = round.(ppred_agg.growth.t_day, sigdigits=2)
    ppred_agg.survival.t_day = round.(ppred_agg.survival.t_day, sigdigits=2)

    comp_repro = rightjoin(ppred_agg.repro, output.observed.repro, on=[:t_day, Symbol(treatment_label)], makeunique=true)
    comp_growth = rightjoin(ppred_agg.growth, output.observed.growth, on=[:t_day, Symbol(treatment_label)], makeunique=true)
    comp_survival = rightjoin(ppred_agg.survival, output.observed.survival, on=[:t_day, Symbol(treatment_label)], makeunique=true)

    comp_repro = comp_repro[:,[:cum_repro_norm_mean, :cum_repro_norm]] |> drop_na
    comp_repro = comp_repro[isnan.(comp_repro.cum_repro_norm_mean).==false,:]
    comp_growth = comp_growth[:,[:carapace_length_norm_mean, :carapace_length_norm]] |> drop_na
    comp_growth = comp_growth[isnan.(comp_growth.carapace_length_norm_mean).==false,:]
    comp_growth.carapace_length_norm_mean = Float64.(comp_growth.carapace_length_norm_mean)

    benchmark_bayesian = DataFrame(
        Variable = ["Cumulative reproduction", "Carapace length", "Survival"],
        Prediction = ["BMA", "BMA", "BMA"],
        NSE = [
            NSE(comp_repro.cum_repro_norm_mean, comp_repro.cum_repro_norm), 
            NSE(comp_growth.carapace_length_norm_mean, comp_growth.carapace_length_norm), 
            NSE(comp_survival.survival_mean, comp_survival.survival)
            ],
        MAE = [
            MAE(comp_repro.cum_repro_norm_mean, comp_repro.cum_repro_norm), 
            MAE(comp_growth.carapace_length_norm_mean, comp_growth.carapace_length_norm), 
            MAE(comp_survival.survival_mean, comp_survival.survival)
            ],
        RMSD = [
            Distances.rmsd(comp_repro.cum_repro_norm_mean, comp_repro.cum_repro_norm), 
            Distances.rmsd(comp_growth.carapace_length_norm_mean, comp_growth.carapace_length_norm), 
            Distances.rmsd(comp_survival.survival_mean, comp_survival.survival)
            ]
    )

    output.ptpred.repro.t_day = round.(output.ptpred.repro.t_day, sigdigits=2)
    output.ptpred.growth.t_day = round.(output.ptpred.growth.t_day, sigdigits=2)
    output.ptpred.survival.t_day = round.(output.ptpred.survival.t_day, sigdigits=2)

    comp_repro = rightjoin(output.ptpred.repro, output.observed.repro, on=[:t_day, Symbol(treatment_label)], makeunique=true)
    comp_growth = rightjoin(output.ptpred.growth, output.observed.growth, on=[:t_day, Symbol(treatment_label)], makeunique=true)
    comp_survival = rightjoin(output.ptpred.survival, output.observed.survival, on=[:t_day, Symbol(treatment_label)], makeunique=true)

    comp_repro = comp_repro[:,[:cum_repro_norm, :cum_repro_norm_1]] |> drop_na
    comp_repro = comp_repro[isnan.(comp_repro.cum_repro_norm).==false,:]
    comp_growth = comp_growth[:,[:carapace_length_norm, :carapace_length_norm_1]] |> drop_na
    comp_growth = comp_growth[isnan.(comp_growth.carapace_length_norm).==false,:]
    comp_growth.carapace_length_norm = Float64.(comp_growth.carapace_length_norm)

    benchmark_pointest = DataFrame(
        Variable = ["Cumulative reproduction", "Carapace length", "Survival"],
        Prediction = ["Point estimate", "Point estimate", "Point estimate"],
        NSE = [
            NSE(comp_repro.cum_repro_norm, comp_repro.cum_repro_norm_1), 
            NSE(comp_growth.carapace_length_norm, comp_growth.carapace_length_norm_1), 
            NSE(comp_survival.survival, comp_survival.survival_1)
            ],
        MAE = [
            MAE(comp_repro.cum_repro_norm, comp_repro.cum_repro_norm_1), 
            MAE(comp_growth.carapace_length_norm, comp_growth.carapace_length_norm_1), 
            MAE(comp_survival.survival, comp_survival.survival_1)
            ],
        RMSD = [
            Distances.rmsd(comp_repro.cum_repro_norm, comp_repro.cum_repro_norm_1), 
            Distances.rmsd(comp_growth.carapace_length_norm, comp_growth.carapace_length_norm_1), 
            Distances.rmsd(comp_survival.survival, comp_survival.survival_1)
            ]
    )

    benchmark_df = vcat(benchmark_bayesian, benchmark_pointest)
    sort!(benchmark_df, :Variable)

    return benchmark_df
end


"""
Evaluate relative contribution of PMoAs (individual sample).
$(TYPEDSIGNATURES)
"""
function evaluate_pmoas(
    posterior::DataFrame, 
    accepted_DEB::DataFrame,
    n_sample::Int64
    )
    let pred_G, pred_M, pred_A, pred_R
        params = posterior_sample(posterior)
        DEB = posterior_sample(accepted_DEB)

        # dataframe to store r-values
        r_df = DataFrame(
            PMoA = String[],
            r_tot = Float64[],
            r_repro = Float64[],
            r_growth = Float64[],
            r_survival = Float64[]
        )
        
        # run complete model + prepare output for use in error function
        pred_full = simulator_function(params, DEB)
        add_id_col!(pred_full, :model, "GMAR")
        add_id_col!(pred_full, :n_sample, n_sample)
        #pred_full = aggregate(pred_full, :metal)
        pred_full.repro.t_day = round.(pred_full.repro.t_day, sigdigits=2)
        pred_full.growth.t_day = round.(pred_full.growth.t_day, sigdigits=2)
        pred_full.survival.t_day = round.(pred_full.survival.t_day, sigdigits=2)
        #rename!(pred_full.repro, Dict(:cum_repro_norm_mean=>:cum_repro_norm))
        #rename!(pred_full.growth, Dict(:carapace_length_norm_mean=>:carapace_length_norm))
        #rename!(pred_full.survival, Dict(:survival_mean=>:survival))

        # run simulation with individual PMoAs engaged
        if abs(params[:S_max_G])>1e-8
            params_G = copy(params)
            params_G[:S_max_M] = params_G[:S_max_A] = params_G[:S_max_R] = params_G[:S_max_EO] = params_G[:h_max] = 0. 
            pred_G = simulator_function(params_G, DEB)
            add_id_col!(pred_G, :model, "G")
            add_id_col!(pred_G, :n_sample, n_sample)

            ρ = error_function_TKTD(
                pred_G,
                aggregate(pred_full, :metal);
                return_partial_error=true
                )
            append!(r_df, DataFrame(
                PMoA="G",
                r_tot = 1/ρ[1],
                r_repro = 1/ρ[2][1],
                r_growth = 1/ρ[2][2],
                r_survival = 1/ρ[2][3]
            ))
        end

        if abs(params[:S_max_M])>1e-8
            params_M = copy(params)
            params_M[:S_max_G] = params_M[:S_max_A] = params_M[:S_max_R] = params_M[:S_max_EO] = params_M[:h_max] = 0.
            pred_M = simulator_function(params_M, DEB)
            add_id_col!(pred_M, :model, "M")
            add_id_col!(pred_M, :n_sample, n_sample)

            ρ = error_function_TKTD(
                pred_M,
                aggregate(pred_full, :metal);
                return_partial_error=true
                )
            append!(r_df, DataFrame(
                PMoA="M",
                r_tot = 1/ρ[1],
                r_repro = 1/ρ[2][1],
                r_growth = 1/ρ[2][2],
                r_survival = 1/ρ[2][3]
            ))
        end

        if abs(params[:S_max_A])>1e-8
            params_A = copy(params)
            params_A[:S_max_G] = params_A[:S_max_M] = params_A[:S_max_R] = params_A[:S_max_EO] = params_A[:h_max] = 0.
            pred_A = simulator_function(params_A, DEB)
            add_id_col!(pred_A, :model, "A")
            add_id_col!(pred_A, :n_sample, n_sample)

            ρ = error_function_TKTD(
                pred_A,
                aggregate(pred_full, :metal);
                return_partial_error=true
                )
            append!(r_df, DataFrame(
                PMoA="A",
                r_tot = 1/ρ[1],
                r_repro = 1/ρ[2][1],
                r_growth = 1/ρ[2][2],
                r_survival = 1/ρ[2][3]
            ))
        end

        if abs(params[:S_max_R])>1e-8
            params_R = copy(params)
            params_R[:S_max_G] = params_R[:S_max_M] = params_R[:S_max_A] = params_R[:S_max_EO] = params_R[:h_max] = 0.
            pred_R = simulator_function(params_R, DEB)
            add_id_col!(pred_R, :model, "R")
            add_id_col!(pred_R, :n_sample, n_sample)

            ρ = error_function_TKTD(
                pred_R,
                aggregate(pred_full, :metal);
                return_partial_error=true
                )
            append!(r_df, DataFrame(
                PMoA="R",
                r_tot = 1/ρ[1],
                r_repro = 1/ρ[2][1],
                r_growth = 1/ρ[2][2],
                r_survival = 1/ρ[2][3]
            ))
        end

        if abs(params[:S_max_EO])>1e-8
            params_EO = copy(params)
            params_EO[:S_max_G] = params_EO[:S_max_M] = params_EO[:S_max_A] = params_EO[:S_max_R] = params_EO[:h_max] = 0.
            pred_EO = simulator_function(params_EO, DEB)
            add_id_col!(pred_EO, :model, "EO")
            add_id_col!(pred_EO, :n_sample, n_sample)

            ρ = error_function_TKTD(
                pred_EO,
                aggregate(pred_full, :metal);
                return_partial_error=true
                )
            append!(r_df, DataFrame(
                PMoA="EO",
                dist_tot = ρ[1],
                dist_repro = ρ[2][1],
                dist_growth = ρ[2][2],
                dist_survival = ρ[2][3]
            ))
        end

        #if params[:h_max]>1e-8
        #    params_H = copy(params)
        #    params_H[:S_max_G] = params_H[:S_max_M] = params_H[:S_max_A] = params_H[:S_max_R] = 0.
        #    pred_H = simulator_function(params_H, DEB)
        #    add_id_col!(pred_H, :model, "H")
        #    add_id_col!(pred_H, :n_sample, n_sample)

        #    ρ = error_function_TKTD(
        #        pred_H,
        #        pred_full;
        #        return_partial_error=true
        #        )
        #    append!(r_df, DataFrame(
        #        PMoA="H",
        #        r_tot = 1/ρ[1],
        #        r_repro = 1/ρ[2][1],
        #        r_growth = 1/ρ[2][2],
        #        r_survival = 1/ρ[2][3]
        #    ))
        #end
        
        r_df[:,:n_sample] .= n_sample
        return r_df
    end
end
"""
Evaluate relative contribution of PMoAs for `n_eval` samples.
$(TYPEDSIGNATURES)
"""
function evaluate_pmoas(posterior::DataFrame, accepted_DEB::DataFrame; n_eval::Int64=100)
    r_df = DataFrame()
    @showprogress for n_sample in 1:n_eval
        try
            r_df_i = evaluate_pmoas(posterior, accepted_DEB, n_sample)
            append!(r_df, r_df_i)
        catch
            @info("Failed to evaluate sample $n_sample.")
        end
    end
    return r_df
end

"""
Plot effect on brood sizes.
$(TYPEDSIGNATURES)
"""
function plot_broods(
    broods_obs::DataFrame,
    broods_pred_agg::DataFrame;
    conc_label="Concentration",
    colors=[1,2],
    ylim=(0,1.5)
    )
    x = broods_obs[:,treatment_label]
    y = broods_obs[:,:brood_size_norm]
    z = broods_obs[:,:instar]

    maxx = maximum(x)

    local p = groupedlineplot(
        x, y, z, (.25, .75),
        xscale=:log10, fillalpha=.25,
        layout=(1,length(unique(z))), 
        marker=true,lw=1.5,
        title=hcat(["AI $x" for x in unique(z)]...),
        xlabel=conc_label,
        ylabel=hcat(vcat(["Control-normalized \n brood size"], repeat([""], length(unique(z))-1))...),
        xrotation=45, ylim=ylim,
        legend=hcat(vcat(repeat([false], length(unique(z))-1), [:topright])...),
        label="Observed",  foreground_color_legend = nothing, color=:black
        ) 
    x = @subset(broods_pred_agg, [b in unique(broods_obs.instar) for b in :instar])[:,treatment_label]
    y = @subset(broods_pred_agg, [b in unique(broods_obs.instar) for b in :instar])[:,:predicted_mean]
    z = @subset(broods_pred_agg, [b in unique(broods_obs.instar) for b in :instar])[:,:instar]
    
    maxx = max(maxx, maximum(x))

    plot!(
        x, y, group=z, 
        xlim=(minimum(x)*.5, maxx*1.5),
        xticks=(10 .^ round.(log10.(unique(x)), sigdigits=2), fround.(unique(x), sigdigits=2)),
        marker=:diamond, lw=1.5, label="Fitted", color=colors[2]
        )
    p = p |> pad |> x->plot(x, size=(1500,350), bottommargin=5mm, leftmargin=2mm, background=:transparent)
    return p
end

"""
Plot effect on brood sizes, observed + predicted from two model versions.
$(TYPEDSIGNATURES)
"""
function plot_broods(
    broods_obs::DataFrame,
    broods_pred_agg_parsim::DataFrame,
    broods_pred_agg_gmar::DataFrame;
    conc_label="Concentration",
    model_tags=["parsimonious", "combined-PMoA"],
    ylim=(0,1.5)
    )

    ml1, ml2 = model_tags

    # plot observed
    x = broods_obs[:,treatment_label]
    y = broods_obs[:,:brood_size_norm]
    z = broods_obs[:,:instar]

    maxx = maximum(x)

    local p = groupedlineplot(
        x, y, z, (.25, .75),
        xscale=:log10, fillalpha=.25,
        layout=(1,length(unique(z))), 
        marker=true,lw=1.5,
        title=hcat(["AI $x" for x in unique(z)]...),
        xlabel=conc_label,
        ylabel=hcat(vcat(["Control-normalized \n brood size"], repeat([""], length(unique(z))-1))...),
        xrotation=45, ylim=ylim,
        legend=hcat(vcat(repeat([false], length(unique(z))-1), [:topright])...),
        label="Observed",  foreground_color_legend = nothing, color=:black
        ) 
    
    # plot predicted (parsimonious)
    x = @subset(broods_pred_agg_parsim, [b in unique(broods_obs.instar) for b in :instar])[:,treatment_label]
    y = @subset(broods_pred_agg_parsim, [b in unique(broods_obs.instar) for b in :instar])[:,:predicted_mean]
    z = @subset(broods_pred_agg_parsim, [b in unique(broods_obs.instar) for b in :instar])[:,:instar]

    maxx = max(maxx, maximum(x))

    plot!(
        x, y, group=z, 
        marker=:diamond, lw=1.5, label="Fitted ($ml1)", color=2
        )

    # plot predicted (gmar)
    x = @subset(broods_pred_agg_gmar, [b in unique(broods_obs.instar) for b in :instar])[:,treatment_label]
    y = @subset(broods_pred_agg_gmar, [b in unique(broods_obs.instar) for b in :instar])[:,:predicted_mean]
    z = @subset(broods_pred_agg_gmar, [b in unique(broods_obs.instar) for b in :instar])[:,:instar]

    maxx = max(maxx, maximum(x))

    plot!(
        x, y, group=z, 
        xlim=(minimum(x)*.5, maxx*1.5),
        xticks=(10 .^ round.(log10.(unique(x)), sigdigits=2), fround.(unique(x), sigdigits=2)),
        marker=:diamond, lw=1.5, label=ml2, color=3
        )
    p = p |> pad |> x->plot(x, size=(1500,350), bottommargin=5mm, leftmargin=2mm, background=:transparent)

    return p
end

function aggregate_broods(broods_pred::D, treatment_label::Union{Symbol,String}) where D<:AbstractDataFrame

    broods_pred_agg = combine(groupby(broods_pred, [:instar, Symbol(treatment_label)])) do df
        
        DataFrame(
            predicted_mean = mean(df.brood_size_norm),
            predicted_median = median(df.brood_size_norm),
            predicted_q25 = quantile(df.brood_size_norm, 0.25),
            predicted_q75 = quantile(df.brood_size_norm, 0.75)
        )
    end |> add_exposure_levels
    return broods_pred_agg
end

"""
Evaluate a parsimonious model containing a subset of mutiple engaged PMoAs. 
`pmoas` is a Vector containing one or pmoa suffixes. Returns simulation output, parameter samples and associated errors.
"""
function evaluate_pmoas(posterior, accepted_DEB, pmoas::Vector{String}; n_eval=100)
    # Input validation
    if sum([(x in ["G", "M", "A", "R", "EO", "KP", "H"])==false for x in pmoas])>0
        error("Uknown PMoA contained in $pmoas.")
    end

    simulations = LifeTableDataset[]
    samples = DataFrame()
    errors = []
    
    disengaged_pmoas = ["G", "M", "A", "R", "EO", "KP", "H"] |> x-> filter(v->occursin(pmoas)==false, x) #x->x[occursin.(pmoas).==false]

    for n in 1:n_eval
        try
            params = posterior_sample(posterior)
            [turn_off!(params, x) for x in disengaged_pmoas]
            pred = simulator_function(params, posterior_sample(accepted_DEB))
            
            pred.repro[!,:n_sample] .= n
            pred.growth[!,:n_sample] .= n
            pred.survival[!,:n_sample] .= n
            
            push!(simulations, pred)
            append!(samples, DataFrame(hcat(params.vals...), params.keys))
            push!(errors, error_function_TKTD(pred, calibration_data))

        catch
            
            @info("Failed to generate prediction no. $n.")
        end
    end
    return concat(simulations), samples, errors
end

function identify_pmoa(sample::OrderedDict{Symbol,Any})
    pmoa=""
    for (sbstr, smaxpar) in zip([
        "G", "M", "A", "R", "EO", "KP"
    ], [
        :S_max_G, :S_max_M, :S_max_A, :S_max_R, :S_max_EO
    ])
        if sample[smaxpar]>1e-3
            pmoa=pmoa*sbstr*"-"
        end
    end
    if (length(pmoa)>0) && (pmoa[end]=='-')
        pmoa=pmoa[1:end-1]
    end
    if length(pmoa)==0
        pmoa="none"
    end

    return pmoa
end

function identify_pmoa(simulation_output::DataFrame)
    pmoa_vec = []
    for row in eachrow(simulation_output)
        pmoa=""
        for (sbstr, smaxpar) in zip([
            "G", "M", "A", "R", "EO", "KP"
        ], [
            :S_max_G, :S_max_M, :S_max_A, :S_max_R, :S_max_EO
        ])
            if row[smaxpar]>1e-3
                pmoa=pmoa*sbstr*"-"
            end
        end
        if (length(pmoa)>0) && (pmoa[end]=='-')
            pmoa=pmoa[1:end-1]
        end
        if length(pmoa)==0
            pmoa="none"
        end
        push!(pmoa_vec, pmoa)
    end
    return pmoa_vec
end

function plot_dists_gmar(
    simulation_output::DataFrame, 
    accepted::DataFrame;
    pmoa_labs=OrderedDict(
        "G"=>"Growth \n costs",
        "M"=>"Maintenance \n costs",
        "A"=>"Assimilation \n efficiency",
        "R"=>"Reproduction \n efficiency",
        "EO"=>"Embryonal \n investment"
    )
    )
    
    kdplots = []
    kdplottitles = []
    countvar=0
    for pmoa in ["G", "M", "A", "R", "EO"]
        let p=pad()
            if nrow(@subset(simulation_output, occursin.(pmoa, :pmoa)))>0
                countvar+=1
                push!(kdplottitles, pmoa_labs[pmoa])
                par = Symbol("k_d_"*pmoa)
                p = histogram(
                    simulation_output[:,par], 
                    normalize=:pdf, fillalpha=.25, lw=0, fill=true,
                    xlabel=L"\dot{k}_D\ (d^{-1})", 
                    ylabel=countvar==1 ? "Probability \n density" : "",
                    legend=countvar==1 ? true : false,
                    label="Prior"
                    )
                histogram!(accepted[:,par], normalize=:pdf, fillalpha=.25, lw=0, label="Posterior")
            end
            push!(kdplots, p)
        end
    end
    filter!(x->x.n>0, kdplots)
    ph = histogram(
        simulation_output[:,:k_d_H], normalize=:pdf, lw=0, fill=true, fillalpha=.25, 
        xlabel=L"\dot{k}_D\ (d^{-1})", 
        ylabel=countvar==0 ? "Probability \n density" : "",
        leg=countvar==0 ? true : false, label="Prior"
        )
    histogram!(accepted[:,:k_d_H], weights=Weights(accepted[:,:weights]), lw=0, normalize=:pdf, fillalpha=.25, label="Posterior")
    push!(kdplots, ph)
    push!(kdplottitles, "Survival")

    kdplot = plot(
        kdplots..., 
        layout=(1,length(kdplots)), 
        size=(1200*(length(kdplots)/4),250), 
        leftmargin=hcat(vcat([10mm],repeat([5mm], length(kdplots)-1))...),
        title=hcat(kdplottitles...), titlefontsize=12
    )

    ed50plots = []
    countvar=0
    for pmoa in ["G", "M", "A", "R", "EO"]
        let p=pad()
            if nrow(@subset(simulation_output, occursin.(pmoa, :pmoa)))>0
                countvar+=1
                par = Symbol("S_ED50_"*pmoa)
                p = histogram(
                    log10.(simulation_output[:,par]), 
                    normalize=:pdf, fillalpha=.25, lw=0, fill=true,
                    xlabel=L"log_{10}(ED50)\ (nM)", 
                    ylabel=countvar==1 ? "Probability \n density" : ""
                    )
                histogram!(log10.(accepted[:,par]), weights=Weights(accepted[:,:weights]), normalize=:pdf, fillalpha=.25, lw=0)
            end
            push!(ed50plots, p)
        end
    end
    filter!(x->x.n>0, ed50plots)
    ph = histogram(
        log10.(simulation_output[:,:h_ED50]), normalize=:pdf, lw=0, fill=true, fillalpha=.25, 
        xlabel=L"log_{10}(ED50)\ (nM)",
        ylabel=countvar==0 ? "Probability \n density" : ""
        )
    histogram!(log10.(accepted[:,:h_ED50]), weights=Weights(accepted[:,:weights]), lw=0, normalize=:pdf, fillalpha=.25)
    push!(ed50plots, ph)

    ed50plot = plot(
        ed50plots..., 
        layout=(1,length(ed50plots)), 
        size=(1200*(length(ed50plots)/4),250),  
        leftmargin=hcat(vcat([10mm],repeat([5mm], length(ed50plots)-1))...)
    )

    betaplots = []
    countvar=0
    for pmoa in ["G", "M", "A", "R", "EO"]
        let p=pad()
            if nrow(@subset(simulation_output, occursin.(pmoa, :pmoa)))>0
                countvar+=1
                par = Symbol("S_beta_"*pmoa)
                p = histogram(
                    simulation_output[:,par], 
                    normalize=:pdf, fillalpha=.25, lw=0, fill=true,
                    xlabel=L"\beta", 
                    ylabel=countvar==1 ? "Probability \n density" : ""
                    )
                histogram!(accepted[:,par], weights=Weights(accepted[:,:weights]), normalize=:pdf, fillalpha=.25, lw=0)
            end
            push!(betaplots, p)
        end
    end
    
    filter!(x->x.n>0, betaplots)
    ph = histogram(
        simulation_output[:,:h_beta], normalize=:pdf, lw=0, fill=true, fillalpha=.25, 
        xlabel=L"\beta",
        ylabel=countvar==0 ? "Probability \n density" : ""
        )
    histogram!(accepted[:,:h_beta], weights=Weights(accepted[:,:weights]), lw=0, normalize=:pdf, fillalpha=.25)
    push!(betaplots, ph)

    betaplot = plot(
        betaplots..., 
        layout=(1,length(betaplots)), 
        size=(1200*(length(betaplots)/4),250), 
        leftmargin=hcat(vcat([10mm],repeat([5mm], length(betaplots)-1))...)
    )
    
    smaxplots = []
    countvar=0
    for pmoa in ["G", "M", "A", "R", "EO"]
        let p=pad()
            if nrow(@subset(simulation_output, occursin.(pmoa, :pmoa)))>0
                countvar+=1
                par = Symbol("S_max_"*pmoa)
                if species=="dp3"
                    @warn("Plotting s_max distributions only for D. pulex (dp3). Species name hardcoded!")
                    p = histogram(
                        @subset(simulation_output, occursin.(pmoa, :pmoa))[:,par], 
                        normalize=:pdf, fillalpha=.25, lw=0, fill=true,
                        xlabel=L"s_{max}", 
                        ylabel=countvar==1 ? "Probability \n density" : ""
                        )
                    histogram!(@subset(accepted, occursin.(pmoa, :pmoa))[:,par], weights=Weights(1 ./@subset(accepted, occursin.(pmoa, :pmoa))[:,:distance]), normalize=:pdf, fillalpha=.25, lw=0)   
                end
                push!(smaxplots, p)
            end
        end
    end

    #filter!(x->x.n>0, smaxplots) # remove empty plots
    ph = histogram(
        simulation_output[:,:h_max], 
        normalize=:pdf, lw=0, fill=true, fillalpha=.25, 
        xlabel=L"\dot{h}_{max}\ (d^{-1})", 
        ylabel=countvar==0 ? "Probability \n density" : "",
        )
    countvar=0
    for pmoa in ["G", "M", "A", "R", "EO", "KP"]
        if nrow(@subset(accepted, :pmoa.==pmoa))>0
            countvar+=1
            histogram!(
                @subset(accepted, occursin.(pmoa, :pmoa))[:,:h_max],
                 weights=Weights(1 ./@subset(accepted, occursin.(pmoa, :pmoa))[:,:distance]), 
                 lw=0, normalize=:pdf, fillalpha=.15, label=pmoa
                 )
        end
    end
    if countvar==0
        histogram!(accepted[:,:h_max],
                 weights=Weights(1 ./accepted[:,:distance]), 
                 lw=0, normalize=:pdf, fillalpha=.15, leg=false
                 )
    end
    push!(smaxplots, ph)

    smaxplot = plot(
        smaxplots..., 
        layout=(1,length(smaxplots)), 
        size=(1200*(length(betaplots)/4),250), 
        leftmargin=hcat(vcat([10mm],repeat([5mm], length(smaxplots)-1))...)
    )

    distplot = plot(
        kdplot, ed50plot, betaplot, smaxplot, layout=(4,1), 
        size=(1200*(length(betaplots)/4),650),
        bottommargin=12mm
        )

    return distplot

end

function plot_dists_parsim(
    simulation_output::DataFrame, 
    accepted::DataFrame;
    pmoa_labs=OrderedDict(
        "G"=>"Growth \n costs",
        "M"=>"Maintenance \n costs",
        "A"=>"Assimilation \n efficiency",
        "R"=>"Reproduction \n efficiency",
        "EO"=>"Embryonal \n investment"
    )
    )
    kdplots = []
    kdplottitles = []
    countvar=0
    for pmoa in ["G", "M", "A", "R", "EO"]
        let p=pad(), simout_df=@subset(simulation_output, occursin.(pmoa, :pmoa)), acc_df=@subset(accepted, occursin.(pmoa, :pmoa))
            if nrow(simout_df)>0
                countvar+=1
                push!(kdplottitles, pmoa_labs[pmoa])
                par = Symbol("k_d_"*pmoa)
                # plot prior samples
                p = histogram(
                    simout_df[:,par], 
                    normalize=:pdf, fillalpha=.25, lw=0, fill=true,
                    xlabel=L"\dot{k}_D\ (d^{-1})", 
                    ylabel=countvar==1 ? "Probability \n density" : ""
                    )
                # plot accepted parameter values
                if nrow(acc_df)>0
                    histogram!(acc_df[:,par], weights=Weights(1 ./acc_df[:,:distance]), normalize=:pdf, fillalpha=.25, lw=0)
                end
            end
            push!(kdplots, p)
        end
    end
    filter!(x->x.n>0, kdplots)
    ph = histogram(
        simulation_output[:,:k_d_H], normalize=:pdf, lw=0, fill=true, fillalpha=.25, 
        xlabel=L"\dot{k}_D\ (d^{-1})", 
        ylabel=countvar==0 ? "Probability \n density" : "",
        leg=countvar==0 ? true : false, label="Prior"
        )
    for pmoa in ["G", "M", "A", "R", "EO", "KP"]
        if nrow(@subset(accepted, occursin.(pmoa, :pmoa)))>0
            histogram!(@subset(accepted, occursin.(pmoa, :pmoa))[:,:k_d_H], weights=Weights(1 ./@subset(accepted, occursin.(pmoa, :pmoa))[:,:distance]), lw=0, normalize=:pdf, fillalpha=.25, label=pmoa)
        end
    end
    push!(kdplots, ph)
    push!(kdplottitles, "Survival")

    kdplot = plot(
        kdplots..., 
        layout=(1,length(kdplots)), 
        size=(1200*(length(kdplots)/4),250), 
        leftmargin=hcat(vcat([10mm],repeat([5mm], length(kdplots)-1))...),
        title=hcat(kdplottitles...), 
        titlefontsize=12
    )

    ed50plots = []
    countvar=0
    for pmoa in ["G", "M", "A", "R", "EO"]
        let p=pad()
            if nrow(@subset(simulation_output, occursin.(pmoa, :pmoa)))>0
                countvar+=1
                par = Symbol("S_ED50_"*pmoa)
                p = histogram(
                    log10.(@subset(simulation_output, occursin.(pmoa, :pmoa))[:,par]), 
                    normalize=:pdf, fillalpha=.25, lw=0, fill=true,
                    xlabel=L"log_{10}(ED50)\ (nM)", 
                    ylabel=countvar==1 ? "Probability \n density" : ""
                    )
                if nrow(@subset(accepted, occursin.(pmoa, :pmoa)))>0
                    histogram!(log10.(@subset(accepted, occursin.(pmoa, :pmoa))[:,par]), weights=Weights(1 ./@subset(accepted, occursin.(pmoa, :pmoa))[:,:distance]), normalize=:pdf, fillalpha=.25, lw=0)
                end
            end
            push!(ed50plots, p)
        end
    end
    filter!(x->x.n>0, ed50plots)
    ph = histogram(
        log10.(simulation_output[:,:h_ED50]), normalize=:pdf, lw=0, fill=true, fillalpha=.25, leg=true, label="",
        xlabel=L"log_{10}(ED50)\ (nM)", 
        ylabel=countvar==0 ? "Probability \n density" : ""
        )
    for pmoa in ["G", "M", "A", "R", "EO", "KP"]
        if nrow(@subset(accepted, occursin.(pmoa, :pmoa)))>0
            histogram!(log10.(@subset(accepted, occursin.(pmoa, :pmoa))[:,:h_ED50]), weights=Weights(1 ./@subset(accepted, occursin.(pmoa, :pmoa))[:,:distance]), lw=0, normalize=:pdf, fillalpha=.25, label=pmoa)
        end
    end
    push!(ed50plots, ph)

    ed50plot = plot(
        ed50plots..., 
        layout=(1,length(ed50plots)), 
        size=(1200*(length(ed50plots)/4),250), 
        leftmargin=hcat(vcat([10mm],repeat([5mm], length(ed50plots)-1))...)
    )

    betaplots = []
    countvar=0
    for pmoa in ["G", "M", "A", "R", "EO"]
        let p=pad()
            if nrow(@subset(simulation_output, occursin.(pmoa, :pmoa)))>0
                countvar+=1
                par = Symbol("S_beta_"*pmoa)
                p = histogram(
                    @subset(simulation_output, occursin.(pmoa, :pmoa))[:,par], 
                    normalize=:pdf, fillalpha=.25, lw=0, fill=true,
                    xlabel=L"\beta", 
                    ylabel=countvar==1 ? "Probability \n density" : ""
                    )
                if nrow(@subset(accepted, occursin.(pmoa, :pmoa)))>0
                    histogram!(@subset(accepted, occursin.(pmoa, :pmoa))[:,par], weights=Weights(1 ./@subset(accepted, occursin.(pmoa, :pmoa))[:,:distance]), normalize=:pdf, fillalpha=.25, lw=0)
                end
            end
            push!(betaplots, p)
        end
    end
    filter!(x->x.n>0, betaplots)
    ph = histogram(
        simulation_output[:,:h_beta], normalize=:pdf, lw=0, fill=true, fillalpha=.25, leg=true, label="",
        xlabel=L"\beta", 
        ylabel=countvar==0 ? "Probability \n density" : ""
        )
    for pmoa in ["G", "M", "A", "R", "EO", "KP"]
        if nrow(@subset(accepted, occursin.(pmoa, :pmoa)))>0
            histogram!(@subset(accepted, occursin.(pmoa, :pmoa))[:,:h_beta], weights=Weights(1 ./@subset(accepted, occursin.(pmoa, :pmoa))[:,:distance]), lw=0, normalize=:pdf, fillalpha=.25, label=pmoa)
        end
    end
    push!(betaplots, ph)

    betaplot = plot(
        betaplots..., 
        layout=(1,length(betaplots)), 
        size=(1200*(length(betaplots)/4),250), 
        leftmargin=hcat(vcat([10mm],repeat([5mm], length(betaplots)-1))...)
    )
    
    smaxplots = []
    countvar=0
    for pmoa in ["G", "M", "A", "R", "EO"]
        let p=pad()
            if nrow(@subset(simulation_output, occursin.(pmoa, :pmoa)))>0
                countvar+=1
                par = Symbol("S_max_"*pmoa)
                if species=="dp3"
                    @warn("Plotting s_max distributions only for D. pulex (dp3). Species name hardcoded!")
                    p = histogram(
                        @subset(simulation_output, occursin.(pmoa, :pmoa))[:,par], 
                        normalize=:pdf, fillalpha=.25, lw=0, fill=true,
                        xlabel=L"s_{max}", 
                        ylabel=countvar==1 ? "Probability \n density" : ""
                        )
                    if nrow(@subset(accepted, occursin.(pmoa, :pmoa)))>0
                        histogram!(@subset(accepted, occursin.(pmoa, :pmoa))[:,par], weights=Weights(1 ./@subset(accepted, occursin.(pmoa, :pmoa))[:,:distance]), normalize=:pdf, fillalpha=.25, lw=0)   
                    end
                end
                push!(smaxplots, p)
            end
        end
    end

    #filter!(x->x.n>0, smaxplots) # remove empty plots
    ph = histogram(
        simulation_output[:,:h_max], 
        normalize=:pdf, lw=0, fill=true, fillalpha=.25, 
        xlabel=L"\dot{h}_{max}\ (d^{-1})", leg=true, label="",
        ylabel=countvar==0 ? "Probability \n density" : "",
        )
    for pmoa in ["G", "M", "A", "R", "EO", "KP"]
        if nrow(@subset(accepted, occursin.(pmoa, :pmoa)))>0
            histogram!(@subset(accepted, occursin.(pmoa, :pmoa))[:,:h_max], weights=Weights(1 ./@subset(accepted, occursin.(pmoa, :pmoa))[:,:distance]), lw=0, normalize=:pdf, fillalpha=.25, label=pmoa)
        end
    end
    push!(smaxplots, ph)

    smaxplot = plot(
        smaxplots..., 
        layout=(1,length(smaxplots)), 
        size=(1200*(length(betaplots)/4),250), 
        leftmargin=hcat(vcat([10mm],repeat([5mm], length(smaxplots)-1))...)
    )

    distplot = plot(
        kdplot, ed50plot, betaplot, smaxplot, layout=(4,1), 
        size=(1200*(length(betaplots)/4),650), 
        bottommargin=12mm
        )
    return distplot
end

"""
Quntitative evaluation based on nNSE, PPC and MAE.
$(TYPEDSIGNATURES)
"""
function quanteval(ppred::LifeTableDataset, calibration_data::LifeTableDataset)
    eval_repro = @pipe aggregate(ppred, :metal).repro |> 
    leftjoin(_, aggregate(calibration_data, :metal).repro, on=[:t_day, :metal], makeunique=true) |> drop_na

    cum_repro_nNSE = nNSE(Vector{Float64}(eval_repro.cum_repro_mean), Vector{Float64}(eval_repro.cum_repro_mean_1))
    cum_repro_MAE = MAE(eval_repro.cum_repro_mean, eval_repro.cum_repro_mean_1)
    cum_repro_PPC = PPC(ppred, calibration_data, :repro, :cum_repro_lo, :cum_repro_hi, :cum_repro)

    repro_norm_nNSE = nNSE(Vector{Float64}(eval_repro.repro_norm_mean), Vector{Float64}(eval_repro.repro_norm_mean_1))
    repro_norm_MAE = MAE(Vector{Float64}(eval_repro.repro_norm_mean), Vector{Float64}(eval_repro.repro_norm_mean_1))
    repro_norm_PPC = PPC(ppred, calibration_data, :repro, :repro_norm_lo, :repro_norm_hi, :repro_norm)

    ppred.growth = @subset(ppred.growth, isnan.(:carapace_length).==false)
    eval_growth = @pipe aggregate(ppred, :metal).growth |> 
    leftjoin(_, aggregate(calibration_data, :metal).growth, on=[:t_day, :metal], makeunique=true) |> drop_na
    if nrow(@subset(eval_growth, :metal.!="Co"))>0
        global growth_nNSE = nNSE(Vector{Float64}(eval_growth.carapace_length_norm_mean), Vector{Float64}(eval_growth.carapace_length_norm_mean))
        global growth_MAE = MAE(eval_growth.carapace_length_norm_mean, eval_growth.carapace_length_norm_mean_1)
        global growth_PPC = PPC(ppred, calibration_data, :growth, :carapace_length_norm_lo, :carapace_length_norm_hi, :carapace_length_norm)
    else
        global growth_nNSE = "N.M."
        global growth_MAE = "N.M."
        global growth_PPC = "N.M."
    end

    eval_surv = @pipe aggregate(ppred, :metal).survival |> 
    leftjoin(_, calibration_data.survival, on=[:t_day, :metal], makeunique=true) |> drop_na 

    surv_nNSE = nNSE(eval_surv.survival_mean, eval_surv.survival)
    surv_MAE = MAE(eval_surv.survival_mean, eval_surv.survival)
    surv_PPC = PPC(ppred, calibration_data, :survival, :survival_lo, :survival_hi, :survival)

    speclabs = Dict(
        "dl"=>"D. longispina",
        "dm"=>"D. magna",
        "dp"=>"D. pulex"
    )

    quant_eval = DataFrame(
        species = speclabs[species[1:2]],
        metal = treatment_label[1:2],
        variable = ["Cumulative reproduction", "Control-normalized brood size", "Control-normalized carapace length", "Survival"],
        nNSE = [cum_repro_nNSE, repro_norm_nNSE, growth_nNSE, surv_nNSE],
        MAE = [cum_repro_MAE, repro_norm_MAE, growth_MAE, surv_MAE],
        PPC = [cum_repro_PPC, repro_norm_PPC, growth_PPC, surv_PPC]
    )
    return quant_eval
end