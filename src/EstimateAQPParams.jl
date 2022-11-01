# EstimateAQPParams.jl
# High-level functionality to estimate AQP parameters from 10-day batchtest data.
# Simon Hansul
# 2021-06-04

using DataFrames, DataFramesMeta
using Pipe


function plot_timeseries(valid_df, calibration_data)
    fig, ax = plt.subplots(ncols=3, figsize=(10,4))
    @with valid_df begin
        sns.lineplot(:t_day, :A, ci="sd", label="predicted", ax=ax[1])
        sns.lineplot(:t_day, :Q, ci="sd", ax=ax[2])
        sns.lineplot(:t_day, :P, ci="sd", ax=ax[3])
    end

    # only for visualization purposes: replace missing values with Inf
    calibration_data[ismissing.(calibration_data.P),:P] .= Inf
    @with calibration_data begin
        sns.scatterplot(:t_day, :A, ci="sd", label="observed", marker="o", ax=ax[1])
        sns.scatterplot(:t_day, :Q, ci="sd", ax=ax[2])
        sns.scatterplot(:t_day, :P, ci="sd", ax=ax[3])
    end
    ax[1].set(xlabel="time (days)", ylabel="A (mg dwt)")
    ax[2].set(xlabel="time (days)", ylabel="Q (mg dwt)")
    ax[3].set(xlabel="time (days)", ylabel="P (mg dwt)")
    sns.despine()
    plt.tight_layout()
end

function loss_function_AQP(predicted, observed)   
    valid_df = rightjoin(predicted, observed, on=[:t_day], makeunique=true)
    A_df = @select(valid_df, :A, :A_1) |> drop_na |> x->@transform(x, :weight=1/nrow(x))
    Q_df = @select(valid_df, :Q, :Q_1) |> drop_na |> x->@transform(x, :weight=1/nrow(x))

    scales = [
                maximum(skipmissing(vcat(A_df.A, A_df.A_1))),
                maximum(skipmissing(vcat(Q_df.Q, Q_df.Q_1)))
            ];

    A_pred = A_df.A   ./ scales[1]
    A_obs  = A_df.A_1 ./ scales[1]
    Q_pred = Q_df.Q   ./ scales[2]
    Q_obs  = Q_df.Q_1 ./ scales[2];

    loss_A = log(wsqeuclidean(A_pred, A_obs, A_df.weight ./ maximum(A_df.weight)))
    loss_Q = log(wsqeuclidean(Q_pred, Q_obs, Q_df.weight ./ maximum(Q_df.weight)))
    total_loss = loss_A + loss_Q

    return total_loss, false
end

"""
Function to estimate AQP parameters from algal drymass and Phosphorus, 
using a Bayesian, likelihood-free approach. <br>
Minimum input required are
- the species name (has to be identical to the name used in reference data)
- path to csv file with reference data
- dictionaries with default global and algal parameters
- dictionary with prior distributions <br>
Output will be stored in a separate folder within `output_dir`, which will carry the species label and date of execution.
$(TYPEDSIGNATURES)
"""

function estimate_AQP_params(
        species_label::String,
        calibration_data_path::String,
        default_global_params::OrderedDict{String,Any},
        default_phyto_params::OrderedDict{String,Any},
        priors::OrderedDict{String,Any};
        output_dir="posteriors",
        k=10,
        n_samples=5000,
        n_samples_init=10000,
        n_samples_validation=100,
        boundaries="none",
        param_labels=OrderedDict(zip(priors.keys, priors.keys)),
        correlation_cutoff=0.3
    )
    
    calibration_data = @pipe CSV.read(calibration_data_path, 
        DataFrame, missingstrings=["missing"], types=Dict("t_day"=>Float64)
        ) |>
    @subset(_, :species.==species_label)
    
    calibration_data = @transform(calibration_data,
        :P=:P_mg_per_L .* default_global_params["volume"],
        :A=:A_mg_per_L .* default_global_params["volume"],
        :Q=:Q_mg_per_L .* default_global_params["volume"],
        :q=:Q_mg_per_L ./ :A_mg_per_L
    )
    
    default_global_params["P0"] = calibration_data.P[1]
    default_global_params["A0"] = [calibration_data.A[1]]
    default_global_params["t_max"] = maximum(calibration_data.t_day)
    
    default_phyto_params["q_min"] = minimum(calibration_data[calibration_data.q.>0,:q])
    default_phyto_params["q_max"] = maximum(calibration_data.q)
    default_phyto_params["sinking_rate"] = 0.
    default_phyto_params["resuspension_factor"] = 0.
    default_phyto_params["m_max"] = 0.
    default_phyto_params["species"] = species_label

    function simulator_function(var_params)
        pdict = copy(default_phyto_params)
        gdict = copy(default_global_params)
        for (par,sampled_value) in zip(priors.keys, var_params)
            pdict[par] = sampled_value
        end
        
        predicted = Run(default_global_params, [pdict], Array{OrderedDict{String,Any}}([]))
        predicted = predicted[:,[
            :t_day, 
            :P, 
            Symbol("A_"*species_label), 
            Symbol("Q_"*species_label)
            ]]
        rename!(predicted, Dict(
                "A_"*species_label=>"A", 
                "Q_"*species_label=>"Q")
        )
        predicted = @transform(predicted, :q = :Q ./ :A)
        predicted.t_day = round.(predicted.t_day, digits=2)
    
        return predicted
    end
    
    @info("#---- Running SMC-ABC ----#")
    t1 = now()
    result = SMC_ABC(
        priors, calibration_data, simulator_function, loss_function_AQP;
        k=k, n_samples=n_samples, n_samples_init=n_samples_init,
        boundaries=boundaries
    )
    t2 = now()
    elapsed_time = t2-t1
    
    valid_df = DataFrame()
    for i in 1:n_samples_validation
        m = simulator_function(posterior_sample(result.accepted))
        m[!,:outer_rep] .= i
        append!(valid_df, m)
    end
    
    resultsdir = output_dir*"\\SMC_ABC_result_"*species_label*"_"*split(string(result.date),"T")[1]
    ispath(resultsdir)==false ? mkdir(resultsdir) : nothing
    ispath(resultsdir*"\\tables")==false ? mkdir(resultsdir*"\\tables") : nothing
    ispath(resultsdir*"\\figures")==false ? mkdir(resultsdir*"\\figures") : nothing
    
    #### tables
    @info("#---- Generating output tables ----#")
    # posterior, priors, trace
    

    result.accepted[!,:q_min] .= default_phyto_params["q_min"]
    result.accepted[!,:q_max] .= default_phyto_params["q_max"]

    CSV.write(resultsdir*"\\tables\\accepted.csv", result.accepted)
    CSV.write(resultsdir*"\\tables\\priors.csv", DataFrame(param=result.priors.keys, dist=result.priors.vals))
    CSV.write(resultsdir*"\\tables\\trace.csv", result.trace)
    
    # posterior prediction
    CSV.write(resultsdir*"\\tables\\posterior_predictions.csv", valid_df)
    
    # point estimates
    
    estimates = DataFrame(
        param = priors.keys,
        estimate =(@subset(result.accepted, :rho.==minimum(:rho)) |> Array)[1:length(priors)],
        CL05 = [quantile(result.accepted[:,i], 0.05) for i in 1:length(priors)],
        CL95 = [quantile(result.accepted[:,i], 0.95) for i in 1:length(priors)]
    )
    CSV.write(resultsdir*"\\tables\\point_estimates.csv", estimates)
    
    # calibration info
   
    calibration_info = DataFrame(
        k = k,
        n_samples = n_samples,
        n_samples_init = n_samples_init,
        boundaries = boundaries,
        acceptance_quantile = result.acceptance_quantile,
        perturbation_factor = result.perturbation_factor,
        date = result.date,
        elapsed_time = elapsed_time,
        final_loss = minimum(result.accepted.rho)
    )
    CSV.write(resultsdir*"\\tables\\calibration_info.csv", calibration_info)
    
    #### figures
    @info("#---- Generating output figures ----#")
    # predicted + observed time-series
    plot_timeseries(valid_df, calibration_data)
    plt.savefig(resultsdir*"\\figures\\posterior_predictions_timeseries.png", dpi=300, bbox_inches="tight")
    #
    ## predicted vs observed scatterplots
    #plot_predicted_vs_observed(valid_df, calibration_data)
    #plt.savefig(resultsdir*"\\figures\\posterior_predictive_check.png", dpi=300, bbox_inches="tight")
    
    # priors vs posteriors
    plot_priors_vs_posteriors(result.accepted, priors; param_labels=param_labels)
    plt.savefig(resultsdir*"\\figures\\priors_and_posteriors.png", dpi=300, bbox_inches="tight")
    
    # correlation heatmap
    plot_correlation_heatmap(result.accepted; cutoff=correlation_cutoff)
    plt.savefig(resultsdir*"\\figures\\posterior_correlations.png", dpi=400, bbox_inches="tight")

    plt.close()
end