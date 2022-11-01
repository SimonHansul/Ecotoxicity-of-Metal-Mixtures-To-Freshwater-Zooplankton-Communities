# Inference.jl
# Functions for parameter inference, 
# Centered around the Sequential Monte Carlo Approximate Bayesian Computation (SMC-ABC) algorithm.
# Simon Hansul
# 2021

using DataFrames, DataStructures
using Distributions
using Distances
using DocStringExtensions
using Pipe
using LsqFit
import StatsBase.sample
using LinearAlgebra
using Dates
using ProgressMeter
using Base.Threads: @threads
import Base.rand

include("Utils.jl")
include("Conversions.jl")

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


"""
Read prior distributions from csv file. For now, function assumes that all distributions are truncated Normal or LogNormal distributions
$(TYPEDSIGNATURES)
"""
function priors_from_file(file::String)
    priors_df = CSV.read(file, DataFrame)
    priors = OrderedDict{Symbol,Any}()
    for row in eachrow(priors_df)
        let dist
            params = eval(Meta.parse(row.params))
            if occursin("{LogNormal{", row.dist)
                dist = Truncated(LogNormal(params...), row.lower, row.upper)
            elseif occursin("{Normal{", row.dist)
                dist = Truncated(Normal(params...), row.lower, row.upper)
            else
                error("Distribution of type "*dist*" not implemented.")
            end
            priors[Symbol(row.param)] = dist 
        end
    end
    return priors
end

"""
Sample from priors, return as dictionary.
$(TYPEDSIGNATURES)
"""
function rand(priors::OrderedDict{Symbol,Any})
    return OrderedDict{Symbol,Any}(
        zip(priors.keys, rand.(priors.vals))
    )
end


infer_inv_transforms(priors::OrderedDict{Symbol,Any}) = Function[occursin("log_", string(k)) ? expinv : identity for k in priors.keys]
infer_inv_transforms(priors::Vector{OrderedDict{Symbol,Any}}) = [infer_inv_transforms(p) for p in priors]

"""
Induce rank correlations to independent samples, using Iman-Conover algorithm. <br>
This function was originally implemented in the MCHammer package. 
MCHammer appears to be insufficiently maintained, so the function was adopted to maintain compatability.
$(TYPEDSIGNATURES)
"""
function iman_conover(ar, cor_mat)
    n_trials = size(ar)[1]
    if typeof(ar) == Array{Float64,2}
          ar = DataFrame(ar)
    end

    #Define how many columns of ISNs are required
    array_dims = size(ar,2)

    #calc cholesky transform to transfer correlation to ISNs
    P = cholesky(cor_mat)

    #Normal() Returns Standard Normals (ISNs)
    R = rand(Normal(),n_trials,array_dims)
    ISN_Matrix = R*P.U
    ISN_Matrix_DF = DataFrame(ISN_Matrix, :auto)

    #apply ranks to create independant correlation rankings matrix
    ISN_Ranked = []
    for i = 1:array_dims
          temp_ranks = ordinalrank(ISN_Matrix_DF[!,i])
          push!(ISN_Ranked, temp_ranks)
    end

    #Reindex the array of samples using the ISN_Ranks. Sort(Array)[OrderingVector]
    final_array=[]
    for i = 1:array_dims
          sorted_array = sort(ar[:,i])[ISN_Ranked[i]]
          push!(final_array, sorted_array)
    end
    return hcat(final_array...)
end

"""
Infer parameter names from prior dictionary.
$(TYPEDSIGNATURES)
"""
function get_par_names(priors::OrderedDict{String,Any})
    return priors.keys
end

"""
Infer parameter names from DataFrame.
"""
function get_par_names(priors::D) where D<:AbstractDataFrame
    colnames = names(priors)
    mask = (colnames.!="distance").&(colnames.!="weights").&(colnames.!="model").&(colnames.!="pmoa")
    par_names = colnames[mask]
    return par_names
end

"""
Draw a posterior sample.
$(TYPEDSIGNATURES)
"""
function posterior_sample(posterior::D; weight_var=:weights) where D<:AbstractDataFrame
    par_names = Symbol.(get_par_names(posterior))
    params = posterior[sample(1:nrow(posterior), Weights(posterior[:,weight_var])), par_names]
    OrderedDict{Symbol,Any}(zip(par_names, params))
end

dictify(par_sample, par_names) = OrderedDict{Symbol,Any}(zip(par_names, par_sample[1:length(par_names)]))

"""
Convert vector of samples into tuple of `DataFrame`s
$(TYPEDSIGNATURES)
"""
function frameify(samples::Vector{Tuple{OrderedDict{Symbol, Any}, OrderedDict{Symbol, Any}}})
    posterior_acc_df = DataFrame(hcat([s[1].vals for s in samples]...)', samples[1][1].keys)
    posterior_prn_df = DataFrame(hcat([s[2].vals for s in samples]...)', samples[1][2].keys)
    return posterior_acc_df, posterior_prn_df
end

"""
Read model selection output dataframes and store in ModelSelectionOutput object.
$(TYPEDSIGNATURES)
"""
mutable struct ModelSelectionOutput
    model_probs::DataFrame
    posteriors::Vector{DataFrame}

    function ModelSelectionOutput(outdir::String; nmodels=4)
        outfiles = glob(outdir*"/SMCModel-outputmodel[1-$nmodels].txt")

        model_probs = CSV.read(outdir*"/SMCModel-outputmodelprobabilities.txt", DataFrame)
        model_probs.Model = Int64.(model_probs.Model)
        model_probs
        acc_dfs = DataFrame[]

        for (i,of) in enumerate(outfiles)
            acc = CSV.read(outfiles[i], DataFrame, skipto=7, header=6)
            push!(acc_dfs, acc)
        end
        new(model_probs, acc_dfs)
    end
end

"""
Posterior sample from model selection output.
$(TYPEDSIGNATURES)
"""
function posterior_sample(msout::ModelSelectionOutput)
    model = sample(out.model_probs.Model, Weights(out.model_probs.Probability))
    prms = posterior_sample(msout.posteriors[model])
    return model, prms
end

offdist = Truncated(Normal(0, 1e-8), 0, 1e-8) # used to turn off effects

"""
Disengage a PMoA by setting all associated parameters to 0.
$(TYPEDSIGNATURES)
"""
function turn_off!(θ::OrderedDict{Symbol,Any}, m)
    for par in θ.keys
        if occursin("_"*m, String(par))
            if typeof(θ[par])<:Distribution
                θ[par] = offdist
            elseif typeof(θ[par])<:Real
                θ[par] = 0.
            else
                "Don't know how to turn off parameter of type "*string(typeof(θ[par]))
            end
        end
    end
end

function turn_off!(df::DataFrame, m)
    df[:,Symbol("S_max_"*m)] .= 0.
end

"""
Check if PMoA is engaged in prior distributions.
$(TYPEDSIGNATURES)
"""
function is_engaged(θ::OrderedDict{Symbol,Any}, PMoA::String)
    if θ[Symbol("S_max_"*PMoA)].upper>1e-8
        return true
    else
        return false
    end
end

function assign_params!(
    params::OrderedDict{Symbol,Any}, 
    parent_sample::OrderedDict{Symbol,Any},
    flea_params::OrderedDict{Symbol,Any}
    )
    assign_params!(flea_params, parent_sample)
    assign_params!(flea_params[:TKTD][STRESSOR_IDX], params)
end

