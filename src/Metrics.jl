# Some common metrics for model evaluation.
# Simon Hansul
# Laboratory of Environmental Toxicology and Chemistry
# Environmental Toxicology unit
# Ghent University
# 2021-09-21

const T = Union{Vector{Union{Missing,R}}, Vector{R}} where R<:Real
include("LifeTableUtils.jl")

"""
Nash-Sutcliffe model efficiency.
$(TYPEDSIGNATURES)
"""
function NSE(predicted, observed)
    if (length(collect(skipmissing(predicted)))>0)&(length(collect(skipmissing(observed)))>0)
        obs_mean = mean(skipmissing(observed))
        numerator = sum(skipmissing((predicted .- observed).^2))
        denominator = sum(skipmissing((observed .- obs_mean).^2))
        return 1 - (numerator/denominator)
    else
        return missing
    end
end

"""
Normalized Nash-Sutcliffe model efficiency.
$(TYPEDSIGNATURES)
"""
function nNSE(predicted::T, observed::T)
    return 1/(2-NSE(predicted, observed))
end

"""
Mean Absolute Error.
$(TYPEDSIGNATURES)
"""
function MAE(predicted::Any, observed::Any)
    return mean(skipmissing(abs.(predicted .- observed)))
end

"""
Standard Deviation of the Absolute Error. 
$(TYPEDSIGNATURES)
"""
function SDAE(predicted::Any, observed::Any)
    return std(skipmissing(abs.(predicted .- observed)))
end

"""
Mean Relative Error.
$(TYPEDSIGNATURES)
"""
function MRE(predicted::Any, observed::Any)
    return mean(skipmissing(abs.(predicted .- observed)./observed))
end

function SDRE(predicted::Any, observed::Any)
    return std(skipmissing(abs.(predicted .- observed)./observed))
end

"""
Mean Absolute Percent Error.
$(TYPEDSIGNATURES)
"""
function MAPE(predicted::Any, observed::Any)
    msk = isfinite.(predicted) .& isfinite.(observed) .& (ismissing.(predicted).==false) .& (ismissing.(observed).==false)
    prd = predicted[msk] .+ 1e-10
    obs = observed[msk] .+ 1e-10
    return (100/(length(obs))) * mean(skipmissing(abs.((obs .- prd)./obs)))
end

"""
Posterior predictive check.
$(TYPEDSIGNATURES)
"""
function PPC(lims, obs)
    lims = hcat(lims...)
    pccount = []

    for (lim,ob) in zip(eachrow(lims), obs)
        if lim[1]<=ob<=lim[2]
            push!(pccount, true)
        else
            push!(pccount, false)
        end
    end

    return sum(pccount)/length(pccount)
end

"""
Posterior predictive check of predicted life-table data. 
Specify part of the life-table data (growth, repro, survival) with `dataset`. 
Specify table columns with lower and upper limits of predictions and observations with 
`col_lo`, `col_hi`, `col_obs`. 
$(TYPEDSIGNATURES)
"""
function PPC(
    predicted::LifeTableDataset, 
    observed::LifeTableDataset, 
    dataset::Symbol, 
    col_lo::Symbol, 
    col_hi::Symbol, 
    col_obs::Symbol
    )
    return @pipe predicted |>
    aggregate(_, :metal) |>
    getfield(_, dataset) |> 
    leftjoin(_, getfield(observed, dataset), on=[:t_day, :metal], makeunique=true) |> 
    drop_na |>
    PPC([_[:,col_lo], _[:,col_hi]], _[:,col_obs])
end

"""
Bayesian Information Criterion. 
Requires number of parameters k, number of data points n, 
maximum likelihood pDH.
"""
function BIC(k::Int64, n::Int64, pDH::Float64)
    return (k*log(n))-(2log(pDH))
end

"""
Akaike Information Criterion. 
Requires number of parameters k and  maximum likelihood pDH.
"""
function AIC(k::Int64, pDH::Float64)
    return 2k - (2log(pDH))
end