# Data structures and functions for Dose-response analysis
# Simon Hansul
# 2022-02-07

using DataStructures
using DocStringExtensions
import Base:copy
import StatsBase:response
import Base:getfield

mutable struct DRCModel
    model::Function
    params::NTuple
end
Base.copy(drc::DRCModel) = DRCModel([copy(getfield(drc, field)) for field in fieldnames(DRCModel)]...)

function response(model::DRCModel, concs::Vector{Float64})
    return model.model(concs, model.params)
end

function response(model::DRCModel, conc::Float64)
    return model.model(conc, model.params)
end

"""
Two-parameter log-logistic function.
$(TYPEDSIGNATURES)
"""
function LL2(x::Float64, p::NTuple{2,Float64})
    return 1/(1+(x/p[1])^p[2])
end


"""
Broadcasted version of LL2. 
$(TYPEDSIGNATURES)
"""
function LL2(x::Vector{Float64}, p::NTuple{2,Float64})
    return @. 1/(1+(x/p[1])^p[2])

end

function LL2inv(y::Float64, p::NTuple{2,Float64})
    return p[1]*(((1/y)-1)^(1/p[2]))
end

#f(x) = c + (d-c) \exp(-\exp(b(\log(x)-\log(e))))

"""
Two-parameter Weibull function.
$(TYPEDSIGNATURES)
"""
function WB2(x::Float64, p::NTuple{2,Float64})
    return exp(-exp(p[2]*(log(x)-log(p[1]))))
end


function WB2(x::Vector{Float64}, p::NTuple{2,Float64})
    return [WB2(xi, p) for xi in x]
end

"""
Bi-phasic log-logistic function, where each phase is described by a two-parameter log-logistic function. 
Requires additional "breakpoint" parameter. <br>
- p1 = EC50_1
- p2 = beta_1
- p3 = EC50_2
- p4 = beta_2
- p5 = breakpoint = relative response at which second phase starts
$(TYPEDSIGNATURES)
"""
function LLBP5(x::Float64, p::NTuple{5,Float64})
    y1 = p[5]/(1+(x/p[1])^p[2]) # first-phase response
    y2 = p[5]/(1+(x/p[3])^p[4]) # second-phase response
    y = y1+y2 # total response
    return y
end

function LLBP5(x::Vector{Float64}, p::NTuple{5,Float64})
    return [LLBP5(xi, p) for xi in x]
end


"""
Asymmetric log-logistic function with additional slope parameter.
- p1 = EC50
- p2 = beta
- p3 = beta_2
$(TYPEDSIGNATURES)
"""
function LLAS3(x::Float64, p::NTuple{3,Float64})
    return 1/((1+((x/p[1])^p[2]))^p[3])
end

function LLAS3(x::Vector{Float64}, p::NTuple{3,Float64})
    return [LLAS3(xi, p) for xi in x]
end

function LL3(x::Float64, p::NTuple{3,Float64})
    return p[3]/(1+(x/p[1])^p[2])
end


"""
3-parameter log-logistic function for application in GUTS model. <br>
Same as LL3, but with negative slope.
$(TYPEDSIGNATURES)
"""
function LL3GUTS(x::Float64, p::NTuple{3,Float64})
    return p[3]/(1+(x/p[1])^-p[2])
end  

function LL3GUTS(x::Vector{Float64}, p::NTuple{3,Float64})
    return @. p[3]/(1+(x/p[1])^-p[2])
end  

function LL3(x::Vector{Float64}, p::NTuple{3,Float64})
    return @. p[3]/(1+(x/p[1])^p[2])
end

"""
Alternative three-parameter log-logistic function used for modelling of increasing effects (costs for structure and maintenance costs). <br>
Corresponding to a four-parameter log-logistic function with negative slope and lower limit fixed to 1. <br>
Set S_max<0 for decreasing response with increasing concentration.
$(TYPEDSIGNATURES).
"""
function LL3GM(x::Float64, p::NTuple{3,Float64})
    return 1 + ((p[3])/(1+(x/p[1])^-p[2]))
end

function LL3GM(x::Vector{Float64}, p::NTuple{3,Float64})
    return [LL3GM(xi, p) for xi in x]
end

function LL3AR(x::Float64, p::NTuple{3,Float64})
    return 1 - ((p[3])/(1+(x/p[1])^-p[2]))
end

function LL3AR(x::Vector{Float64}, p::NTuple{3,Float64})
    return [LL3AR(xi, p) for xi in x]
end

LL4(x, p) = @. p[4] + (p[3]/(1+(x/p[1])^p[2]))

"""
Cedergreend-Ritz-Streibig model. <br>
Parameters are <br>
- p[1] = α = rate of hormetic increase
- p[2] = b = quasi-slope
- p[3] = c = lower limit
- p[4] = d = response in the control
- p[5] = e = inflection point
- p[6] = f = hormesis parameter
"""
CRS6(x::Float64, p::NTuple{6,Float64}) = @.(p[3]+((p[4]-p[3]+p[6]*exp(-1/(x^p[1])))/(1+exp(p[2]*(log(x)-log(p[5]))))))

"""
4-parameter CRS model. Lower limit and response in the control are fixed to 0 and 1, respectively.
"""
CRS4(x::Float64, p::NTuple{4,Float64}) = @.(((1+p[4]*exp(-1/(x^p[1])))/(1+exp(p[2]*(log(x)-log(p[3]))))))

"""
4-parameter CRS model transformed to u-shape, where the response in the control is 1, the lower limit is 0 and the maximum is 2.0.
"""
function CRS4U(x::Float64, p::NTuple{4,Float64})
    y = @. 2.0 .- (((1+p[4]*exp(-1/(x^p[1])))/(1+exp(p[2]*(log(x)-log(p[3]))))))
    return max.(0, y)
end

"""
6-parameter U-shaped CRS model.
"""
function CRS6U(x::Float64, p::NTuple{6,Float64})
    y = @. p[4]-(p[3]+((p[4]-p[3]+p[6]*exp(-1/(x^p[1])))/(1+exp(p[2]*(log(x)-log(p[5]))))))
end

"""
Re-scaled 6-parameter U-shaped CRS model.
"""
function CRS6US(x::Float64, p::NTuple{6,Float64})
    y = @. 1 + (p[4]-(p[3]+((p[4]-p[3]+p[6]*exp(-1/(x^p[1])))/(1+exp(p[2]*(log(x)-log(p[5])))))))
    return max.(0, y)
end

"""
Re-scaled U-shaped CRS model with parameter C fixed to 0. <br>
Parameters are <br>
- alpha: rate of hormetic increase
- b: slope of the inclining part of the curve
- d: maximum stress 
- e: inflection point of the inclining part of the curve
- f: hormesis parameter
"""
function CRS5US(x::Float64, p::NTuple{5,Float64})
    y = 1 + (p[3]-(((p[3]+p[5]*exp(-1/(x^p[1])))/(1+exp(p[2]*(log(x)-log(p[4])))))))
    return max.(0, y)
end

"""
$(TYPEDSIGNATURES)
"""
function CRS5US(x::Vector{Float64}, p::NTuple{5,Float64})
    y = @. 1 + (p[3]-(((p[3]+p[5]*exp(-1/(x^p[1])))/(1+exp(p[2]*(log(x)-log(p[4])))))))
    return max.(0, y)
end

"""
Hockey-stick model.
$(TYPEDSIGNATURES)
"""
HS(x::Vector{Float64}, p::NTuple{2,Float64}) = @. (1/p[1]) * max(0, x-p[2])


"""
Define prior distributions for various dose-response models. 
Assumes that `data` is a DataFrame with two columns called x and y, 
corresponding to exposure concentration and relative response, respectively.
$(TYPEDSIGNATURES)
"""
function set_priors(model_label::String, data::D) where D<:DataFrame
    let edist, bdist
        # if effect >= 50% as been observed
        if maximum(data.y)>=0.5
            # test concentration closest to EC50 = prior mean of EC50
            mean_response = combine(groupby(data, :x), :y=>mean)
            e = mean_response[:,1][argmin(abs.(mean_response[:,2].-0.5))]
            
            # β as Uniform distribution
        # otherwise, take 2-times the highest test concentration as mean
        else
            e = maximum(data.x)*2
        end
        # prior distribution for EC50 as log-Normal with fixed sd on a log-scale
        edist = LogNormal(log(e), 1)
        bdist = Truncated(Normal(2, 20), 1, Inf)
        # log-logistic and weibull model
        if model_label=="LL2"
            return [edist, bdist]
        elseif model_label=="WB2"
            return [edist, bdist]
        # biphasic model: priors for both inflection points at 25% at 75% effect
        elseif model_label=="LLBP5"
            e1 = data[:,1][argmin(abs.(data[:,2].-0.25))]
            e2 = data[:,1][argmin(abs.(data[:,2].-0.75))]
            return [
                LogNormal(log(e1), 1),
                bdist, 
                LogNormal(log(e2), 1),
                bdist,
                Uniform(0.1,0.9)   
            ]
        # asymmetric: same as for log-logistic, but two slope parameters
        elseif model_label=="LLAS3"
            return [
                edist,
                bdist,
                bdist
            ]
        # cedergreen-ritz streibig: same as for log-logistic, 
        # parameter α between 0.25 and 1.0 (cf. original Cedergreen-Ritz-Streibig publication),
        # parameter f between 0 and 1.5-times the maximum observed hormetic effect
        elseif model_label=="CRS4"
            return [
                Uniform(0.25, 1.0),
                bdist,
                edist,
                Uniform(0, (maximum(data[:,2])-1)*1.5)
            ]
        else
            error("Model $model_label not implemented.")
        end
    end
end

#### We may want the getfield/setfield functions to be used by default
#### Traceur however seems to disagree, but this option is still slightly faster for now.
"""
function Base.getproperty(drc::DRCModel, f::Symbol)
    return getfield(drc, f)
end
"""
