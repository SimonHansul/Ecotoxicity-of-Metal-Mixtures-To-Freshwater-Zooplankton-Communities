# Utils.jl
# All kinds of convenvience functions.
# Simon Hansul
# 2021-02-21

using DataFrames, Pipe
using LaTeXStrings
using DocStringExtensions
using Unitful
import Base:copy

# weight-length relationship from urabe & watanabe (1991)
weight_from_length(length_cm) = ismissing(length_cm) ? missing : (2*1.6e-3u"mg/mm^3"*((length_cm*10)^3)*u"mm^3").val


import Distributions.scale


"""
Coefficient of variation.
$(TYPEDSIGNATURES)
"""
function cv(x::Vector)
    return std(skipmissing(x))/mean(skipmissing(x))
end

"""
Scale values in Dictionary.
$(TYPEDSIGNATURES)
"""
function scale(x::D) where D<:AbstractDict
    x.vals = x.vals ./ sum(x.vals)
end

"""
Scales values in Vector.
$(TYPEDSIGNATURES)
"""
function scale(x::V) where V<:AbstractVector
    x = x ./ sum(x)
end


function drop_na(df::AbstractDataFrame; verbose=false)
    n0 = nrow(df)
    df2 = df[completecases(df),:]
    dn = nrow(df2)-n0
    if verbose
        @info("Dropped $dn of $n0 rows containing missing values.")
    end
    return df2
end

function drop_na!(df::DataFrame)
    df = drop_na(df)
end

"""
Replace missing values in DataFrame with zeros.
$(TYPEDSIGNATURES)
"""
function replace_na(
    df::AbstractDataFrame;
    cols::Array{Symbol,1}=[:carapace_length],
    replace_val=0.0
    )
    for col in cols
        df[ismissing.(df[:,col]),col] .= replace_val
    end
    return df
end

"""
Skip inifite values.
$(TYPEDSIGNATURES)
"""
skipinf(x::Array) = x[isfinite.(x)]

"""
Return `missing` if none of the data is finite or non-nothing or non-missing. Otherwise, return mean.
$(TYPEDSIGNATURES)
"""
function robust_mean(x)
    x = filter(!isnothing, x) |> Array{Any,1}
    x = filter(!ismissing, x) |> Array{Any,1}
    x = filter(isfinite, x) |> Array{Any,1} 
    
    if length(x)==0
        return missing
    else
        return mean(x)
    end
end

function robust_sum(x)
    x = filter(!isnothing, x) |> Array{Any,1}
    x = filter(!ismissing, x) |> Array{Any,1}
    x = filter(isfinite, x) |> Array{Any,1} 
    
    if length(x)==0
        return missing
    else
        return sum(x)
    end
end

"""
Calculate sum, correcting for initial length if elements are nothing, missing or infinite.
$(TYPEDSIGNATURES)
"""
function length_corrected_sum(x)
    num_elems = length(x)
    x = filter(!isnothing, x) |> Array{Any,1}
    x = filter(!ismissing, x) |> Array{Any,1}
    x = filter(isfinite, x) |> Array{Any,1} 
    
    if length(x)==0
        return missing
    else
        return sum(x) * (num_elems/length(x))
    end
end

"""
Infer treatment types (categorical), levels (ordinal) and names (type+level) from Array of exposure concentrations.
$(TYPEDSIGNATURES)
"""
function get_treatment_names(exposure::Array{Array{Float64,1},1}, stressors::Array{Symbol,1})
    treatment_type = ["co"]
    treatment_level = [0]
    treatment = ["co"]

    treatment_level_counter = 0
    stressor_pre = "co"
    for (i,x) in enumerate(exposure[2:end])
        sum(x.>0)==1 ? stressor = string(stressors[x.>0][1]) : stressor = "mix"
        if stressor==stressor_pre
            treatment_level_counter += 1
        else
            treatment_level_counter = 1
        end
        push!(treatment_type, string(stressor))
        push!(treatment_level, treatment_level_counter)
        push!(treatment, stressor*string(treatment_level_counter))
        stressor_pre = stressor
    end
    return treatment_type, treatment_level, treatment
end

"""
Geometric series created from range of values within a vector.
$(TYPEDSIGNATURES)
"""
geomrange(v::Vector{Float64}; length=50) = 10 .^ range(log10(minimum(v)), log10(maximum(v)), length=length)

"""
Geometric series created from two extreme values.
$(TYPEDSIGNATURES)
"""
geomrange(a::Real, b::Real; length=50) = 10 .^ range(log10(a), log10(b); length=length)

"""
Parse string to Vector of Floats.
"""
vectify(x) = parse.(Float64, split(split(split(x, "[")[end], "]")[1]," "))

"""
Calculate difference along a vector, inserting 0. as first element.
$(TYPEDSIGNATURES)
"""
diffvec(x) = vcat([0.], diff(x))


copy(s::String)=s
copy(f::Function)=f


"""
Function to write DataFrames to disc during loop. Will overwrite existing file if step==1 and append if step>1.
$(TYPEDSIGNATURES)
"""
function wrappend(file::String, data::DataFrame, step::Int64)
    if (isfile(file)==false)&(step>1)
        error("Attempt to append to non-existing file: step>1 but file does not exist.")
    end
    if step==1
        CSV.write(file, data, append=false)
    else
        CSV.write(file, data, append=true)
    end
end

"""
Formatted rounding to significant digits (omitting decimal point when appropriate). 
Returns rounded number as string.
$(TYPEDSIGNATURES)
"""
function fround(x; sigdigits=2)
    xround = string(round(x, sigdigits=sigdigits))
    if xround[end-1:end]==".0"
        xround = string(xround[1:end-2])
    end
    return xround
end

"""
Identify element of `possibilities` occuring in `x`. \\
E.g., which_in("I like apples", ["apples", "bananas"]) returns "apples". \\
If multiple possibilities occur, first detection is returned. 
$(TYPEDSIGNATURES)
"""
function which_in(x::String, possibilities::Vector{String}; none_found_return_val="") 
    idxs = findall(z->occursin(z, x), possibilities)
    if length(idxs)>0
        return possibilities[idxs[1]]
    else
        return none_found_return_val
    end
end







