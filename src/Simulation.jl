using Base: Real
using StatsBase: normalize!
using Random
using DataFrames, DataFramesMeta
using CSV
using DocStringExtensions
import Base.run
include("DEB.jl")
include("Utils.jl")
include("Inference.jl")
include("Schedules.jl")
include("DataRecorders.jl")

"""
Assign params from one dictionary to another.
$(TYPEDSIGNATURES)
"""
function assign_params!(flea_params::OrderedDict{Symbol,Any}, var_params::OrderedDict{Symbol,Any})
    for (param,value) in zip(var_params.keys, var_params.vals)
        flea_params[param] = value
    end
end

"""
Assign parameters from dictionary to Daphnid instance.
"""
function assign_params!(f::Daphnid, params::OrderedDict{Symbol,Any})
    for (param, value) in zip(params.keys, params.vals)
        setproperty!(f, param, value)
    end
end

"""
Assign parameters from dictionary to Phyto instance.
$(TYPEDSIGNATURES)
"""
function assign_params!(p::Phyto, phyto_params::OrderedDict{Symbol,Any})
    for (param, value) in zip(phyto_params.keys, phyto_params.vals)
        setproperty!(p, param, value)
    end
end

"""
Assign parameters from dictionary to Model instance.
$(TYPEDSIGNATURES)
"""
function assign_params!(m::Model, global_params::OrderedDict{Symbol,Any})
    for (param, value) in zip(global_params.keys, global_params.vals)
        try
            setproperty!(m, param, value)
        catch
            println((param, value))
        end
    end
end

"""
Run the model...
$(TYPEDSIGNATURES)
"""
function run!(m::Model)
    t_max_steps = Int64(m.t_max*m.timestep_phyto)
    for _ in 0:t_max_steps
        step!(m)
    end
end

"""
Run the model with (time-variable) exposure concentrations per time-point passed on as Matrix. 
$(TYPEDSIGNATURES)
"""
function run!(m::Model; C::Union{Nothing,Matrix{Float64}}=nothing)
    if C!=nothing
        if size(C)[2]!=length(m.C)
            error("Number of columns in matrix C has to match length of Model's concentration Vector.")
        end
        t_max_steps = Int64(m.t_max*m.timestep_phyto)
        for i in 0:t_max_steps
            m.C = C[i+1,:] # concentrations are the ith row of Matrix C
            step!(m)
        end
    else
        t_max_steps = Int64(m.t_max*m.timestep_phyto)
        for _ in 0:t_max_steps
            step!(m)
        end
    end
end

"""
Infer survival data from global and daphnid output DataFrames, 
using `global_record` to make sure that all time-points are included.
$(TYPEDSIGNATURES)
"""
function infer_survival(f_out::DataFrame, n_reps, t_max, data_recording_interval)
    if data_recording_interval!=1.
        @warn("Inference of survival will currently fail for data_recording_interval !=1.")
    end
    # count individuals and divide by number of replicates
    survival = combine(
    groupby(f_out, :t_day), 
    x->DataFrame(num_surviving=nrow(x), survival=nrow(x)/n_reps)
    )
    survival.t_day = round.(survival.t_day, sigdigits=2)
    # fill up with zeros
    if maximum(survival.t_day)!=t_max
        missing_days = maximum(survival.t_day)+1:t_max
        append!(survival, DataFrame(
            t_day=missing_days,
            num_surviving=repeat([0], length(missing_days)),
            survival=repeat([0], length(missing_days))
        ))
    end
    return survival
end

"""
Simulate a life-table experiment with `n_reps` equal to the number of simulated replicates.
$(TYPEDSIGNATURES)
"""
function life_table(
    global_params::OrderedDict{Symbol,Any}, 
    phyto_params::Vector{OrderedDict{Symbol,Any}},
    flea_params::OrderedDict{Symbol,Any};
    n_reps::Int64=10
    )
    output = DataFrame()
    output_global_vars = DataFrame()
    let m::Model, f_out::DataFrame
        for rep in 1:n_reps
            m = initialize_model(global_params, phyto_params, [flea_params])
            run!(m)
            f_out = daphnid_output(m, 1)
            f_out[!,:rep] .= rep
            append!(output, f_out)
            append!(output_global_vars, global_output(m))
        end
    end
    repro = output
    growth = @select(output, :t_day, :rep, :L)
    growth[!,:carapace_length] = growth.L ./ flea_params[:shape_factor]
    survival = infer_survival(repro, n_reps, global_params[:t_max], global_params[:data_recording_interval])

    return LifeTableDataset(repro, growth, survival)
end


