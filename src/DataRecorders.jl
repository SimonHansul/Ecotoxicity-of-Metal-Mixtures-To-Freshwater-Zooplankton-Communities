using DataFrames, DataFramesMeta
include("Structures.jl")

"""
Check if end of time interval is reached. `interval` is a property of `Model`.
(TYPEDSIGNATURES)
"""
function time_interval(m::Model, interval::Symbol)
    if typeof(getfield(m, interval))<:Real
        let interval_step=getfield(m, interval)*m.timestep_phyto
            return m.t_step%interval_step==0
        end
    elseif getfield(m, interval)=="renewals"
        return day_of_renewal(m)
    else
        error("Invalid value for data recording interval, expect scalar or \"renewals\".")
    end
end

"""
Check if end of time interval is reached. `interval` is a time interval in days.
(TYPEDSIGNATURES)
"""
function time_interval(m::Model, interval::Float64)
    let interval_step=interval*m.timestep_phyto
        return m.t_step%interval_step==0
    end
end

function getproperties(f::Daphnid)
    let props=Vector{Any}(zeros(length(RECORDED_FLEA_ATTRIBS)))
        for (i,attrib) in enumerate(RECORDED_FLEA_ATTRIBS)
            try
                props[i]=getfield(f, attrib)
            catch
                println(attrib)
                error()
            end
        end
        return props
    end    
    #[getfield((f, attrib) for attrib in RECORDED_FLEA_ATTRIBS]
end

"""
Record zooplankton individual state variables.
$(TYPEDSIGNATURES)
"""
function record!(f::Daphnid, m::Model)
    if (m.record_flea_ind)&&(f.juvenile)&&(time_interval(m, :data_recording_interval))
        v = vcat(
            [m.t_day], 
            getproperties(f)
            ) |>
            x->vcat(x, [typeof(tktd.TK)==Vector{TKModel} ?  [tk.C_int for tk in tktd.TK] : tktd.TK.C_int for tktd in f.TKTD])
        push!(m.fleas_record_ind, v)
    end
end

"""
Record phytoplankton state variables.
$(TYPEDSIGNATURES)
"""
function record!(p::Phyto, m::Model)
    if time_interval(m, :data_recording_interval)
        v = [m.t_day, p.species_idx, p.A, p.Q]
        push!(m.phyto_record, v)
    end
end

"""
Record global state variables.
$(TYPEDSIGNATURES)
"""
function record!(m::Model)
    if time_interval(m, :data_recording_interval)
        v = vcat([m.t_day, m.P, m.T], m.C)
        push!(m.global_record, v)
    end
end

"""
Retrieve daphnid output as DataFrame from model object.
$(TYPEDSIGNATURES)
"""
function daphnid_output(m::Model, include_inds::Union{Float64,Int64})
    f_out = DataFrame(hcat(m.fleas_record_ind...)', :auto)

    if length(names(f_out))>0
        rename!(f_out, vcat(
            [:t_day], 
            RECORDED_FLEA_ATTRIBS, 
            [Symbol("C_int_"*string(i)) for i in 1:length(m.C)]
            ))
        f_out = @subset(f_out, :unique_id.<=include_inds)
        return f_out
    else
        nms = vcat([:t_day], RECORDED_FLEA_ATTRIBS, [Symbol("C_int_"*string(i)) for i in 1:length(m.C)])
        f_out = DataFrame(hcat(repeat([[]], length(nms))...), :auto)
        rename!(f_out, nms)
        return f_out
    end

    DataFrame(hcat(m.fleas_record_ind...)', :auto)
end

"""
Retrieve global output as DataFrame from model object.
$(TYPEDSIGNATURES)
"""
function global_output(m::Model)
    g_out = DataFrame(hcat(m.global_record...)', :auto)
    rename!(g_out, vcat([:t_day, :P, :T], [Symbol("C_"*string(i)) for i in 1:length(m.C)]))
    return g_out
end
