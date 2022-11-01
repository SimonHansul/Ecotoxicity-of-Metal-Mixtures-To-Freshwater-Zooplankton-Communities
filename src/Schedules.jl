include("Structures.jl")
include("Environment.jl")
include("AQP.jl")
include("DEB.jl")
include("DataRecorders.jl")
include("Environment.jl")

import Base.step

"""
Check whether current timestep is a timestep in which zooplankton model is called,
based on temporal resolutions for zooplankton and phytoplankton.
$(TYPEDSIGNATURES)
"""
current_step_is_flea_step(m::Model)=(m.t_step%(m.timestep_phyto/m.timestep_fleas))==0

"""
Execute a single phytoplankton model step.
$(TYPEDSIGNATURES)
"""
function step_phyto!(m::Model)
    shuffle!(m.phyto)
    # iterate over phytoplankton populations
    let is_full_day::Bool = time_interval(m, 1.), 
        p::Phyto
        for (i,p) in enumerate(m.phyto)
            # daily events
            if is_full_day
                medium_renewal!(p, m) # medium renewal at certain weekdays
                inoculate!(p,m) # incoluation with algae occurs daily or not at all
            end
            # if phytoplankton is simulated dynamically
            if m.phyto_dynamic
                # calculate limitation terms
                fT_t = fT(m.T, p.T_opt, p.T_min, p.T_max)
                fI_t = fI(m.I, p.I_opt)
                fQ_t = fQ(p.A, p.Q, p.q_min)
                fQP_t = fQP(p.A, p.Q, m.P, p.q_min, p.q_max, p.k_s, m.volume)

                # calculate (mixture) toxicity
                fC_t = prod([fC(conc, p.toxicity_params[i]) for (i,conc) in enumerate(m.C)])
                
                # calculate change in major state variables
                dA_t, dA_sink_t = dA(p.A, fT_t, fI_t, fQ_t, fC_t, p.mu_max, p.m_max, m.D, p.sinking_rate, m.phyto_dynamic) ./ m.timestep_phyto
                dQ_t, dQ_sink_t = dQ(p.A, p.Q, m.P, fQP_t, p.v_max, p.m_max, m.D, p.sinking_rate, m.phyto_dynamic) ./ m.timestep_phyto
                dP_t = dP(p.A, p.Q, m.P, fQP_t, p.v_max, p.m_max) / m.timestep_phyto

                # update phytoplankton state variables
                p.A = p.A + dA_t - dA_sink_t
                p.Q = p.Q + dQ_t - dQ_sink_t
                p.A_sink = p.A_sink + dA_sink_t
                p.Q_sink = p.Q_sink + dQ_sink_t
                
                m.P = max(0, m.P + dP_t)
                m.total_P = max(0, m.P + sum([p.Q for p in m.phyto]))
            else
                dA_t, dA_sink_t = dA(p.A, 1., 1., 1., 1., 0., 0., m.D, p.sinking_rate, m.phyto_dynamic) ./ m.timestep_phyto
                dQ_t, dQ_sink_t = dQ(p.A, 1., 0., 1., 0., 0., m.D, p.sinking_rate, m.phyto_dynamic) ./ m.timestep_phyto

                p.A = p.A - dA_sink_t
                p.Q = p.Q - dQ_sink_t

                p.A_sink = p.A_sink + dA_sink_t
                p.Q_sink = p.Q_sink + dQ_sink_t
            end
            record!(p, m)
        end
    end
end

"""
Execute a single model time-step for fleas.
$(TYPEDSIGNATURES)
"""
function step!(f::Daphnid, m::Model)
    reset!(f)
    f.age += 1/m.timestep_fleas
    derived_attribs!(f)
    resource_choice!(f, m)
    stress!(f, m)
    apply_effects!(f)
    feed!(f, m)
    p_A!(f)
    p_M!(f)
    p_J!(f)
    p_C!(f)
    dE!(f)
    dV!(f)
    dE_H!(f)
    dE_R!(f)
    dq_accel!(f)
    dh_a!(f)
    transition!(f)
    update!(f, m)
    reproduce!(f, m)
    death!(f, m)
    record!(f, m)
end

"""
Execute a single model step.
"""
function step_fleas!(m::Model)
    # randomize order in which individuals are called
    shuffle!(m.fleas)
    V_sum = 0.
    m_sum = 0.

    # reset population sizes
    if time_interval(m, :data_recording_interval)
        for i in 1:length(m.popsizes)
            m.popsizes = setindex(m.popsizes, 0, i)
        end
    end

    for f in m.fleas
        step!(f, m)
        if (f.juvenile)
            V_sum += f.V
            m_sum += f.mass
        end
        # count number of individuals per species
        if time_interval(m, :data_recording_interval)
            m.popsizes = setindex(m.popsizes, m.popsizes[f.species_idx]+1, f.species_idx)
        end
    end
    m.V_density = V_sum / m.volume
    # record population data
    if time_interval(m, :data_recording_interval)
        m.fleas_record_pop = vcat(
            m.fleas_record_pop, 
            hcat(repeat([m.t_day], length(m.popsizes)), 1:length(m.popsizes), m.popsizes)
            )
    end

    # remove dead individuals
    filter!(f->f.die==false, m.fleas)

    return m
end


"""
Execute a model step.
$(TYPEDSIGNATURES)
"""
function step!(m::Model)
    P_inoculation!(m) # P may be inoculated daily
    step_phyto!(m) # execute phytoplankton model step
    # execute zooplankton model step
    # (possibly at lower resolution than phytoplankton model step)
    if current_step_is_flea_step(m)
        step_fleas!(m)
    end
    m.t_day += 1/m.timestep_phyto
    m.t_step += 1
    record!(m)
    return m
end
