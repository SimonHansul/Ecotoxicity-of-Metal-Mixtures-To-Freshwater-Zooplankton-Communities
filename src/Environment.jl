# Environment.jl
# Handling of environmental and experimental processes.
# Simon Hansul
# 2021-02-16

"""
Inoculate with new algal biomass.
Algal P is updated based on the size-weighted mean between old P quota and P quota in inoculum.
$(TYPEDSIGNATURES)
"""
function inoculate!(p::Phyto, m::Model)
    if m.phyto_constant
        p.A = m.A_inoculate[p.species_idx]
    else
        let P_quota_pre::Float64, P_quota_post::Float64
            P_quota_pre = p.A>0 ? p.Q / p.A : 0. # record current P quota + biomass
            A_pre = p.A
            p.A += m.A_inoculate[p.species_idx] # update algal biomass
            P_quota = [P_quota_pre, p.q_init] # new P quota are weighted average between P quota in inoculate and old P quota
            P_quota_post = mean(P_quota, Weights([A_pre, m.A_inoculate[p.species_idx]]))
            p.Q = p.A > 0 ? p.A * P_quota_post : 0. # update algal P correspondingly
        end
    end
end


"""
Update environmental state variables. \n
$(SIGNATURES)
"""
function P_inoculation!(m::Model)
    if time_interval(m, 1.)
        m.P += m.P_inoculate
        if true in [floor(m.t_day)%6==x for x in m.renew_medium]
            m.P = m.P * (1-m.frct_medium_renewal)
        end
    end
    return m
end

function day_of_renewal(m::Model)
    floor(m.t_step)%(7*m.timestep_phyto) in m.renew_medium .* m.timestep_phyto
end

"""
Simulate medium renewal. 
Results in a corresponding removal of algal biomass while maintaining P quota. <br>
A fraction of sedimented algal biomass may be re-suspended based on `resuspension_factor`.
$(TYPEDSIGNATURES)
"""
function medium_renewal!(p::Phyto, m::Model)
    if day_of_renewal(m)
        P_quota_pre = p.A > 0 ? p.Q / p.A : 0. # record P quota
        p.A = max(0, p.A * (1-m.frct_medium_renewal)) # remove a fraction of biomass, according to the fraction of medium that is renewed
        p.Q = p.A > 0 ? p.A * P_quota_pre : 0. # update algal P to maintain P quota
        # a fraction (resuspension_factor) of sedimented biomass goes back into suspension
        p.A += p.A_sink * p.resuspension_factor
        p.Q += p.Q_sink * p.resuspension_factor
        p.A_sink = 0.
        p.Q_sink = 0.
    end
end

