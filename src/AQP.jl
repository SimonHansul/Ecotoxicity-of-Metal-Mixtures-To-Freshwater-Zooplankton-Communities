# AQP.jl
# Central data structures and model equations of the AQP model adopted from Weber et al.
# Author: Simon Hansul
# Date (first version): 2020-09-23

using DocStringExtensions, DataStructures
include("Structures.jl")

"""
Light limitation.
$(TYPEDSIGNATURES)
"""
function fI(I::Float64, I_opt::Float64)
    return (I/I_opt)*exp(1-(I/I_opt))
end

"""
Temperature limitation.
$(TYPEDSIGNATURES)
"""
function fT(T::Float64, T_opt::Float64, T_min::Float64, T_max::Float64)
    let T_x
        if T < T_opt
            T_x = T_min
        elseif T >= T_opt
            T_x = T_max
        end
        return exp(-2.3*(((T-T_opt)/(T_x-T_opt)))^2)
    end
end

"""
Limitation by P quota.
$(TYPEDSIGNATURES)
"""
function fQ(A::Float64, Q::Float64, q_min::Float64)
    return 1 - exp(-0.693*((Q/(q_min*A))-1))
end

"""
Limitation of P uptake by dissolved P and P Quota.
$(TYPEDSIGNATURES)
"""
function fQP(A::Float64, Q::Float64, P::Float64, q_min::Float64, q_max::Float64, k_s::Float64, V::Float64)
    #f = (P/V)/(k_s+(P/V))
    #if f>1
    #    error("Scaled functional response cannot be > 1, is $f.")
    #end
    return (((q_max*A)-Q)/((q_max-q_min)*A))*((P/V)/(k_s+(P/V)))
end

"""
Chemical stress.
$(TYPEDSIGNATURES)
"""
function fC(C::Float64, params::Array{Float64,1})
    return 1/(1+(C/params[1])^params[2])
end

"""
Change in phytoplankton population density.
$(TYPEDSIGNATURES)
"""
function dA(A::Float64, fT::Float64, fI::Float64, fQ::Float64, fC::Float64, mu_max::Float64, m_max::Float64, D::Float64, sinking_rate::Float64, phyto_dynamic::Bool)
    if phyto_dynamic
        return ((mu_max*fT*fI*fQ*fC)-m_max-D)*A, sinking_rate*A
    else
        return 0., sinking_rate*A
    end
end

"""
Change in algal P.
$(TYPEDSIGNATURES)
"""
function dQ(
    A::Float64, 
    Q::Float64, 
    P::Float64, 
    fQP::Float64, 
    v_max::Float64,
    m_max::Float64, 
    D::Float64, 
    sinking_rate::Float64, 
    phyto_dynamic::Bool
    )
    if phyto_dynamic
        return (v_max*fQP*A) - (m_max*Q) - (D*Q), sinking_rate*Q
    else
        return 0., sinking_rate*Q
    end
end

"""
Change in dissolved P due to uptake and release by algae.
Inflow/outflow of the flow through system is done in `update_environment()`.
$(TYPEDSIGNATURES)
"""
function dP(A::Float64, Q::Float64, P::Float64, fQP::Float64, v_max::Float64, m_max::Float64)
    return (m_max*Q)-(v_max*fQP*A)
end

"""
Change in stressor concentration.
(TYPEDSIGNATURES)
"""
function dC(C_0::Array{Float64,1}, C::Array{Float64,1}, k::Array{Float64,1}, D::Float64)
    return @. (C_0*D)-(k*C)-(D*C)
end

"""
Inflow/outflow of P in flowthrough system.
$(TYPEDSIGNATURES)
"""
function flowthrough(m::Model)
    # inflow of P via fresh medium
    inflow = (m.D*m.P0) / m.timestep_phyto
    # outflow of P
    outflow = (m.D*m.P) / m.timestep_phyto
    m.P = max(0, m.P + (inflow-outflow))
    return m
end