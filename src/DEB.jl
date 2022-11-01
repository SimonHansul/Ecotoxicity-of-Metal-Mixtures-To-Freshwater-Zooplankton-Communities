# DEB.jl
# Model equations for the standard DEB model, plus ageing, TKTD and crowding submodels.
# Note: This script only contains function definitions, but does not execute any computations. 
# To follow back how functions are applied in each timestep, see Schedules.jl
# Author: Simon Hansul
# Laboratory of Environmental Toxicology and Aquatic Ecology
# Ghent University
# 2021-10-29

using Random, Distributions
using DocStringExtensions
include("Structures.jl")
include("InitStructures.jl")
include("Constants.jl")

if isdefined(Main, :EMPTVEC)==false
    const EMPTVEC = Float64[]
end

#### Function definitions
"""
Adaptive plasticity.
($TYPEDSIGNATURES)
"""
function adaptive_plasticity(f::Daphnid)
    error("Adaptive plasticity still needs to be implemented.")
    return f
end

""" 
Initialize a juvenile with age <= age _ range _ init. 
Age is given in days and refers to time since hatching 
(other than flea attribute f.age which refers to time since oviposition).
$(SIGNATURES)
"""
function init_random_age!(f::Daphnid, m::Model)
    let init_age_from_hatching = rand(Uniform(0, m.age_init_max)), age_counter=0.
        while age_counter <= init_age_from_hatching
            [p.A = 100. for p in m.phyto] # unlimited food before onset of experiment
            step!(f, m)
            if f.juvenile
                age_counter += 1/m.timestep_fleas
            end
            if (f.age>ABORT_INIT_EMB)&&(f.juvenile==false)
                error("Development of initial individuals lasts too long. Coding error or implausible parameters?")
            end
            f.die = false # initial individuals did, by definition, not die
        end
        # make sure individuals start with correct reserve levels
        f.e = E_SCALED_INIT_MOTHERS
        f.E = (f.e*f.p_A_m_0*f.V)/f.v
        [p.A = m.A_0[i] for (i,p) in enumerate(m.phyto)] # reset food level to correct value
    end
end

"""
Calculation of ingestion rate. Returns ingestion rate `J_X` (J d^-1), current food density `X_p` (J d^-1) and scaled functional response `f` (-).
"""
function ingestion_rate(
    juvenile::Bool, 
    A_p::Float64, 
    volume::Float64, 
    L::Float64, 
    J_X_A_m::Float64, 
    F_m::Float64,
    timestep_fleas::Int64
    )
    # X_p is current food density converted from dry mass to energy
    let J_X, X_p, f, K=J_X_A_m/F_m
        # Embryos don't feed
        if juvenile==false
            f = 0.
            X_p = 0.
            J_X = 0.
        # Juveniles and adults do
        else
            X_p = (A_p/volume)*MG_ALG_TO_J # conversion from mg dry mass to energy L^-1
            f = X_p/(K + X_p) # scaled functional response
            # control for possible mass balance errors
            J_X = min(X_p, (f * L^2 * J_X_A_m))
        end
        return J_X, X_p, f
    end
end

"""
Simulation of feeding process.
Returns updated model-level (food density) and individual-level (conditional food density, functional response) state variables.
$(SIGNATURES)
"""
function feed(;
    juvenile::Bool,
    F_m::Float64, 
    L::Float64,  
    J_X_A_m::Float64,
    timestep_fleas::Int64, 
    A_p::Float64,
    volume::Float64
    )
    J_X, X_p, f = ingestion_rate(juvenile, A_p, volume, L, J_X_A_m, F_m, timestep_fleas)
    if (X_p - J_X/timestep_fleas)<-1e-3
        error("Mass balance error above tolerance occured during feeding.")
    end
    A_p = max(0, A_p - (J_X/MG_ALG_TO_J)/timestep_fleas)
    return J_X, X_p, f, A_p
end

"""
Sample a resource from vector of available resource types. <br>
Sampling weights are proportional resource biomass. 
$(TYPEDSIGNATURES)
"""
function resource_choice!(f::Daphnid, m::Model)
    f.resource_idx = sample(1:length(m.phyto), Weights([p.A for p in m.phyto]))
    f.resource_correction = f.resource_corr_facts[f.resource_idx]
end

"""
$(TYPEDSIGNATURES)
"""
function feed!(
        f::Daphnid, 
        m::Model
    )
    let A_p = m.phyto[f.resource_idx].A, Q_p = m.phyto[f.resource_idx].Q, i=f.resource_idx
        q_t = A_p > 0 ? Q_p / A_p : 1.
        f.q_p_rel = f.juvenile ? q_t : 1.

        f.J_X, f.X_p, f.f, m.phyto[i].A = feed(; 
            juvenile=f.juvenile,
            F_m=f.F_m, 
            L=f.L, 
            J_X_A_m=f.J_X_A_m, 
            timestep_fleas=m.timestep_fleas, 
            A_p=A_p,
            volume=m.volume
            )
        
        # algal P is updated to maintain internal P concentration
        m.phyto[i].Q = m.phyto[i].A > 0 ? m.phyto[i].A * q_t : 0.
    end
end

function init_statevars!(f::Daphnid)
    let ind_var = exp(rand(Normal(0, f.cv)))
        f.J_X_A_m_0 *= ind_var
        f.p_A_m_0 = f.J_X_A_m_0 * f.kappa_EX_max

        f.juvenile = false
        f.adult = false
        f.age = 0.
        f.L = L_O
        f.V = f.L^3
        f.V_m_hist = f.V

        f.f = 0.
        f.X_p = 0.
        f.J_X = 0.
        f.p_M = 0.
        f.p_J = 0.
        f.p_C = 0.
        f.q_accel = 0.
        f.h_a = 0.
        f.clutch_size = 0
        f.E_H = 0.
        f.E_R = 0.
        f.embryos = Daphnid[]
        f.cum_repro = 0.

        f.J_X_A_m = f.J_X_A_m_0
        f.kappa = f.kappa_0
        f.kappa_R = f.kappa_R_0
        f.kappa_EX = f.kappa_EX_max
        f.E_G = f.E_G_0
        f.p_M_V = f.p_M_V_0
        f.E_H_b = f.E_H_b_0

        #f.TKTD.h_ED50 .*= @. exp(rand(Normal(0, f.TKTD.h_cv)))
        f.s_G_z = f.s_M_z = f.s_A_z = f.s_R_z = f.s_E_o_z = f.s_kappa_z = 1.
        f.s_G_cr = f.s_M_cr = f.s_A_cr = f.s_R_cr = 1.

        reset_internal_concentrations!(f.TKTD)
        
        f.q_p_rel = f.resource_correction = 1.
        f.h_z = 0.
        f.cause_of_death = 0.
        f.die = false
        f.stress_crowding = ones(4)
    end
end

"""
Initialization of an embryo. <br>
Initial reserves are calculated based on the scaled reserve density of the mother.
$(TYPEDSIGNATURES)
"""
function initialize_embryo!(f::Daphnid, e_mother::Float64, E_o::Float64)
    init_statevars!(f)
    f.E = E_o
    derived_attribs!(f)
end

"""
Calculation of assimilation rate p_A, based on ingestion rate J_X (J/d) and 
assimilation efficiency kappa_EX (-). For embryos, p_A=0. <br>
This formulation assumes that surface area L^2 is already taken into account in J_X!
$(TYPEDSIGNATURES)
"""
function p_A(
    J_X::Float64,
    kappa_EX::Float64,
    juvenile::Bool
    )
    if juvenile
        return J_X * kappa_EX
    else
        return 0.
    end
end

"""
Calculation of assimilation rate for individual.
$(TYPEDSIGNATURES)
"""
function p_A!(f::Daphnid)
    f.p_A = p_A(f.J_X, f.kappa_EX, f.juvenile)
end

"""
Somatic maintenance p_M (J/d), given structural length L (cm) 
and volume-specific somatic maintenance constant p_M_V
$(TYPEDSIGNATURES)
"""
function p_M(V::Float64, p_M_V::Float64)
    return p_M_V*V
end

"""
$TYPEDSIGNATURES
"""
function p_M!(f::Daphnid)
    f.p_M = p_M(f.V, f.p_M_V)
end

"""
Maturity maintenance (J/d).
$(TYPEDSIGNATURES)
"""
function p_J(k_J::Float64, E_H::Float64)
    return k_J * E_H
end

"""
$(TYPEDSIGNATURES)
"""
function p_J!(f::Daphnid)
    f.p_J = p_J(f.k_J, f.E_H)
end

"""
Mobilization rate (J/d).
$(TYPEDSIGNATURES)
"""
function p_C(
    E::Float64, 
    v::Float64, 
    L::Float64,
    p_M_V::Float64,
    E_G::Float64,
    V::Float64,
    kappa::Float64
    )
    return E*(E_G*(v/L)+p_M_V)/(E_G+((E/V)*kappa))
end

"""
$(TYPEDSIGNATURES)
"""
function p_C!(f::Daphnid)
    f.p_C = p_C(f.E, f.v, f.L, f.p_M_V, f.E_G, f.V, f.kappa)
end

"""
Change in reserves (J/d).
$(TYPEDSIGNATURES)
"""
function dE(p_A::Float64, p_C::Float64, juvenile::Bool)
    if juvenile
        return p_A-p_C
    else
        return -p_C
    end
end

"""
$(TYPEDSIGNATURES)
"""
function dE!(f::Daphnid)
    f.dE = dE(f.p_A, f.p_C, f.juvenile)
end

"""
Change in structural volume (1/d).
$(TYPEDSIGNATURES)
"""
function dV(E_G::Float64, kappa::Float64, p_C::Float64, p_M::Float64)
    return ((kappa*p_C) - p_M)/E_G
end

"""
$(TYPEDSIGNATURES)
"""
function dV!(f::Daphnid)
   f.dV = dV(f.E_G, f.kappa, f.p_C, f.p_M) 
end

"""
Change in maturity.
$(TYPEDSIGNATURES)
"""
function dE_H(kappa::Float64, p_C::Float64, p_J::Float64, E_H::Float64, E_H_p::Float64)
    if E_H < E_H_p
        return ((1-kappa)*p_C) - p_J
    else 
        return 0.
    end
end

"""
$(TYPEDSIGNATURES)
"""
function dE_H!(f::Daphnid)
    f.dE_H = dE_H(f.kappa, f.p_C, f.p_J, f.E_H, f.E_H_p)
end

""" 
$(TYPEDSIGNATURES)
"""
function dE_R(kappa::Float64, p_C::Float64, p_J::Float64, E_H::Float64, E_H_p::Float64)
    if E_H>=E_H_p
        return (1-kappa)*p_C - p_J
    else
        return 0.
    end
end

"""
$(TYPEDSIGNATURES)
"""
function dE_R!(f::Daphnid)
    f.dE_R = dE_R(f.kappa, f.p_C, f.p_J, f.E_H, f.E_H_p)
end


"""
Change in ageing acceleration.
($TYPEDSIGNATURES)
"""
function dq_accel(
    q_accel::Float64,
    V::Float64,
    V_m::Float64,
    s_G::Float64,
    h_a_accel::Float64,
    e::Float64,
    v::Float64,
    L::Float64,
    dL::Float64,
    juvenile::Bool
)
    if juvenile
        r = (3/L)*dL
        return ((q_accel*(V/V_m)*s_G+h_a_accel))*e*(((v/L)-r)-r*q_accel)
    else
        return 0.
    end
end

"""
$(TYPEDSIGNATURES)
"""
function dq_accel!(f::Daphnid)
    f.dq_accel = dq_accel(f.q_accel, f.V, f.V_m, f.s_G, f.h_a_accel, f.e, f.v, f.L, f.dL, f.juvenile)
end

"""
Change in ageing hazard rate.
$(TYPEDSIGNATURES)
"""
function dh_a(q_accel::Float64, L::Float64, dL::Float64, h_a::Float64)
    r = (3/L)*dL
    return q_accel - r*h_a
end

"""
$(TYPEDSIGNATURES)
"""
function dh_a!(f::Daphnid)
    f.dh_a = dh_a(f.q_accel, f.L, f.dL, f.h_a)
end


"""
Manage transition between life stages.
$(TYPEDSIGNATURES)
"""
function transition!(f::Daphnid)
    if f.E_H >= f.E_H_b
        f.juvenile = true
    end
    if f.E_H > f.E_H_p
        f.adult = true
    end
end

"""
Calculate compound attributes, such as maturity rate constant, 
scaled reserve density, etc. <br>
$(TYPEDSIGNATURES)
"""
function compound_attribs!(f::Daphnid)
    f.k_M = f.p_M_V / f.E_G
    f.k_J_0 = f.p_M_V_0 / f.E_G_0
    f.g = (f.E_G*f.v) / (f.kappa*f.p_A_m_0)
    f.U_H_b = f.E_H_b / f.p_A_m_0
    f.U_H_p = f.E_H_p / f.p_A_m_0 
    f.p_A_m = f.J_X_A_m * f.kappa_EX_max
    f.K = f.J_X_A_m / f.F_m
    f.e = (f.E/f.V)*(f.v/f.p_A_m_0)
    f.L_m = f.kappa*(f.p_A_m_0/f.p_M_V_0)
    f.V_m = f.L_m^3
    f.mass = gellermass(f.L/f.shape_factor) #(f.m_v*f.V) + ((f.E+f.E_R)*RESERVE_TO_MASS)

end

"""
"Particulate incomplete beta function" as reported by Kooijman et al. 
Used to calculate initial reserves.
$(TYPEDSIGNATURES)
"""
kooijmans_beta(x) = sqrt(3)*(atan((1+(2x)^(1/3))/sqrt(3))-atan(1/sqrt(3)))+(0.5*log10(1+x^(1/3))+x^(2/3))-log10(1-x^(1/3))

"""
Calculation of initial reserves E_o. Strictly speaking, this is only valid if k_J==k_M. <br>
However, we assume that toxicant stress does not affect investment per embryo. <br> 
Compared to Kooijman et al. (2008), `f` was replaced by `e` since food density is not assumed constant.
$(TYPEDSIGNATURES)
"""
function initial_reserves(L_b::Float64, E_G_0::Float64, kappa_0::Float64, e::Float64, v::Float64, k_M::Float64, p_A_m_0::Float64)
    let g_0 = (E_G_0*v)/(kappa_0*p_A_m_0)
        β = kooijmans_beta(g_0/(g_0+e))
        term_1 = (1/(L_b*(g_0+e)^(1/3)))
        term_2 = β/((3g_0^(1/3)v)/(k_M))
        U_E_o = ((term_1-term_2)^-3)/v
        return max(0, U_E_o * p_A_m_0)
    end
end

"""
Update initial reserves/energy investment per egg.
$(TYPEDSIGNATURES)
"""
function initial_reserves!(f::Daphnid)
    f.E_o = initial_reserves(f.L_b, f.E_G_0, f.kappa_0, f.e, f.v, f.p_M_V_0 / f.E_G_0, f.p_A_m_0)*f.s_E_o_z
end

"""
Update embryonal parametes. <br>
$(TYPEDSIGNATURES)
"""
function emb_attribs!(f::Daphnid)
    let g_0 = (f.E_G_0*f.v)/(f.kappa_0*f.p_A_m_0),
        U_H_b_0=f.E_H_b_0/f.p_A_m_0

        f.L_o = L_O # initial length is fixed
        f.V_o = V_O # initial volume is fixed
        f.E_H_o = 0. # initial maturity is 0E
        f.L_b = ((U_H_b_0*f.v)/(g_0*(1-f.kappa_0)))^(1/3)
        f.V_b = f.L_b^3
        initial_reserves!(f)
    end
end

"""
Calculate derived attributes, being compound and embryonal attributes.
$(TYPEDSIGNATURES)
"""
function derived_attribs!(f::Daphnid)
    compound_attribs!(f)
    emb_attribs!(f)
end

"""
Update state variables.
$(TYPEDSIGNATURES)
"""
function update!(f::Daphnid, m::Model)
    let L_t=f.L
        f.E = max(0., f.E+(f.dE/m.timestep_fleas))
        f.V = max(0., f.V+(f.dV/m.timestep_fleas))
        f.V_m_hist = max(f.V_m_hist, f.V)
        f.L = f.V^(1/3)
        f.dL = (f.L-L_t)*m.timestep_fleas
        f.E_H = max(0., f.E_H+(f.dE_H/m.timestep_fleas))
        f.E_R = max(0., f.E_R+(f.dE_R/m.timestep_fleas))
        f.q_accel = max(0., f.q_accel+(f.dq_accel/m.timestep_fleas))
        f.h_a = max(0., f.h_a+(f.dh_a/m.timestep_fleas))
        transition!(f)
    end
    derived_attribs!(f)
end

"""
Inheritance of attributes.
$(SIGNATURES)
"""
function inherit!(embryo::Daphnid, mother::Daphnid)
    # iterate over attributes
    for attrib in INHERITED_ATTRIBUTES
        # transfer value from mother to embryo
        new_value = deepcopy(getfield(mother, attrib))
        setfield!(embryo, attrib, new_value)
    end
end

"""
Simulate embryonal development.
$(TYPEDSIGNATURES)
"""
function embryonal_development!(f::Daphnid, m::Model)
    # keep track of which embryos will transition to juvenile
    let unhatched_emb = BitVector(ones(length(f.embryos)))
        # for every embryo
        for (i,e) in enumerate(f.embryos)
            # execute basic schedule
            step!(e, m)
            # if maturity at birth is reached
            if e.E_H>=e.E_H_b
                unhatched_emb[i] = 0 # embryo will be removed from brood pouch
                e.E_H_b = f.E_H_b_0
                push!(m.fleas, e) # individual is added to population
                # reproduction is observed at the time of transition to juvenile
                m.cum_repro += 1
                f.cum_repro += 1
            end
        end
        # keep embryos in development
        f.embryos = f.embryos[unhatched_emb]
        f.clutch_size = length(f.embryos)
    end
    f.embryos = filter(e->e.die==false, f.embryos)
end

"""
Reproduction during molting.
$(TYPEDSIGNATURES)
"""
function reproduce!(f::Daphnid, m::Model)
    let E_R_spent::Float64, 
        num_eggs::Float64
        # molting controls timing of reproduction
        f.molt_time += (1/m.timestep_fleas)
        if f.adult&(f.molt_time>f.tau_molt)
            f.molt_time = 0. # reset time since last molt
            # convert reproduction buffer to number of eggs, 
            # using reserve at oviposition
            num_eggs = f.E_o>0 ? (f.kappa_R*f.E_R)/f.E_o : 0.
            # convert to integer
            if isnan(num_eggs)
                f.clutch_size = 0
            else
                f.clutch_size = Int64(floor(num_eggs))
            end
            # back-calculate calculate energy spent on reproduction
            E_R_spent = (f.clutch_size * f.E_o)/f.kappa_R
            # throw error if mass balance does not add up...
            if  max(0, f.E_R-E_R_spent)<-1e-3
                error("Mass balance error above tolerance in `reproduce!(f,m)`")
            else
                f.E_R = max(0, f.E_R-E_R_spent)
            end
            # create offspring individuals as embryos
            for _ in 1:f.clutch_size
                embryo = Daphnid()
                m.unique_id_count += 1
                embryo.unique_id = m.unique_id_count
                inherit!(embryo, f)
                initialize_embryo!(embryo, f.e, f.E_o)
                embryo.E_H_b = embryo.E_H_b_0 * f.s_E_o_z
                push!(f.embryos, embryo)
            end
        else
            f.clutch_size = 0
        end
        embryonal_development!(f, m)
    end
end


"""
Determine whether an individual dies by background mortality, starvation ageing or chemical stress.
$(TYPEDSIGNATURES)
"""
function death!(f::Daphnid, m::Model)
    # background mortality
    if (f.juvenile)&(rand() > exp(- f.h_b / m.timestep_fleas))
        f.cause_of_death = 1.
        m.total_mortality += 1
        f.die = true
        @goto dead
    end
    
    # stochastic starvation mortality occurs if e <= e_starve
    if f.juvenile
        let h_e = f.e<=f.e_starve ? -(log(exp(log(f.e/f.e_starve)*f.beta_e))) : 0.
            if rand() > exp(-h_e / m.timestep_fleas)
                f.die=true
                f.cause_of_death = 2.1
                @goto dead
            end
        end
    end
	
	# deterministic starvation mortality occurs if individual lost a certain fraction of its structure
    #println([
    #    m.t_day, f.unique_id, f.V, f.V_m_hist, f.V_rel_crit
    #])
    if (f.juvenile)&(f.V<(f.V_m_hist*f.V_rel_crit))
        # TODO remove magic number
        if rand() < (1-(1-.35)^(1/m.timestep_fleas))
            f.cause_of_death = 2.2
            m.total_mortality += 1
            f.die = true
            @goto dead
        end
	end
	
    # ageing
    if rand() > ((1-min(f.h_a, 1))^(1/m.timestep_fleas))
        f.cause_of_death = 3.
        m.total_mortality += 1
        m.death_by_age += 1
        f.die = true
        @goto dead
    end

    # lethal toxicity
    if rand() > exp(-f.h_z/m.timestep_fleas)
        f.die = true
        m.total_mortality += 1
        m.death_by_stressor += 1
        f.cause_of_death = 4.
        @goto dead
    end

    # embryos will be killed at some point 
    # if they don't manage to reach the juvenile stage
    if (f.age>ABORT_INIT_EMB)&(f.juvenile==false)
        f.die=true
    end
    
    # experimental removal
    if  (f.juvenile) && (true in [floor(round(m.t_day, sigdigits=5))%7==x for x in m.remove_juveniles]) && (f.unique_id > m.N_0[f.species_idx])
        f.cause_of_death = 5.
        f.die = true
    end

    @label dead
end


function crowding_stress!(f::Daphnid)
    # TODO: validate implementation of crowding dose-responses
    # until then, crowding submodel is disontinued
    f.s_G_cr = max(0, 1 .+ f.stress_crowding[1])
    f.s_M_cr = max(0, 1 .+ f.stress_crowding[2])
    f.s_A_cr = max(0, f.stress_crowding[3])
    f.s_R_cr = max(0, f.stress_crowding[4])
end

"""
Application of stress (chemical, crowding) to the individual.
$(SIGNATURES)
"""
function stress!(f::Daphnid, m::Model)
    if f.juvenile
        # TKTD submodel takes some computation time 
        # -> execute only if any concentrations>0
        if sum(m.C)>0
            TKTD!(f, m)
        end
        """
        TODO: make crowding submodel use tuples to store parameter, not vector
        potential 100-fold difference in performance!
        f.stress_crowding = [LL3(
            m.V_density, 
            [
                f.crowding.s_EN50[i], 
                f.crowding.s_beta[i],
                f.crowding.s_max[i]
                ]) 
                for i in 1:length(f.crowding.s_EN50)
                    ]
         
        crowding_stress!(f)
        """
    end
end

"""
Execute a step of a vector of TK modes for a single stressor.
$(TYPEDSIGNATURES)
"""
function step!(TK::Vector{TKModel}, f::Daphnid, m::Model, C_ext::Float64)
    for tkmodel in TK
        step!(tkmodel, f, m, C_ext)
    end
end

"""
Execute a TK model step for a single stressor.
$(TYPEDSIGNATURES)
"""
function step!(TK::TKModel, f::Daphnid, m::Model, C_ext::Float64)
    # if needed retrieve additional internal and external states
    if (length(TK.add_states_int)>0)&&(length(TK.add_states_ext)>0)
        let add_states_int = Float64[getfield(f, x) for x in TK.add_states_int],
            add_states_ext = Float64[getfield(m, x) for x in TK.add_states_ext]
            # calculate change in internal concentration
            setfield!(
                TK, 
                :dC_int, 
                TK.model(
                    TK.C_int, 
                    C_ext, 
                    TK.params, 
                    add_states_int, 
                    add_states_ext
                    )
            )

        end
    # if not needed, save computation time by avoiding the array comprehensions
    else
        setfield!(
            TK, 
            :dC_int, 
            TK.model(
                TK.C_int, 
                C_ext, 
                TK.params, 
                EMPTVEC, 
                EMPTVEC
                )
        )
    end

    setfield!(TK, :C_int, max(0, TK.C_int + TK.dC_int/m.timestep_fleas))
end

"""
Execute a single TD model step.
$(TYPEDSIGNATURES)
"""
function step!(TD::TDModel, C_int::Float64, f::Daphnid)
    let current_state=getfield(f, TD.affected_state),
        y=response(TD.DRC, C_int)
        # the response is multiplied with the current value of the affected state
        setfield!(f, TD.affected_state, TD.apply_response((current_state, y)))
    end
    """
    # for all affected state variables
    for (app_resp,state) in zip(TD.apply_response, TD.affected_state)
        # retrieve the current state
        let current_state = getfield(f, state), new_state::Float64
            # apply response to current state
            setfield!(f, state, app_resp([current_state, y]))
        end
    end
    """
end

"""
Execute a TKTD model step for all stressors.
$(TYPEDSIGNATURES)
"""
function TKTD!(f::Daphnid, m::Model)
    # for every stressor
    for (z, TKTD) in enumerate(f.TKTD)
        if m.C[z]>0
            # calculate damage dynamics for each model
            step!(TKTD.TK, f, m, m.C[z])
        end
        # for every toxicodynamic component
        for (j,TD) in enumerate(TKTD.TD)
            # calculate the physiological response
            step!(TD, TKTD.TK[j].C_int, f)
        end
    end
end

"""
Some state variables have to be reset at the beginning of every time-step.
$(TYPEDSIGNATURES)
"""
function reset!(f::Daphnid)
    f.s_G_z = f.s_M_z = f.s_A_z = f.s_R_z = f.s_E_o_z = f.s_kappa_z = 1.0
    f.h_z = 0.0
end


"""
Apply different effects on metabolic fluxes (toxicity, crowding, adpative plasticity, resource P quota, resource correction factors).
$(TYPEDSIGNATURES)
"""
function apply_effects!(f::Daphnid)
    f.E_G = max(0, f.E_G_0 * f.s_G_z * f.s_G_cr)
    f.p_M_V = max(0, f.p_M_V_0 * f.s_M_z * f.s_M_cr)
    f.k_J = f.k_J_0 #max(0, f.k_J_0 * f.s_E_o_z * f.s_M_z)
    f.kappa_EX = max(0, min(1, (f.kappa_EX_max - (f.kappa_EX_max - f.kappa_EX_min) * f.f) * f.s_A_z * f.s_A_cr * f.q_p_rel * f.resource_correction))
    f.kappa_R = max(0, min(1, f.kappa_R_0 * f.s_R_z * f.s_R_cr))
    f.p_A_m = f.p_A_m_0
    f.J_X_A_m = f.J_X_A_m_0
    f.F_m = f.F_m_0

    if f.juvenile
        #f.kappa = (min(1, f.kappa_0*(1+(1/f.beta_e_adpt)*max(0, (1-(f.e))-(1-f.e_50_adpt)))))*f.s_kappa_z # adaptive plasticity following Gergs et al.
        f.kappa = (f.kappa_0 + (1-f.kappa_0)/(1+(f.e/f.e_50_adpt)^f.beta_e_adpt))*f.s_kappa_z
    else
        f.kappa = f.kappa_0
    end
end
