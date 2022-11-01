
# Initialization functions for mutable structs within DEB-AQP.
# Simon Hansul
# 2021-11-11
# Laboratory of Environmental Toxicology and Aquatic Ecology
# Ghent University

include("Structures.jl")

"""
Initialize phytoplankton community in model `m` based on species-specific phytoplankton parameters.
$(TYPEDSIGNATURES)
"""
function initialize_phytoplankton_community!(m::Model, phyto_params::Vector{OrderedDict{Symbol,Any}})
    for (species_idx, phyto_params) in enumerate(phyto_params)
        p = Phyto()
        assign_params!(p, phyto_params)
        p.species_idx = species_idx
        p.A = m.A_0[species_idx]
        p.Q = p.A*p.q_init
        push!(m.phyto, p)
    end
end

"""
Initialize zooplankton community in model `m` based on species-specific zooplankton (DEB/TKTD) parameters.
$(TYPEDSIGNATURES)
"""
function initialize_zooplankton_community!(m::Model, flea_params::Vector{OrderedDict{Symbol,Any}})
    for (species_idx, flp) in enumerate(flea_params)
        for _ in 1:m.N_0[species_idx]
            f = Daphnid()
            assign_params!(f, copy(flp))
            init_TKTD!(f.TKTD) # make sure that TKTD models are correctly initialized
            m.unique_id_count += 1
            f.unique_id = m.unique_id_count
            f.species_idx = species_idx
            # calculate initial reserves from primary parameters
            g = (f.E_G_0*f.v)/(f.kappa_0*f.p_A_m_0)
            k_M = f.p_M_V_0/f.E_G_0
            L_b = (((f.E_H_b_0/f.p_A_m_0)*f.v)/(g*(1-f.kappa_0)))^(1/3)
            E_o = f.p_A_m_0*(((((1/(L_b*(g+E_SCALED_INIT_MOTHERS)^(1/3))))-(kooijmans_beta(g/(g+E_SCALED_INIT_MOTHERS))/((3g^(1/3)f.v)/(k_M))))^(-3))/f.v)
            # intialize embryo with E_o
            initialize_embryo!(f, E_SCALED_INIT_MOTHERS, E_o)
            init_random_age!(f, m)
            push!(m.fleas, f)
        end
    end
end

"""
Initialize a model object. <br>
Requires dictionary of global parameters and vectors of phytoplankton and zooplankton parameters.
$(TYPEDSIGNATURES)
"""
function initialize_model(
    global_params::OrderedDict{Symbol,Any}, 
    phyto_params::Vector{OrderedDict{Symbol,Any}},
    flea_params::Vector{OrderedDict{Symbol,Any}}
    )
    m = Model()
    assign_params!(m, global_params)
    
    m.t_step = 0
    m.t_day = 0
    m.unique_id_count = 0
    m.n_spec = length(flea_params)
    m.V_density = 0.
    m.cum_repro = 0
    m.total_mortality = 0
    m.death_by_starvation = 0
    m.death_by_age = 0
    m.death_by_stressor = 0
    m.fleas = Daphnid[]
    m.phyto = Phyto[]
    m.fleas_record_ind = Vector{Vector{Any}}[]
    m.fleas_record_pop = Matrix(undef, 0, 3)
    m.phyto_record = Vector{Vector{Any}}[]

    initialize_phytoplankton_community!(m, phyto_params)
    initialize_zooplankton_community!(m, flea_params)
    
    # initize records
    m.fleas_record_ind = Vector{Vector{Any}}[]
    m.fleas_record_pop = Matrix(undef, 0, 3)
    m.phyto_record = Vector{Vector{Any}}[]
    m.global_record = Vector{Vector{Any}}[]

    m.popsizes = SVector{m.n_spec,Int64}(repeat([0], m.n_spec))

    return m
end