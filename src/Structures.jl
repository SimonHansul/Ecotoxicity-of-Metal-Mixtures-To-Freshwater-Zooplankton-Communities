# Structures.jl
# Definition of data structures for Zooplankton individuals, Phytoplankton populations and model objects.
# Simon Hansul
# 2021-02-16

using DataStructures, DataFrames
using StaticArrays
import Base:copy
include("TKTD.jl")

mutable struct SimulationOutput
    global_params::OrderedDict{String,Any}
    phyto_params_list::Array{OrderedDict{String,Any},1}
    flea_params_vec::Array{OrderedDict{String,Any},1}
    fleas_output_df::DataFrame
    phyto_output_df::DataFrame
end
  
mutable struct RecordInstruct
    output_type::String
    fields::Array{Symbol, 1} # fields to record
    how::Union{Nothing,Array{Function, 1}} # function to apply to array of values in field
    who::T where T <: Real # "whose" data to report; set to Inf to record data of all individual
    colnames::Array{String,1} # column names in final DataFrame
end


mutable struct Crowding
    s_max::Vector{Float64}
    s_EN50::Vector{Float64}
    s_beta::Vector{Float64}
end
copy(c::Crowding) = Crowding([copy(getproperty(c, field)) for field in fieldnames(Crowding)]...)


"""
Phytoplankton species oject.
$(TYPEDFIELDS)
"""
mutable struct Phyto
    species::String
    species_idx::Int64
    A::Float64
    A_sink::Float64
    Q::Float64
    Q_sink::Float64
    mu_max::Float64
    q_min::Float64
    q_max::Float64
    q_init::Float64
    v_max::Float64
    k_s::Float64
    m_max::Float64
    I_opt::Float64
    T_opt::Float64
    T_min::Float64
    T_max::Float64
    toxicity_params::Array{Vector{Float64}}
    sinking_rate::Float64
    resuspension_factor::Float64
    Phyto() = new()
end

"""
A mutable structure for individuals. <br>
Contains individual-level parameters and state variables. Initialize with values from `fleas_params` using function `make_flea`.
$(TYPEDEF)
"""
mutable struct Daphnid
    species::String # name of the species
    species_idx::Int64 # index of the species
    cv::Float64 # coefficient of variation to induce individual variability (-)
    p_A_m_0::Float64 # maximum surface-area specific ingestion rate (J cm^-2 d^-1)
    p_A_m::Float64 # maximum surface-area specific assimilation rate (J cm^-2 d^-1)
    F_m_0::Float64 # maximum surface-area specific filtration rate under reference conditions (L cm^-2 d^-1)
    F_m::Float64 # maximum surface-area specific filtration rate (L cm^-2 d^-1)
    p_M_V_0::Float64 # volume-specific somatic maintenance rate under reference conditions (J cm^-3 d^-1)
    p_M_V::Float64 # volume-specific somatic maintenance rate (J cm^-3 d^-1)
    E_G_0::Float64 # energy investment per unit of structural volume under reference conditions (J cm^-3)
    E_G::Float64 # energy investment per unit of structural volume (J cm^-3)
    v::Float64 # conductance rate (cm d^-1)
    kappa_0::Float64 # investment ratio to soma under reference conditions (-)
    kappa::Float64 # investment ratio to soma (-)
    kappa_R_0::Float64 # reproduction efficiency under reference conditions (-)
    kappa_R::Float64 # reproduction efficiency (-)
    kappa_EX_max::Float64 # maximum assimilation efficiency (-)
    kappa_EX_min::Float64 # minimum assimilation efficiency under reference conditions (-)
    kappa_EX::Float64 # assimilation efficiency (-)
    resource_corr_facts::Vector{Float64} # assimilation correction factors for different resources (-)
    h_a_accel::Float64 # ageing acceleration (d^-2)
    beta_e::Float64 # maximum hazard rate at reduced energy density (d^-1)
    e_starve::Float64 # level of scaled reserve density that will trigger starvation mortality
    shape_factor::Float64 # conversion from carapace length to structural length (-)
    tau_molt::Float64 # molting intervals (d)
    h_b::Float64 # hazard rate for background mortality (d^-1)
    TKTD::Vector{TKTDModel} # TKTD models
    crowding::Crowding # Crowding parameters (Crowding)
    unique_id::Int64 # unique identifier (-)
    resource_idx::Int64 # index of the resource currently used (-)
    resource_correction::Float64 # correction factor of the resource currently used (-)
    q_p_rel::Float64 # relative phosphorus saturation of the resource currently used (-)
    J_X::Float64 # ingestion rate (J d^-1)
    J_X_A_m_0::Float64 # maximum surface-area-specific ingestion rate under reference conditions (J cm^-2 d^-1)
    J_X_A_m::Float64 # maximum surface-area specific ingestion rate
    X_p::Float64 # density of the resource currently used (before feeding) (J/L)
    f::Float64 # scaled functional response (-)
    juvenile::Bool # indicates whether juvenile stage has been reached
    adult::Bool # indicates whether adult stage has been reached 
    L::Float64 # structural length (cm)
    dL::Float64 # change in structural length (cm d^-1)
    L_m::Float64 # maximum structural length (cm)
    L_b::Float64 # structural length at birth (cm)
    V::Float64 # structural volume (cm^3)
    dV::Float64 # change in structural volume (cm^3 d^-1)
    V_m::Float64 # maximum structural volume (cm^3)
    V_m_hist::Float64 # historical maximum structural volume (cm^3), i.e. maximum achieved by this individual
    V_rel_crit::Float64 # critical structural volume, relative to historical maximum (-)
    V_b::Float64 # structural volume at birth (cm^3)
    E::Float64 # reserve (J)
    dE::Float64 # change in reserve (J d^-1)
    E_H::Float64 # maturity (J)
    dE_H::Float64 # change in maturity (J d^-1)
    E_H_b_0::Float64 # maturity at birth under reference conditions (J)
    E_H_b::Float64 # maturity at birth (J)
    E_H_p::Float64 # maturity at puberty
    E_R::Float64 # reproduction buffer (J)
    dE_R::Float64 # change in reproduction buffer (J d^-1)
    p_M::Float64 # somatic maintenance rate (J d^-1)
    p_J::Float64 # maturity maintenance rate (J d^-1)
    p_A::Float64 # assimilation rate (J d^-1)
    p_C::Float64 # mobilization rate (J d^-1)
    q_accel::Float64 # ageing acceleration 
    dq_accel::Float64 # change in ageing acceleration (d^-1)
    h_a::Float64 # ageing hazard rate (d^-1)
    dh_a::Float64 # change in ageing hazard rate (d^-2)
    s_G::Float64 # gompertz stress coefficient (-)
    clutch_size::Int64 # cluth size (#)
    k_M::Float64 # somatic maintenance rate constant (d^-1)
    k_J_0::Float64 # maturity maintenance rate constant, 
    k_J::Float64 # maturity maintenance rate constant, calculated as k_J = k_M (d^-1)
    g::Float64 # energy investment ratio, calculated internally (-)
    U_H_b::Float64 # scaled maturity at birth, calculated internally (cm^2 d)
    U_H_p::Float64 # scaled maturity at puberty, calculated internally (cm^2 d)
    K::Float64 # half-saturation constant, calculated internally (J/L)
    e::Float64 # # scaled reserve density (-)
    L_o::Float64 # structural length at oviposition (cm)
    V_o::Float64 # structural volume at oviposition (cm^3)
    E_H_o::Float64 # maturity at oviposition (J)
    E_o::Float64 # reserves at oviposition (J)
    molt_time::Float64 # time since last molting (d)
    cum_repro::Float64 # cumulative reproduction (#)
    die::Bool # indicates whether individual will die at this timestep (-)
    cause_of_death::Float64 # encoding of what has caused death; 1=background, 2=starvation, 3=ageing, 4=lethal toxicity, 5=experimental removal
    h_z::Float64 # combined hazard rate (1/d) caused by chemical stressors z
    s_G_z::Float64 # combined somatic growth stress coefficient caused by stressors z
    s_M_z::Float64 # combined maintenance stress coefficient caused by stressors z
    s_A_z::Float64 # combined assimilation stress coefficient caused by stressors z
    s_R_z::Float64 # combined reproduction stress coefficient caused by stressors z
    s_E_o_z::Float64 # combined effect on maturity maintenance
    s_kappa_z::Float64 # combined effect on energy allocation
    s_G_cr::Float64 # somatic growth stress coefficient caused by crowding
    s_M_cr::Float64 # maintenance stress coefficient caused by crowding
    s_A_cr::Float64 # assimilation stress coefficient caused by crowding
    s_R_cr::Float64 # reproduction stress coefficient caused by crowding
    stress_crowding::Vector{Float64} # crowding stress per mode of action
    stress::Matrix{Float64} # chemical stres per mode of action and chemical
    age::Float64 # current age (since oviposition) in days
    embryos::Vector{Daphnid} # embryos in food pouch
    m_v::Float64 # dry mass of structure (g cm^-3)
    mass::Float64 # dry mass of the individual (g)
    e_50_adpt::Float64 # reserve density threshold for adaptive plasticity
    beta_e_adpt::Float64 # slope for adaptive plasticity
    Daphnid() = new()
end

"""
A mutable structure for models.
Contains model-level parameters and state variables. Intitialize with values from `global_params`. <br>
Note that all food levels are given in absolute amounts, not in concentrations. This is necessary to acurrately represent effects of resource limitation.
$(TYPEDFIELDS)
"""
mutable struct Model
    t_max::Float64 # maximum simulation time (d)
    data_recording_interval::Union{String,Float64} # intervals in which data is recorded (d). alternatively, use keyword "renewals" to record on data on renewal day
    record_flea_ind::Bool # record individual-level zooplankton data?
    record_flea_pop::Bool # record population/community-level zooplankton data?
    record_phyto::Bool # record phytoplankton data?
    remove_juveniles::Vector{Int64} # weekdays on which juveniles are removed from the experiment
    renew_medium::Vector{Int64} # weekdays on which medium is renewed (in most cases identical to remove_juveniles)
    frct_medium_renewal::Float64 # fraction of medium that is replaced at each renewal
    A_inoculate::Vector{Float64} # incoulation with algal biomass (mg dwt d^-1); absolute value, not per volume!
    P_inoculate::Float64 # inoculation with P (mg P d^-1); absolute value, not per volume!
    volume::Float64 # habitat volume in Liters
    timestep_fleas::Int64 # temporal resolution for Daphnids (timesteps per day)
    timestep_phyto::Int64 # temporal resolution for Phyto (timesteps per day); typically ca 10-times higher than for fleas due to faster metabolism
    t_step::Int64 # elapsed simulated time in model timesteps
    t_day::Float64 # elapsed simulated time in days
    age_init_max::Float64 # maximum flea's age at onset of experiment
    unique_id_count::Int64 # unique number of individuals throughout simulation
    fleas::Vector{Daphnid} # community of fleas
    n_spec::Int64 # number of species
    N_0::Vector{Int64} # initial number of fleas
    V_density::Float64 # structural volume density (cm^3/L), used for calculation of crowding effects
    cum_repro::Float64 # cumulative reproduction of all animals in the simulation (#)
    phyto::Vector{Phyto} # community of primary producers
    phyto_dynamic::Bool # indicates whether primary production is simulated dynamically or as a forcing variable
    phyto_constant::Bool # indicates wheter phytoplankton density is kept constant or ingestion by zooplankton is accounted for
    C::Vector{Float64} # vector of environmental chemical concentrations
    total_mortality::Int64 # cumulative total mortality of fleas
    death_by_starvation::Int64 # cumulative total mortality due to starvation
    death_by_age::Int64 # cumulative total mortality by ageing
    death_by_stressor::Int64 # cumulative total mortality by chemical stress
    A_0::Vector{Float64} # initial algal densities (mg)
    P_0::Float64 # initial Phosphorus in solution (mg)
    D::Float64  # dilution rate (day^-1)
    k::Vector{Float64} # toxicant degradation rates (day^-1)
    P::Float64 # P in solution (mg)
    T::Float64 # Temperature (°C)
    I::Float64 # light (µmol m^-2 s^-1)
    total_P::Float64 # sum of algal and dissolved P
    fleas_record_ind::Vector{Vector{Any}} # recorded individual state variables
    fleas_record_pop::Matrix{Any} # recorded population sizes
    popsizes::SVector # popsizes at the current time-step
    phyto_record::Vector{Vector{Any}} # recorded phytoplantkon state variables
    global_record::Vector{Vector{Any}} # recorded global state variables
    Model() = new()
end