using Base: Float64
# DefaultInput.jl
# Default input parameters for AQP-DEB model. <br>
# Input dictionaries are organized in Named Tuples. <br>
# tpye ´default_´, then use code completion. <br>
# Simon Hansul
# 2021-10-29

using DocStringExtensions
using UnitfulMoles
include("Structures.jl")
include("DEB.jl")
include("TKTD.jl")

const P_INOCULATE_MG_DEFAULT = uconvert(u"mg", 3.0u"µmolP").val # default addition of P

#### ----- default global dictionaries ---- ####

"""
Collection of default model inputs .
"""
default_global = (
    zooplankton = (
        population_renewal = OrderedDict{Symbol,Any}(
            :data_recording_interval=> 1.,            
            :record_flea_ind => false,
            :record_flea_pop => true,
            :record_phyto => true,
            :remove_juveniles => Array{Int64,1}([]),
            :renew_medium => [0,3],
            :frct_medium_renewal => 0.25,
            :A_inoculate => [2.], 
            :P_inoculate =>  0.,
            :volume => 0.5,
            :t_max => 56.,
            :timestep_fleas => 24,
            :timestep_phyto => 24,
            :age_init_max => 1.,
            :N_0 => [12],
            :phyto_dynamic => false,
            :phyto_constant => false,
            :C => [0., 0., 0.],
            :A_0 => [0.],
            :P_0 => 0.,
            :D => 0.,
            :k => [0.],
            :T => 20.,
            :I => 100.
            ),
        lifetable = OrderedDict{Symbol,Any}(
            :data_recording_interval=> 1.,
            :record_flea_ind => true,
            :record_flea_pop => false,
            :record_phyto => true,
            :remove_juveniles => [0,2,5],
            :renew_medium => [0,2,5],
            :frct_medium_renewal => 1.,
            :A_inoculate => [.4], # mg dwt per day per species per vessel
            :P_inoculate => 0., # mgP added per day per vessel
            :volume => 0.04,
            :t_max => 21.,
            :timestep_fleas => 24,
            :timestep_phyto => 24,
            :age_init_max => 1.,
            :N_0 => [1],
            :phyto_dynamic => false,
            :phyto_constant => false,
            :C => [0., 0., 0.],
            :A_0 => [0.],
            :P_0 => 0.,
            :D => 0.,
            :k => [0.],
            :T => 20.,
            :I => 100.
        )),
    phytoplankton = (
        population_renewal = OrderedDict{Symbol,Any}(
            :data_recording_interval=> 1.,
            :record_flea_ind => false,
            :record_flea_pop => false,
            :record_phyto => true,
            :remove_juveniles => Array{Int64,1}([]),
            :renew_medium => [0,3],
            :frct_medium_renewal => 0.25,
            :A_inoculate => [0.], 
            :P_inoculate => P_INOCULATE_MG_DEFAULT,
            :volume => 0.5,
            :t_max => 56.,
            :timestep_fleas => 24,
            :timestep_phyto => 240,
            :age_init_max => 1.,
            :N_0 => [0],
            :phyto_dynamic => true,
            :phyto_constant => false,
            :C => [0., 0., 0.],
            :A_0 => [2.],
            :P_0 => 0.,
            :D => 0.,
            :k => [0.],
            :T => 20.,
            :I => 100.
            ),
        mono_chemostat = OrderedDict{Symbol,Any}(
            :data_recording_interval=> 1.,
            :record_flea_ind => false,
            :record_flea_pop => false,
            :record_phyto => true,
            :remove_juveniles => Array{Int64,1}([]),
            :renew_medium => [],
            :frct_medium_renewal => 0.,
            :A_inoculate => [0.], # mg dwt per day per species per vessel
            :P_inoculate =>  0., # mgP added per day per vessel
            :volume => 0.5,
            :t_max => 56.,
            :timestep_fleas => 24,
            :timestep_phyto => 240,
            :age_init_max => 1.,
            :N_0 => [0],
            :phyto_dynamic => true,
            :phyto_constant => false,
            :C => [0., 0., 0.],
            :A_0 => [2.],
            :P_0 => P_INOCULATE_MG_DEFAULT,
            :D => 0.5,
            :k => [0.],
            :T => 20.,
            :I => 100.
            ),
    ),
    two_trophic = (
        population_renewal = OrderedDict(
            :data_recording_interval=> 1.,
            :record_flea_ind => false,
            :record_flea_pop => true,
            :record_phyto => true,
            :remove_juveniles => Array{Int64,1}([]),
            :renew_medium => [0,3],
            :frct_medium_renewal => 0.25,
            :A_inoculate => [0.], # mg dwt per day per species per vessel
            :P_inoculate => P_INOCULATE_MG_DEFAULT, # mgP added per day per vessel
            :volume => 0.5,
            :t_max => 56.,
            :timestep_fleas => 24,
            :timestep_phyto => 240,
            :age_init_max => 1.,
            :N_0 => [12],
            :phyto_dynamic => true,
            :phyto_constant => false,
            :C => [0., 0., 0.],
            :A_0 => [2.],
            :P_0 => 0.,
            :D => 0.,
            :k => [0.],
            :T => 20.,
            :I => 100.
        ),
        mono_inoculate_algae = OrderedDict(
            :data_recording_interval=> 1.,
            :record_flea_ind => false,
            :record_flea_pop => true,
            :record_phyto => true,
            :remove_juveniles => Array{Int64,1}([]),
            :renew_medium => [0,3],
            :frct_medium_renewal => 0.25,
            :A_inoculate => [2.], # mg dwt per day per species per vessel
            :P_inoculate =>  uconvert(u"mg", 0.3u"µmolP").val, # mgP added per day per vessel
            :volume => 0.5,
            :t_max => 56.,
            :timestep_fleas => 24,
            :timestep_phyto => 240,
            :age_init_max => 1.,
            :N_0 => [12],
            :phyto_dynamic => true,
            :phyto_constant => false,
            :C => [0., 0., 0.],
            :A_0 => [2.],
            :P_0 => 0.,
            :D => 0.,
            :k => [0.],
            :T => 20.,
            :I => 100.
        )     
    ))

#### ---- default flea dictionaries ---- ####
# matrix indices for PMoAs
const G = 1
const M = 2
const A = 3
const R = 4
const H = 5
const EO = 6
const KP = 7

# matrix indices for stressors
const CU = 1
const NI = 2
const ZN = 3

#### default TD models 

if isdefined(Main, :TDG)==false
const TDG = TDModel(
    DRCModel( # a drc model object
        LL3GM, # dose-response function
        (0., 0., 0.) # parameters
        ),
        :s_G_z, # state variable modified by the dose-response function
        prod # function that is ussed to update the affected state
        )
end

if isdefined(Main, :TDM)==false
const TDM = TDModel(
    DRCModel(
        LL3GM,
        (0., 0., 0.)
        ),
        :s_M_z,
        prod
        )
end

if isdefined(Main, :TDA)==false
TDA = TDModel(
DRCModel(
    LL3AR,
    (0., 0., 0.)
    ),
:s_A_z,
prod
)
end

if isdefined(Main, :TDR)==false
const TDR = TDModel(
    DRCModel(
        LL3AR,
        (0., 0., 0.)
        ),
    :s_R_z,
    prod
    )
end

if isdefined(Main, :TDGUTS)==false
const TDGUTS = TDModel(
    DRCModel(
        LL3GUTS,
        (0., 0., 0.)
        ),
        :h_z,
        sum
        )
end

if isdefined(Main, :TDEO)==false
const TDEO = TDModel(
    DRCModel(
        LL3AR,
        (0., 0., 0.)
    ),
    :s_E_o_z,
    prod
)
end

if isdefined(Main, :TDkap)==false
const TDkap = TDModel(
    DRCModel(
        LL3AR,
        (0., 0., 0.)
        ),
        :s_kappa_z,
        prod
    )
end

if isdefined(Main, :default_TK)==false
    const default_TK = repeat([TKModel(0., 0., Symbol[], Symbol[], [0.], minimaltk)], 7)
end

if isdefined(Main, :default_TD)==false
    const default_TD = [TDG, TDM, TDA, TDR, TDGUTS, TDEO, TDkap]
end

if isdefined(Main, :default_TKTD)==false
    const default_TKTD = 
        TKTDModel(
            # TK component
            default_TK,
            # TD components
            default_TD,
            # separate TK by default
            false
            ) |> x -> [deepcopy(x) for _ in 1:3]
end

"""
Quick way to turn off crowding submodel, by setting fleas_params["crowding"] = zero_crowding
"""
zero_crowding = Crowding(
    Vector([0., 0., 0., 0.]),
    Vector([0., 0., 0., 0.]),
    Vector([0., 0., 0., 0.])
)

"""
Parameters for generic crowding submodel, based on Ban et al. data.
"""
parsim_crowding = Crowding(
    Vector([0., 0., 0.823, 0.]),
    Vector([0., 0., 0.212, 0.]),
    Vector([0., 0., 1.023, 0.])
)

"""
Parameters for generic crowding submodel, based on Ban et al. data.
"""
full_crowding = Crowding(
    Vector([2.46, 0.21, 0.65, 0.88]),
    Vector([0.28, 0.15, 0.12, 0.02]),
    Vector([8.1, 8.4, 7.1, 6.8])
)

generic_crowding = parsim_crowding #zero_crowding

default_fleas = (
    DM_AMP = OrderedDict{Symbol,Any}(
        :species => "DM_AMP",
        :cv => 0.05,
        :J_X_A_m_0 => 313.169/0.9, 
        :p_A_m_0 => 313.169,
        :F_m_0 => 30., 
        :E_H_b_0 => 0.05464, 
        :E_H_p => 1.09,
        :p_M_V_0 => 1200.,
        :E_G_0 => 4400.,
        :v => 0.1858,
        :kappa_0 => 0.5809,
        :kappa_R_0 => 0.95,
        :kappa_EX_max => 0.9,
        :kappa_EX_min => 0.5,
        :resource_corr_facts => [1.],
        :s_G => -0.3,
        :h_a_accel => 8e-5,#0.0002794, 
        :beta_e => 0.,
        :e_starve => 0.,
		:V_rel_crit => 0.4,
        :shape_factor => 0.221, 
        :tau_molt => 2.324,
        :h_b => 0.,
        :TKTD => deepcopy(default_TKTD),
        :crowding => zero_crowding,
        :m_v => 0.19, # dry mass of structure (g/cm^3),
        :e_50_adpt=>0.5,
        :beta_e_adpt=>2.0

    ),
    DP_AMP = OrderedDict{Symbol,Any}( 
        :species => "DP_AMP",
        :cv => 0.1,
        :J_X_A_m_0 => 520.825/0.9,
        :p_A_m_0 => 520.825,
        :F_m_0 => 6.5,
        :E_H_b_0 => 0.26, 
        :E_H_p => 3.926,
        :p_M_V_0 => 1599.,
        :E_G_0 => 4448.,
        :v => 0.04889,
        :kappa_0 => 0.4378,
        :kappa_R_0 => 0.95,
        :kappa_EX_max => 0.9,
        :kappa_EX_min => 0.5,
        :resource_corr_facts => [1.],
        :s_G => 0.0001,
        :h_a_accel => 8e-5,
        :beta_e => 0., 
        :e_starve => 0.,
		:V_rel_crit => 0.4,
        :shape_factor => 0.221,
        :tau_molt => 2.324,
        :h_b => 0.,
        :TKTD => deepcopy(default_TKTD),
        :crowding => zero_crowding,
        :m_v => 0.17,
        :e_50_adpt=>0.5,
        :beta_e_adpt=>2.0
    ),

    DL_AMP = OrderedDict{Symbol,Any}(
        :species => "DL_AMP",
        :cv => 0.05,
        :J_X_A_m_0 => 558.735/0.9,
        :p_A_m_0 => 558.735,
        :F_m_0 => 6.5,
        :E_H_b_0 => 0.01746, 
        :E_H_p => 1.401,
        :p_M_V_0 => 5686.83,
        :E_G_0 => 4400.59,
        :v => 0.020455,
        :kappa_0 => 0.5568,
        :kappa_R_0 => 0.95,
        :kappa_EX_max => 0.9,
        :kappa_EX_min => 0.5,
        :resource_corr_facts => [1.],
        :s_G => 0.001,
        :h_a_accel => 8e-5,
        :beta_e => 0,
        :e_starve => 0.,
		:V_rel_crit => 0.4,
        :shape_factor => 0.211,
        :tau_molt => 2.324,
        :h_b=>0.,
        :TKTD => deepcopy(default_TKTD),
        :crowding => zero_crowding,
        :m_v => 0.17, # dry mass of structure (g cm^-3)
        :e_50_adpt=>0.5,
        :beta_e_adpt=>2.0
        )
)
default_fleas[:DL_AMP][:TKTD][ZN].link_TK = false
default_fleas[:DP_AMP][:TKTD][CU].link_TK = false

default_phyto = (
    Rsub = OrderedDict{Symbol,Any}(
    :species => "Rsub",
    :mu_max => 1.78,
    :q_max => 0.018,
    :q_min => 0.0011,
    :q_init => 1., # this cancels out any P corrections
    :v_max => 0.0620,
    :k_s => 0.0625,
    :m_max => 0., 
    :I_opt => 120.,
    :T_opt => 28.,
    :T_max => 42.,
    :T_min => 0.,
    :toxicity_params => [[128.,1.199],[10e5,1.],[10e5,1.]],
    :sinking_rate => 0.,
    :resuspension_factor => 0.
    ),

    Dsub = OrderedDict{Symbol,Any}(
        :species => "Dsub",
        :mu_max => 1.74,
        :q_max => 0.014,
        :q_min => 0.0011,
        :q_init => 1., # this cancels out any P corrections
        :v_max => 0.052,
        :k_s => 0.068,
        :m_max => 0., # no background mortality
        :I_opt => 120.,
        :T_opt => 27.,
        :T_max => 35.,
        :T_min => 0.,
        :toxicity_params => [[115., 1.27],[10e5, 1.],[10e5,1.]],
        :sinking_rate => 0.5,
        :resuspension_factor => 1.
        )
)