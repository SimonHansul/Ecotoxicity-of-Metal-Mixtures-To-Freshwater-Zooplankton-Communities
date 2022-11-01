#### Definition of constants


const L_O = 1e-5 # initial structural length (cm)
const V_O = L_O^3 # initial structural volume (cm^3)
const MG_ALG_TO_J = 21.84 # conversion from algal dry mass (mg) to energy (J) (J mg^-1)
const N_FLEAs_MAX = 5_000 # number of fleas at which simulation is terminated
const E_SCALED_INIT_MOTHERS = 1.0 # scaled reserve density of the mothers of initial embryos (-)
const W_ED = 2.0 # contribution of reserve and reproduction buffer to dry mass (g dry mass g C^-1)
const MU_E = 550_000 # chemical potential of reserve and reproduction buffer (J mol C^-1)
const RESERVE_TO_MASS = W_ED/MU_E # conversion factor from reserve (J) to dry mass (g)
const ABORT_INIT_EMB = 10 # time to abortion of intial embryos (d)
#const P_POTENT = 1221 # chemcical potential of reserve (J/mg P) NOT USED IN THIS MODEL VERSION

# attributes that are transferred to embryos
if isdefined(Main, :INHERITED_ATTRIBUTES)==false 
    const INHERITED_ATTRIBUTES = [
        :species,
        :species_idx,
        :cv,
        :J_X_A_m_0,
        :J_X_A_m,
        :p_A_m_0,
        :p_A_m,
        :F_m_0,
        :p_M_V_0,
        :E_G_0,
        :v,
        :kappa_0,
        :kappa,
        :kappa_R_0,
        :kappa_EX_min,
        :kappa_EX_max,
        :resource_corr_facts,
        :resource_correction,
        :h_a_accel,
        :beta_e,
        :shape_factor,
        :tau_molt,
        :h_b,
        :TKTD,
        :crowding,
        :s_G,
        :L_o,
        :V_o,
        :E_H_o,
        :E_H_b_0,
        :E_H_p,
        :e_50_adpt,
        :beta_e_adpt,
        :V_rel_crit
        ]
end
# attributes that are recorded in a flea_record
if isdefined(Main, :RECORDED_FLEA_ATTRIBS)==false
    const RECORDED_FLEA_ATTRIBS =  [
            :unique_id 
            :species_idx
            :kappa
            :kappa_EX 
            :q_p_rel 
            :J_X 
            :X_p 
            :f 
            :juvenile 
            :adult 
            :L 
            :dL 
            :L_b 
            :V 
            :dV 
            :E 
            :dE 
            :E_H 
            :E_R 
            :dE_R 
            :p_M
            :p_J 
            :p_A 
            :p_C
            :p_M_V
            :E_G
            :q_accel 
            :dq_accel
            :h_a
            :dh_a
            :clutch_size
            :k_M
            :k_J
            :g
            :U_H_b
            :U_H_p
            :K
            :e
            :E_o
            :molt_time
            :cum_repro
            :die
            :cause_of_death
            :h_z
            :s_G_z
            :s_M_z
            :s_A_z
            :s_R_z
            :s_G_cr
            :s_M_cr
            :s_A_cr
            :s_R_cr
            :age
            :mass
            ]
end