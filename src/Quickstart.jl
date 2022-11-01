# Quickstart.jl
# Functins to quickly set up simulations of different types, using default settings.
# Simon Hansul
# 2021-07-03

"""
Set up simulation of a life-table experiment.
$(TYPEDSIGNATURES)
"""
function quickstart_lifetable(;nreps_sim::Int64=10, default_fleas_params=default_fleas.DP_AMP)
    global global_params = copy(default_global.zooplankton.lifetable)
    global phyto_params = copy(default_phyto.Rsub)
    global fleas_params = copy(default_fleas_params)
    fleas_params["h_a_accel"] = 0.
    update_params!(fleas_params; compound_parameters=false)
end

"""
Set up simulation of a community experiment.
$(TYPEDSIGNATURES)
"""
function quickstart_community()
    global global_params = copy(default_global.zooplankton.population_renewal)
    global phyto_params = copy(default_phyto.Rsub)
    global_params["N0_adult"] = [0, 0, 0]
    global_params["N0_juvenile"] = [4, 4, 4]
    global flea_params_vec = [
        copy(default_fleas.DL_AMP),
        copy(default_fleas.DM_AMP),
        copy(default_fleas.DP_AMP)
        ]
    [update_params!(f; compound_parameters=false) for f in flea_params_vec]
end