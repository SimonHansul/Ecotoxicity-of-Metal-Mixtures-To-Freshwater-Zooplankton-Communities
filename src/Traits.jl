# Traits.jl
# Function definitions for various traits in relation to DEB parameters.

using DataStructures

"""
Maximum structural length of an individual.
"""
function L_max(v::Float64, k_M::Float64, g::Float64)
    return v/(k_M*g)
end

"""
Maximum carapace length of an individual.
$(TYPEDSIGNATURES)
"""
function L_max(v::Float64, k_M::Float64, g::Float64)
    return (1/shape_factor)*(v/(k_M*g))
end

"""
Maximum carapace length of an individual.
$(TYPEDSIGNATURES)
"""
function L_max(fleas_params::OrderedDict{String,Any})
    return (1/fleas_params["shape_factor"])*(fleas_params["v"]/(fleas_params["g"]*fleas_params["k_M"]))
end


"""
Returns dL/dt (cm/day) according to Von Bertalanffy function.
$(TYPEDSIGNATURES)
"""
function von_bertalanffy(L::Float64, L_max::Float64, k::Float64)
    return  k*(L_max-L)
end

"""
Returns length (cm) over time (days) according to Von Bertalanffy function.
$(TYPEDSIGNATURES)
"""
function vb_predicted_growth(L0, L_max, k, timepoints)
    Lt = [0]
    for t in timepoints[timepoints.>0]
        append!(Lt, Lt[end]+von_bertalanffy(Lt[end], L_max, k))
    end
    return DataFrame(
        t_day = timepoints,
        Length = Lt
    )
end