# Functions for survival analysis
# Simon Hansul
# Laboratory of Environmental Toxicology And Aquatic Ecology
# Environmental Toxicology Unit
# Ghent University
# 2022-02-17

"""
Calculate survival probability from hazard rate and timespan.
$(TYPEDSIGNATURES)
"""
function survival_prob(h::Float64, t::Float64)
    exp(-h*t)
end

"""
Calculate hazard rate from survival probability and timespan.
$(TYPEDSIGNATURES)
"""
function hazard_rate(p::Float64, t::Float64)
    -(log(p)/t)
end

"""
Calculate ED50 from 21-day LC50 and maximum hazard rate.
$(TYPEDSIGNATURES)
"""
function calc_ED50(LC50, h_max)
    let haz_rate_LC50=hazard_rate(0.5, 21.)
        return LC50/((h_max/(haz_rate_LC50)-1)^(-1/2))
    end
end