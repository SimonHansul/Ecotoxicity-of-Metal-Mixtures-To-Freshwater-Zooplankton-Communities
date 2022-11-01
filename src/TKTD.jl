# Data structures for TKTD modelling 
# Simon Hansul
# 2022-02-07
using Setfield
import Base.deepcopy
import Base.getproperty
include("DRC.jl")

mutable struct TKModel
    dC_int::Float64
    C_int::Float64
    add_states_int::Vector{Symbol}
    add_states_ext::Vector{Symbol}
    params::Vector{Float64}
    model::Function
end

mutable struct TDModel
    DRC::DRCModel
    affected_state::Symbol
    apply_response::Function
end

mutable struct TKTDModel
    TK::Vector{TKModel}
    TD::Vector{TDModel}
    link_TK::Bool
end

function init_TKTD!(tktd::TKTDModel)
    if tktd.link_TK
        [tki.params=tktd.TK[1].params for tki in tktd.TK]
    end
end

function init_TKTD!(tktdvec::Vector{TKTDModel})
    for tktd in tktdvec
        if tktd.link_TK
            [tki.params=tktd.TK[1].params for tki in tktd.TK]
        end
    end
end

"""
Minimal TK model.
$(TYPEDSIGNATURES)
"""
function minimaltk(
    C_int::Float64,
    C_ext::Float64,
    p::Vector{Float64},
    add_states_ind::Vector{Float64},
    add_states_ext::Vector{Float64}
    )
    return p[1]*(C_ext-C_int)
end


"""
One-parameter, size-dependent tokicokinetics.
$(TYPEDSIGNATURES)
"""
function oneparamtk(
    C_int::Float64, 
    C_ext::Float64, 
    p::Vector{Float64},
    add_states_int::Vector{Float64},
    add_states_ext::Vector{Float64}
    )
    let L=add_states_int[1], 
        dL= add_states_int[2], 
        L_m = add_states_int[3]
        return ((L_m/L)*p[1]*(C_ext-C_int))-(C_int*(3/L)*dL)
    end
end

"""
Assign parameters from vector to TKTD model. 
$(TYPEDSIGNATURES)
"""
function assign_params!(tktd::TKTDModel, params::Vector)
    # keep track of which element of params vector is currently being dealt with
    idx_count=1

    # assign TK parameters
    if length(tktd.TK)!=length(tktd.TD)
        error("Number of TK components has to match number of TD components.")
    end

    # for every TK component
    for (i,TK) in enumerate(tktd.TK)
        # if TK components are linked, reset the counter --> same parameters for all TK components
        tktd.link_TK ? idx_count=1 : nothing
        # for every parameter within TK component 
        for (j,_) in enumerate(TK.params)
            TK.params[j]=params[idx_count]
            idx_count += 1
        end
    end

    # assign TD parameters
    for (i,TD) in enumerate(tktd.TD)
        TD_params = zeros(length(TD.DRC.params))
        for(j,_) in enumerate(TD.DRC.params)
            TD_params[j]=params[idx_count]
            idx_count+=1
        end
        TD.DRC.params = tuple(TD_params...)
    end
    
    # if there are still parameters left, something's wrong.
    if idx_count<length(params)
        error("Length of params does not match up with number of TKTD params.")
    end
end

"""
Assign parameters from dictionary to TKTD model.
$(TYPEDSIGNATURES)
"""
function assign_params!(tktd::TKTDModel, params::OrderedDict{Symbol,Any})
    assign_params!(tktd, params.vals)
end

"""
Reset internal concentration of a TK model.
$(TYPEDSIGNATURES)
"""
function reset_internal_concentrations!(tk::TKModel)
    tk.C_int = 0.
end

"""
Reset all internal concentrations of a Vector of TK models.
$(TYPEDSIGNATURES)
"""
function reset_internal_concentrations!(tkvec::Vector{TKModel})
    for tk in tkvec
        reset_internal_concentrations!(tk)
    end
end

"""
Reset all internal concentrations of a Vector of TKTD models.
$(TYPEDSIGNATURES)
"""
function reset_internal_concentrations!(tktdvec::Vector{TKTDModel})
    for tktd in tktdvec
        reset_internal_concentrations!(tktd.TK)
    end
end

# Base function deepcopy is extended to achieve copying of the composite type

Base.deepcopy(tk::TKModel) = TKModel([deepcopy(getproperty(tk, field)) for field in fieldnames(TKModel)]...)
Base.deepcopy(tkvec::Vector{TKModel}) = [deepcopy(tk) for tk in tkvec]
Base.deepcopy(td::TDModel) = TDModel([deepcopy(getproperty(td, field)) for field in fieldnames(TDModel)]...)
Base.deepcopy(tdvec::Vector{TDModel}) = [deepcopy(td) for td in tdvec]
Base.deepcopy(tktd::TKTDModel) = TKTDModel([deepcopy(getproperty(tktd, field)) for field in fieldnames(TKTDModel)]...)
Base.deepcopy(tktdvec::Vector{TKTDModel}) = [deepcopy(tktd) for tktd in tktdvec]

#### We want the getfield/setfield functions to be used by default (faster)

function Base.getproperty(tktd::TKTDModel, f::Symbol)
    return getfield(tktd, f)
end

function Base.getproperty(tk::TKModel, f::Symbol)
    return getfield(tk, f)
end

function Base.getproperty(td::TDModel, f::Symbol)
    return getfield(td, f)
end