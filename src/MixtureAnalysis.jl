# MixtureAnalysis.jl
# Computation of IA/CA predicted effects, either under the null hypothesis or for a given deviation factor.
# Simon Hansul
# Laboratory of Environmental Toxicology and Aquatic Ecology
# Ghent University
# 2021

"""
Get IA-predicted response, given a deviation factor a. <br>
Toxic units calculated internally.
$(TYPEDSIGNATURES)
"""
function IA(
        concentrations::Array{Float64,1},
        EC50s::Array{Float64,1},
        single_dose_responses::Array{Float64,1},
        a::Float64
    )
    
    # Step 1: calculate toxic units and relative contributions to toxicity
  
    TUs = concentrations ./ EC50s
    sumTU = sum(TUs)
    zi = TUs ./ sumTU
    
    # Step 2: calculate the deviation term
    
    G = a * prod(zi)
    # this is how it's written in the hochmuth et al paper
    # not quite the same?? above makes more sense
    # (a * prod(TUs))/(sum(TUs)^2) #

    # Step 3: calculate IA-predicted respose under the null hypothesis that a=0
    IA_null = prod(single_dose_responses)
    
    # Step 4: evaluate full IA equation
    IA_predicted = cdf(Normal(), quantile(Normal(), IA_null)+G)

    return IA_predicted
end

"""
Get IA-predicted response, given a deviation factor a. <br>
Toxic units have to be passed on as argument.
$(TYPEDSIGNATURES)
"""
function IA(
        TUs::Array{Float64,1},
        single_dose_responses::Array{Float64,1},
        a::Float64
    )
    
    # Step 1: calculate relative contributions to toxicity
    TUs = TUs[TUs.>1e-8]
    sumTU = sum(TUs)
    zi = TUs ./ sumTU
    
    # Step 2: calculate the deviation term
    G = a * prod(zi)

    # Step 3: calculate IA-predicted under the null hypothesis that a=0
    IA_null = prod(single_dose_responses)
    
    # for values >1, IA model is not defined because quantile function has to be applicable
    # setting IA_nul to 0.999 solves this problem while still providing accurate results
    if IA_null>=1
        IA_null=0.999
    end
    
    # Step 4: evaluate full IA equation
    IA_predicted = cdf(Normal(), quantile(Normal(), IA_null)+G)
    
    return IA_predicted

end

"""
Get CA-predicted response, given a deviation factor a, 
and given that all stressors are modelled by a two-parameter log-logistic function. 
Can be easily modified to account for different dose-response models.
$(TYPEDSIGNATURES)
"""
function CA(
        concentrations::Array{Float64,1},
        EC50s::Array{Float64,1},
        betas::Array{Float64,1},
        a::Float64
    )
    
    # Step 1: Calculate the deviation term

    TUs = concentrations ./ EC50s
    TUs = TUs[TUs.>1e-8]
    sumTU = sum(TUs)
    zi = TUs ./ sumTU
    G = a * prod(zi)
    
    # Step 2: Solve numerically

    function simulator_function(x::Array{Float64,1})
        return sum(@. concentrations/(EC50s*((x/(1-x))^(1/betas))))-exp(G)
    end
    
    true_value = 0. # we are looking for the root
    loss_function(predicted::Float64, observed::Float64) = (predicted-observed)^2 # use sum of squares as loss term
    
    @warn("Solver has to be updated - using old version of SMC_ABC. Use rejection ABC or SMCABC from ApproxBayes.jl instead.")
    # numerical solver
    _, accepted = SMC_ABC(
        OrderedDict{String,Any}("y"=>Uniform(0.,1.0)),
        true_value,
        simulator_function,
        loss_function;
        k=1,
        n_samples=100,
        n_samples_init=1000,
        progress_every=5
    )
    
    # function returns relative effect -> return relative response instead
    CA_predicted = 1 .- @subset(accepted, :rho.==minimum(:rho)).y[1]
    
    return CA_predicted
end