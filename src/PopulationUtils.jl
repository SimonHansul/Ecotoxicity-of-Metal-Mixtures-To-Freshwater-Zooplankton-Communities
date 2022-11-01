mutable struct PopulationDataset
    environment::DataFrame
    phyto_pop::DataFrame
    fleas_pop::DataFrame
    PopulationDataset() = new()
end