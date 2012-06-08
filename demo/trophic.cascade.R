# Demonstrates a trophic cascade using the Yodzis and Innes model.

# Three seperate simulations are run, each of a different community. 

# Each simulation uses the same model parameters.
# The first is of a single-producer, showing that the producer reaches 
# carrying capacity.
# The second is of a resource-consumer system, showing that the addition of a 
# herbivore results in a density of the producer much lower than carrying 
# capacity. 
# The third is of three-species chain, showing that a predator limits the 
# density of the herbivore and allows the producer to reach carrying capacity.

options(warn=2)
library(gruyere)

RunMySimulation <- function(spec, community, max.time=30, ylim=c(-3, 3), 
                            col=c(3, 1, 2))
{
    # A helper function that runs a model simulation and plots the resulting
    # biomass time series.

    # Collect simulation results in memory
    collector <- CollectChunksObserver()

    # Construct model parameters
    params <- IntermediateModelParams(community, spec)
    params <- BuildModelParams(community, params) # containing rho,x,z etc

    # No extinctions
    simulation <- LSODASimulation(model=YodzisInnesDyDt, 
                                  params=params, 
                                  sampling.interval=0.1, 
                                  atol=1e-20, 
                                  chunk.time=10)

    # Run the simulation until max.time is reached
    res <- RunSimulation(initial.state=Biomass(params$community), 
                         simulation=simulation,
                         controller=MaxTimeController(max.time=max.time), 
                         observers=list(collector))

    # Get the matrix of biomass time series
    tseries <- GetTimeSeries(collector)

    # Plot biomass time series    
    PlotBvT(community, tseries, category=NP(community, 'node'), col=col, 
            ylim=ylim)

    points(rep(0, NumberOfNodes(community)), Log10Biomass(community), pch=19, 
           col=c(3, 1, 2))

    # Mark carrying capacity
    abline(h=log10(spec['K']), col='grey', lty=2)
    mtext('K', side=4, at=log10(spec['K']), las=1, line=0.5)

    legend('topright', legend=NP(community, 'node'), col=col, lty=1)
}


# The same model parameters specification is used by all three simulations
spec <- ModelParamsSpec(f.constants=AllFConstantsEqual(0.1), K=100, B0=0.5, 
                        q=1)

par(mfcol=c(1,3))

# 1. Single producer
single.producer <- Community(nodes=data.frame(node=c('R'), 
                                              M=1, 
                                              N=20, 
                                              category='producer'), 
                                properties=list(title='Single producer', 
                                                M.units='kg', N.units='m^-2'))
RunMySimulation(spec, single.producer)

# 2. Resource-consumer motif
rc <- Community(nodes=data.frame(node=c('R','C'), 
                                 M=c(1,5), 
                                 N=c(20,1), 
                                 category=c('producer', 'invertebrate')),
                trophic.links=data.frame(resource='R', consumer='C'), 
                properties=list(title='Resource-consumer', 
                                M.units='kg', N.units='m^-2'))
RunMySimulation(spec, rc)

# 3. Three-species chain motif
three.species <- Community(nodes=data.frame(node=c('R','C','P'), 
                                            M=c(1,5,500), 
                                            N=c(20,1,0.1), 
                                            category=c('producer', 
                                                       rep('invertebrate', 2))),
                           trophic.links=data.frame(resource=c('R', 'C'), 
                                                    consumer=c('C', 'P')), 
                           properties=list(title='Three-species chain', 
                                           M.units='kg', N.units='m^-2'))
RunMySimulation(spec, three.species)

