# Test that bring together many components
TestYodzisInnesModelBestTL84 <- function()
{
    # Runs a Yodzis & Innes (1992)-style simulation of the dataset from Tuesday 
    # Lake, 1984. Values of 7 free model parameters are the best set found 
    # by the optimizations described in 'A cure for the plague of parameters: 
    # constraining models of complex population dynamics with physiological and 
    # community allometries'

    data(TL84)

    # Take out isolated species
    community <- RemoveIsolatedNodes(TL84)

    # MODEL PARAMETERS

    # The best set of parameters found by the optimizer
    best <- c(fr.producer=0.597603718793107, 
              fJ.invertebrate=0.797592259166822, 
              fJ.vert.ecto=0.746630200690507, 
              d=1.17927240328204, 
              q=2.97244534827081, 
              B0=0.000159780748387173, 
              K=0.0337832616127192)

    # Three steps to assembling model parameters
    #   1) Parameters specification
    params.spec <- ModelParamsSpec(a.constants=YodzisInnes92AConstants(),
                                   f.constants=AllFConstantsEqual(), 
                                   e.producer=0.5, 
                                   e.consumer=0.5, 
                                   fe=1,    # Ingestion efficiency of 1 - 
                                            # typical of gape-feeders
                                   a=1)     # Producer competition

    # Use the best parameters
    params.spec[names(best)] <- best

    #   2) Intermediate params
    #      Williams et al 2007, eqs 2.8 - 2.12
    #      Allows you to make per-species changes
    params <- IntermediateModelParams(community, params.spec)

    #   3) Model params normalised to the growth rate of the primary producer 
    #      with the smallest body mass
    #      Williams et al 2007, eqs 2.14 - 2.18
    #      What the model functions use - x, y, rho etc
    params <- BuildModelParams(community, params)


    # CREATE THE OBJECTS THAT RUN THE SIMULATION

    # The controller object decides when the simulation ends. 
    # RunningAverageController() terminates simulations when either 
        #   a) the running averages over means.to.keep change by less than 
        #      equilibrium.fraction
        #   b) the values within the current chunk change by less than 
        #      equilibrium.fraction
    controller <- RunningAverageController(community, 
                                           equilibrium.fraction=0.01, 
                                           time.limit=1000000000, 
                                           means.to.keep=10, 
                                           burn.in.time=1000, 
                                           print.debug.msgs=FALSE)

    # LSODASimulation() runs a differential equation model, solved using, 
    # lsoda()which is provided by the deSolve package.
    # Per-species extinction threshold are the biomass density of one 
    # individual per volume of water in Tuesday Lake (83337 m^3)
    simulation <- LSODASimulation(model=YodzisInnesDyDt, 
                                  params=params, 
                                  sampling.interval=1, 
                                  chunk.time=100,
                                  extinction.threshold=NP(community, 'M')/83337,
                                  print.debug.msgs=FALSE)

    # observers is a list of observers that are shown the simulation as it runs. 
    # Cheddar contains several and users can write their own. 
    collector <- CollectChunksObserver()
    observers <- list(collector)

    # Start simulation at the empirical biomass densities
    res <- RunSimulation(initial.state=Biomass(community),
                         simulation=simulation,
                         controller=controller, 
                         observers=observers)

    # res is a list
    stopifnot(c("terminate", "chunk.equilibrium", "means.equilibrium", 
                "time", "final.state") == names(res))
    stopifnot(4400==res$time)
    stopifnot(res$chunk.equilibrium)
    stopifnot(!res$means.equilibrium)

    # res$final.state contains biomass abundances at the end of the simulation
    # Convert these into population densities
    N.final <- res$final.state / NP(community, 'M')

    # Compute model-data agreement values
    stopifnot(all.equal(27.26297023896913529484, 
                        SSEModel(NP(community, 'N'), N.final)))
    stopifnot(all.equal(50.03538198538740999766, 
                        ModelScore(NP(community, 'N'), N.final)))
}

TestTrophicCascade <- function()
{
    # Three seperate simulations, each of a different community. 

    # Each simulation uses the same model parameters.
    # The first is of a single-producer, showing that the producer reaches 
    # carrying capacity.
    # The second is of a resource-consumer system, showing that the addition of 
    # a herbivore results in a density of the producer much lower than carrying 
    # capacity. 
    # The third is of three-species chain, showing that a predator limits the 
    # density of the herbivore and allows the producer to reach carrying 
    # capacity.
    RunMySimulation <- function(spec, community, max.time=50)
    {
        # A helper function that runs a model simulation
        # Collect simulation results in memory
        collector <- CollectChunksObserver()

        # Construct model parameters
        params <- IntermediateModelParams(community, spec)
        params <- BuildModelParams(community, params) # containing rho,x,z etc

        # No extinctions
        simulation <- LSODASimulation(model=YodzisInnesDyDt, 
                                      params=params, 
                                      sampling.interval=0.1, 
                                      atol=1e-20)

        # Run the simulation until max.time is reached
        res <- RunSimulation(initial.state=Biomass(params$community), 
                             simulation=simulation,
                             controller=MaxTimeController(max.time=max.time))
        return (res)
    }


    # The same model parameters specification is used by all three simulations
    spec <- ModelParamsSpec(f.constants=AllFConstantsEqual(0.1), K=100, B0=0.5, 
                            q=0.2)


    # 1. Single producer
    single.producer <- Community(nodes=data.frame(node=c('R'), 
                                                  M=1, 
                                                  N=20, 
                                                  category='producer'), 
                                 properties=list(title='Single producer', 
                                                 M.units='kg', N.units='m^-2'))
    res <- RunMySimulation(spec, single.producer)
    stopifnot(100==res$time)
    stopifnot(all.equal(100, res$final.state))

    # 2. Resource-consumer motif
    rc <- Community(nodes=data.frame(node=c('R','C'), 
                                     M=c(1,5), 
                                     N=c(20,1), 
                                     category=c('producer', 'invertebrate')),
                    trophic.links=data.frame(resource='R', consumer='C'), 
                    properties=list(title='Resource-consumer motif', 
                                    M.units='kg', N.units='m^-2'))

    res <- RunMySimulation(spec, rc)
    stopifnot(100==res$time)
    stopifnot(all.equal(c(R=0.04325729349920991867, C=0.02449396363191617618),
                        res$final.state))

    # 3. Three-species chain motif
    three.species <- Community(nodes=data.frame(node=c('R','C','P'), 
                                        M=c(1,5,500), 
                                        N=c(20,1,1), 
                                        category=c('producer', 
                                                   rep('invertebrate', 2))),
                       trophic.links=data.frame(resource=c('R', 'C'), 
                                                consumer=c('C', 'P')), 
                       properties=list(title='Three-species chain motif', 
                                       M.units='kg', N.units='m^-2'))
    res <- RunMySimulation(spec, three.species)
    stopifnot(100==res$time)
    stopifnot(all.equal(c(R=98.32602570066444513941, 
                          C=0.04415248501705855422, 
                          P=2.17963044636725022940), 
                        res$final.state))
}

