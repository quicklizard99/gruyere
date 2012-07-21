# Test that bring together many components
TestYodzisInnesModelBestTL84 <- function()
{
    # Runs a Yodzis & Innes (1992)-style simulation of the dataset from Tuesday 
    # Lake, 1984. Values of 7 free model parameters are the best set found 
    # by the optimizations described in 'A cure for the plague of parameters: 
    # constraining models of complex population dynamics with physiological and 
    # community allometries'

    data(TL84)

    community <- RemoveIsolatedNodes(TL84)

    # The best set of parameters found by the optimizer
    best <- c(fr.producer=0.597603718793107, 
              fJ.invertebrate=0.797592259166822, 
              fJ.vert.ecto=0.746630200690507, 
              d=1.17927240328204, 
              q=2.97244534827081, 
              W=0.000159780748387173, 
              K=0.0337832616127192)

    params.spec <- ModelParamsSpec(a.constants=YodzisInnes92AConstants(),
                                   f.constants=AllFConstantsEqual(), 
                                   e.producer=0.5, 
                                   e.consumer=0.5, 
                                   fe=1,    # Ingestion efficiency of 1 - 
                                            # typical of gape-feeders
                                   a=1)     # Producer competition

    params.spec[names(best)] <- best
    params <- IntermediateModelParams(community, params.spec)
    params <- BuildModelParams(community, params)

    RunSimWithModel <- function(model)
    {
        controller <- RunningAverageController(community, 
                                           equilibrium.fraction=0.01, 
                                           time.limit=1000000000, 
                                           means.to.keep=10, 
                                           burn.in.time=1000, 
                                           print.debug.msgs=FALSE)

        # Per-species extinction threshold are the biomass density of one 
        # individual per volume of water in Tuesday Lake (83337 m^3)
        simulation <- ODESimulation(model=model, 
                                  params=params, 
                                  sampling.interval=1, 
                                  chunk.time=100,
                                  extinction.threshold=NP(community, 'M')/83337,
                                  print.debug.msgs=FALSE)
    
        collector <- CollectChunksObserver()
        observers <- list(collector)

        res <- RunSimulation(initial.state=Biomass(community),
                             simulation=simulation,
                             controller=controller, 
                             observers=observers)

        stopifnot(all.equal(unname(res$final.state), 
                        c(6.324702e-05,7.335908e-05,5.065245e-05,5.320044e-05,
                          4.996170e-05,6.659788e-05,1.117226e-04,5.974146e-05,
                          7.357275e-05,5.831480e-05,6.778327e-05,5.717585e-05,
                          4.602739e-05,9.527862e-05,8.500467e-05,7.382404e-05,
                          4.843098e-05,8.481479e-05,6.240739e-05,5.299807e-05,
                          8.917368e-05,5.943967e-05,4.942115e-05,8.233726e-05,
                          5.128871e-05,1.661020e-06,4.855112e-06,1.142715e-04,
                          6.432128e-05,3.908659e-05,1.439603e-04,1.718612e-04,
                          3.796918e-07,5.660368e-05,3.705669e-06,2.363699e-04,
                          3.692404e-06,3.034923e-06,4.603726e-08,7.611646e-05,
                          5.710653e-05,4.215358e-06,5.515687e-06,5.312567e-06,
                          5.522161e-06,4.284370e-05,1.744401e-04,3.434426e-04,
                          3.392212e-04,3.544342e-04), tolerance=1e-4))

        stopifnot(c("terminate", "chunk.equilibrium", "means.equilibrium", 
                    "time", "final.state") == names(res))
        stopifnot(4400==res$time)
        stopifnot(res$chunk.equilibrium)
        stopifnot(!res$means.equilibrium)
        return (GetTimeSeries(collector))
    }

    a <- RunSimWithModel(YodzisInnesDyDt)
    b <- RunSimWithModel(YodzisInnesDyDt_R)
    stopifnot(all.equal(a, b))
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
        simulation <- ODESimulation(model=YodzisInnesDyDt, 
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
    spec <- ModelParamsSpec(f.constants=AllFConstantsEqual(0.1), K=100, W=0.5, 
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

