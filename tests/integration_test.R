# Test that bring together many components
TestYodzisInnesModelBestTL84 <- function()
{
    # Runs a Yodzis & Innes (1992)-style simulation of the dataset from Tuesday 
    # Lake, 1984. Values of 7 free model parameters are the best set found 
    # by the optimizations described in 'A cure for the plague of parameters: 
    # constraining models of complex population dynamics with physiological and 
    # community allometries'

    expected <- c('Nostoc sp.'=6.32470214001375e-05,
                  'Arthrodesmus sp.'=7.33590847782641e-05,
                  'Cryptomonas sp. 1'=5.06524501920096e-05,
                  'Cryptomonas sp. 2'=5.32004428142976e-05,
                  'Chroococcus dispersus'=4.99616989094858e-05,
                  'Closteriopsis longissimus'=6.659787582353e-05,
                  'Dinobryon bavaricum'=0.000111722581783931,
                  'Dinobryon cylindricum'=5.97414615329074e-05,
                  'Dactylococcopsis fascicularis'=7.3572753771713e-05,
                  'Dictyosphaerium pulchellum'=5.83148019663645e-05,
                  'Dinobryon sertularia'=6.77832725832823e-05,
                  'Dinobryon sociale'=5.71758486968916e-05,
                  'Glenodinium quadridens'=4.60273940737796e-05,
                  'Microcystis aeruginosa'=9.52786152994513e-05,
                  'Mallomonas sp. 1'=8.50046708768603e-05,
                  'Mallomonas sp. 2'=7.38240408633035e-05,
                  'Unclassified flagellates'=4.84309762993981e-05,
                  'Peridinium limbatum'=8.48147895468447e-05,
                  'Peridinium cinctum'=6.24073941736165e-05,
                  'Peridinium pulsillum'=5.29980673306457e-05,
                  'Peridinium wisconsinense'=8.91736840482382e-05,
                  'Chromulina sp.'=5.94396672822349e-05,
                  'Selenastrum minutum'=4.94211547526107e-05,
                  'Synedra sp.'=8.23372609528922e-05,
                  'Trachelomonas sp.'=5.1288707071927e-05,
                  'Ascomorpha eucadis'=1.66102054531371e-06,
                  'Synchaeta sp.'=4.85511176720362e-06,
                  'Bosmina longirostris'=0.000114271522408096,
                  'Conochilus (solitary)'=6.43212829199353e-05,
                  'Cyclops varians rubellus'=3.90865889607138e-05,
                  'Diaphanosoma leuchtenbergianum'=0.000143960319502823,
                  'Daphnia pulex'=0.00017186123740103,
                  'Filinia longispina'=3.79691852590964e-07,
                  'Conochiloides dossuarius'=5.66036794582023e-05,
                  'Gastropus stylifer'=3.7056692872143e-06,
                  'Holopedium gibberum'=0.00023636989366993,
                  'Kellicottia sp.'=3.69240402362384e-06,
                  'Keratella cochlearis'=3.03492353145989e-06,
                  'Keratella testudo'=4.60372628149274e-08,
                  'Leptodiaptomus siciloides'=7.61164584462832e-05,
                  'Orthocyclops modestus'=5.71065314179584e-05,
                  'Ploesoma sp.'=4.2153584057591e-06,
                  'Polyarthra vulgaris'=5.51568641814386e-06,
                  'Trichocerca multicrinis'=5.31256695294234e-06,
                  'Trichocerca cylindrica'=5.52216050847156e-06,
                  'Tropocyclops prasinus'=4.28437025191453e-05,
                  'Chaoborus punctipennis'=0.000174440065710481,
                  'Phoxinus eos'=0.000343442638999802,
                  'Phoxinus neogaeus'=0.000339221181547895,
                  'Umbra limi'=0.000354434227517679)

    RunSimWithModel <- function(community, model)
    {
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

        controller <- RunningAverageController(community, 
                                           equilibrium.fraction=0.01, 
                                           time.limit=1000000000, 
                                           means.to.keep=10, 
                                           burn.in.time=1000, 
                                           print.debug.msgs=FALSE)

        # Per-species extinction threshold are the biomass density of one 
        # individual per volume of water in Tuesday Lake (83337 m^3)
        extinction.threshold <- NP(community, 'M')/83337
        extinction.threshold[Externals(community)] <- 0
        simulation <- ODESimulation(model=model, 
                                  params=params, 
                                  sampling.interval=1, 
                                  chunk.time=100,
                                  extinction.threshold=extinction.threshold,
                                  print.debug.msgs=FALSE)
    
        collector <- CollectChunksObserver()
        observers <- list(collector)

        B0 <- Biomass(community)
        B0[Externals(community)] <- 0
        res <- RunSimulation(initial.state=B0,
                             simulation=simulation,
                             controller=controller, 
                             observers=observers)

        stopifnot(all.equal(res$final.state[names(expected)],
                            expected))
        stopifnot(c("terminate", "chunk.equilibrium", "means.equilibrium", 
                    "time", "final.state") == names(res))
        stopifnot(4400==res$time)
        stopifnot(res$chunk.equilibrium)
        stopifnot(!res$means.equilibrium)
        return (GetTimeSeries(collector))
    }

    data(TL84)
    community <- RemoveIsolatedNodes(TL84)
    a <- RunSimWithModel(community, YodzisInnesDyDt)
    b <- RunSimWithModel(community, YodzisInnesDyDt_R)
    stopifnot(all.equal(a, b))

    # As above but with the community ordered by decreasing body size.
    community <- OrderCommunity(RemoveIsolatedNodes(TL84), 'M', decreasing=TRUE)
    a <- RunSimWithModel(community, YodzisInnesDyDt)
    b <- RunSimWithModel(community, YodzisInnesDyDt_R)
    stopifnot(all.equal(a, b))

    # With an external input that has no input, no decay and no consumers
    properties <- CPS(TL84)
    properties$title <- "Tuesday Lake sampled in 1984, with detritus"
    nodes <- NPS(RemoveIsolatedNodes(TL84))[,c('node','M','N','category')]
    nodes <- rbind(nodes, data.frame(node='Detritus',M=NA,N=NA,category=''))
    tl <- TLPS(TL84)
    community <- Community(nodes=nodes, properties=properties, trophic.links=tl)
    a <- RunSimWithModel(community, YodzisInnesDyDt)
    b <- RunSimWithModel(community, YodzisInnesDyDt_R)
    stopifnot(all.equal(a, b))

    # With an external input that has no input, no decay and many consumers
    properties <- CPS(TL84)
    properties$title <- "Tuesday Lake sampled in 1984, with detritus"
    nodes <- NPS(RemoveIsolatedNodes(TL84))[,c('node','M','N','category')]
    nodes <- rbind(nodes, data.frame(node='Detritus',M=NA,N=NA,category=''))
    tl <- TLPS(TL84)
    tl <- rbind(tl, data.frame(resource='Detritus', 
                               consumer=NonBasalNodes(RemoveIsolatedNodes(TL84))))
    community <- Community(nodes=nodes, properties=properties, trophic.links=tl)
    a <- RunSimWithModel(community, YodzisInnesDyDt)
    b <- RunSimWithModel(community, YodzisInnesDyDt_R)
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

