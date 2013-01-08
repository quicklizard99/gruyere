CalculateFlux <- function(community, tseries, flux.fn, flux.fn.params, 
                          from=1, from.time=NULL, 
                          to=nrow(tseries), to.time=NULL)
{
    # Returns a list of growth, respiration, consumption and assimilation 
    # values for the given time range. tseries should be a time series of 
    # biomass abundances.
    stopifnot(ncol(tseries)==1+NumberOfNodes(community))
    stopifnot(all(colnames(tseries)==c('time', NP(community, 'node'))))
    if(!is.null(from.time))
    {
        from <- which(tseries[,1]==from.time)
    }

    if(!is.null(to.time))
    {
        to <- which(tseries[,1]==to.time)
    }

    stopifnot(all.equal(1, length(to), length(from)))

    n.species <- NumberOfNodes(community)

    time <- from:to

    # Growth per species
    growth <- matrix(NA, ncol=n.species, nrow=length(time))

    # Respiration per species
    respiration <- matrix(NA, ncol=n.species, nrow=length(time))

    # Losses to consumption per trophic link
    consumption <- array(NA, dim=c(length(time), n.species, n.species), 
                         dimnames=c('time', 'Resource', 'Consumer'))
    # Total losses to consumption per species
    total.consumption <- matrix(NA, nrow=length(time), ncol=n.species)

    # Gains from consuming others per trophic link
    assimilation <- array(NA, dim=c(length(time), n.species, n.species), 
                         dimnames=c('time', 'Resource', 'Consumer'))

    # Total gains from consuming others per species
    total.assimilation <- matrix(NA, nrow=length(time), ncol=n.species)

    # Total flux per species
    total <- matrix(NA, nrow=length(time), ncol=n.species)

    node <- unname(NP(community, 'node'))
    colnames(growth) <- colnames(respiration) <- 
        colnames(total) <- 
        colnames(total.consumption) <- colnames(total.assimilation) <- 
        node
    dimnames(consumption) <- dimnames(assimilation) <- list(NULL, node, node)

    for(index in 1:length(time))
    {
        t <- time[index]
        f <- flux.fn(tseries[t,-1], flux.fn.params)

        growth[index,] <- f$growth
        respiration[index,] <- -f$respiration
        consumption[index,,] <- f$consumption
        assimilation[index,,] <- f$assimilation
        total.consumption[index,] <- rowSums(f$consumption, na.rm=TRUE)
        total.assimilation[index,] <- colSums(f$assimilation, na.rm=TRUE)

    }

    producers <- Producers(community)
    consumers <- Consumers(community)

    total[,producers] <- growth[,producers]
    total[,consumers] <- -respiration[,consumers]
    total <- total + total.assimilation - total.consumption
    return (list(time=tseries[from:to, 1], 
                 growth=growth, 
                 respiration=respiration,
                 consumption=consumption,
                 total.consumption=total.consumption,
                 assimilation=assimilation, 
                 total.assimilation=total.assimilation,
                 total=total))
}


