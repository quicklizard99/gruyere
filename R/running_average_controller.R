RunningAverageController <- function(community, equilibrium.fraction, 
                                     time.limit, means.to.keep=10, 
                                     burn.in.time=1000, 
                                     print.debug.msgs=FALSE)
{
    # A controller that keeps track of running averages of each chunk. The 
    # means for up to means.to.keep chunks are kept.

    # The controller terminates the simulation if any of the following are 
    # TRUE:
    #   a) the running averages change by less than equilibrium.fraction
    #   b) the values within the current chunk change by less than 
    #      equilibrium.fraction

    # a) attempts to detect stability in oscillating systems
    # b) attempts to detect stability in system that reach stable equilibria

    # Chunks in which a species has gone extinct (i.e. abundance of 0) are 
    # not considered for b)

    # If time.limit is reached, an error is thrown.

    # The first burn.in.time time steps are not considered in the running 
    # averages. This is makes it possible to ignore large transients 
    # typically seen at the start of simulations.

    stopifnot(time.limit>0)
    stopifnot(means.to.keep>0)
    stopifnot(burn.in.time>=0)

    stopifnot(1==length(equilibrium.fraction) || 
              NumberOfNodes(community)==length(equilibrium.fraction))

    self <- new.env(hash=TRUE, parent=emptyenv())

    # Parameters
    assign('means.to.keep', means.to.keep, self)
    assign('equilibrium.fraction', equilibrium.fraction, self)
    assign('time.limit', time.limit, self)
    assign('burn.in.time', burn.in.time, self)
    assign('debug', print.debug.msgs, self)

    # State
    means <- matrix(NA, ncol=NumberOfNodes(community), nrow=0)
    colnames(means) <- NP(community, 'node')
    assign('means', means, self)    # A matrix of running averages
    assign('means.processed.rows', 0, self)   # The time represented by 'means'

    # A vector of bools - TRUE if a species has gone extinct
    extinct <- rep(FALSE, NumberOfNodes(community))
    names(extinct) <- NP(community, 'node')
    assign('extinct', extinct, self)

    class(self) <- c('RunningAverageController', class(self))
    return(self)
}

.DebugCat <- function(self, msg)
{
    # A private helper function
    if(get('debug', self))
    {
        cat(msg)
        cat('\n')
    }
}

.DebugPrint <- function(self, msg)
{
    # A private helper function
    if(get('debug', self))
    {
        print(msg)
    }
}

.UpdateExtinctions <- function(self, chunk)
{
    # Records which species are extinct. Returns a vector of booleans; one for
    # each species - TRUE for species which have gone extinct in this chunk
    previous <- get('extinct', self)
    current <- 0==tail(chunk, 1)[,-1]
    new <- xor(previous, current)
    if(any(new))
    {
        .DebugCat(self, paste('Extinctions:', 
                              paste(which(new), names(new)[new], 
                                    collapse=', ')))
        assign('extinct', current, self)
    }

    return (new)
}

ControlOnChunk.RunningAverageController <- function(self, chunk)
{
    # The ControlOnChunk override for RunningAverageController
    new.extinctions <- .UpdateExtinctions(self, chunk)

    chunk.equil <- .IsChunkInEquilibrium(self, chunk[,-1], new.extinctions)

    burn.in.time <- get('burn.in.time', self)
    if(head(chunk, 1)[1]>=burn.in.time)
    {
        # We have passed the first burn.in.time steps. This is hopefully 
        # long enough to have gone past the large transients that occur at 
        # at the start of most simulations.
        .AddChunkToMeans(self, chunk)
        means.equil <- .AreMeansInEquilibrium(self)
    }
    else
    {
        .DebugCat(self, paste('Not considering chunk means because', 
                              'first time in chunk of', head(chunk, 1)[1], 
                              'has not yet reached burn in time of', 
                              burn.in.time))
        means.equil <- FALSE
    }

    .DebugCat(self, paste('RunningAverageController ', 
                          'chunk lower time=', head(chunk, 1)[1], ', ', 
                          'chunk upper time=', tail(chunk, 1)[1], ', ', 
                          'means.processed.rows=', 
                          get('means.processed.rows', self), ', ', 
                          'chunk.equil=', chunk.equil, ', ', 
                          'means.equil=', means.equil, ', ', 
                          sep=''))

    time <- tail(chunk, 1)[1]

    # Evaluate stopping conditions
    if(chunk.equil)
    {
        # Return the averages over the this chunk only
        return (list(terminate=TRUE, 
                     chunk.equilibrium=TRUE, 
                     means.equilibrium=FALSE, 
                     time=time, 
                     final.state=apply(chunk[,-1], 2, mean)))
    }
    else if(means.equil)
    {
        # Return the last set of mean values

        # The final state is the current weighted means with extinct species 
        # set to zero
        final.state <- tail(get('means', self), 1)  # Matrix of 1 x N
        final.state <- final.state[1,]  # Convert to vector
        final.state[get('extinct', self)] <- 0

        return (list(terminate=TRUE, 
                     chunk.equilibrium=FALSE, 
                     means.equilibrium=TRUE, 
                     time=time, 
                     final.state=final.state))
    }
    else if(time >= get('time.limit', self))
    {
        stop(paste('Time limit of', get('time.limit', self), 'reached'))
    }
    else
    {
        return (list(terminate=FALSE, 
                     chunk.equilibrium=FALSE, 
                     means.equilibrium=FALSE, 
                     time=time, 
                     current.state=tail(chunk)[,-1]))
    }
}

.AddChunkToMeans <- function(self, chunk)
{
    # A private helper function that updates the matrix of means with the 
    # the given chunk.
    means <- get('means', self)

    chunk.no.time <- chunk[,-1]

    m.chunk <- apply(chunk.no.time, 2, mean)    # Mean for each species
    if(0==nrow(means))
    {
        mip1 <- m.chunk
    }
    else
    {
        # Weights
        # wi - the weighting for the means so far
        wi <- get('means.processed.rows', self)

        # wip1 - weighting for this chunk
        wip1 <- nrow(chunk)

        mi <- tail(means, 1, addrownums=FALSE)
        mip1 <- (wi*mi + wip1*m.chunk) / (wi+wip1)
    }

    # Add mip1 into the matrix of means. deparse.level=0 constructs no labels
    means <- rbind(means, mip1, deparse.level=0)

    # Set the means for extinct species to 0. - DON'T DO THIS - CHANGING 
    # HISTORY!
    # means[,which(0==tail(chunk.no.time,1))] <- 0

    # Limit the size of the matrix
    means <- tail(means, get('means.to.keep', self), addrownums=FALSE)

    stopifnot(all(!is.na(means)))
    stopifnot(all(!is.nan(means)))

    assign('means', means, self)
    assign('means.processed.rows', 
           nrow(chunk) + get('means.processed.rows', self), 
           self)
}

.AreMeansInEquilibrium <- function(self)
{
    # A private helper function that returns TRUE if the means are in 
    # equilibrium. At least means.to.keep rows of means must have been 
    # collected.
    means <- get('means', self)

    if(nrow(means)==get('means.to.keep', self))
    {
        equilibrium.fraction <- get('equilibrium.fraction', self)

        means.min <- apply(means, 2, min)
        means.max <- apply(means, 2, max)
        total.fractions <- (means.max-means.min) / means.max

        # means.max will be 0 for extinct species, giving NaN
        total.fractions[0==means.max] <- 0

        means.equil <- total.fractions < equilibrium.fraction

        # Extinct species are in equilibrium
        means.equil <- means.equil | get('extinct', self)

        if(all(means.equil))
        {
            .DebugCat(self, 'All means are in equilibrium')
        }
        else
        {
            .DebugCat(self, paste('Running averages for', 
                                  paste(which(!means.equil), collapse=','), 
                                  'not in equilibrium'))

            # Print out a table of means and fractions for all species that 
            # are not in equilibrium
            p <- rbind(means[,which(!means.equil),drop=FALSE], 
                       total.fractions[which(!means.equil)])
            rownames(p) <- rownames(p, do.NULL=FALSE, prefix='mean')
            rownames(p)[nrow(p)] <- 'fraction'
            .DebugPrint(self, p)
        }

        return (all(means.equil))
    }
    else
    {
        .DebugCat(self, paste('Not examining equilibrium of means as only', 
                              nrow(means), 
                              'sets of means collected so far'))
        return (FALSE)
    }
}

.IsChunkInEquilibrium <- function(self, chunk, new.extinctions)
{
    # A private helper function that returns TRUE if the chunk is in 
    # equilibrium.

    if(any(new.extinctions))
    {
        .DebugCat(self, paste('Not examining equilibrium for chunk as ', 
                              'extinction occurred; only ', nrow(chunk), 
                              ' rows', sep=''))
        return (FALSE)
    }
    else
    {
        equilibrium.fraction <- get('equilibrium.fraction', self)
        chunk.min <- apply(chunk, 2, min)
        chunk.max <- apply(chunk, 2, max)
        chunk.equil <- (chunk.max-chunk.min) <= equilibrium.fraction*chunk.max

        if(all(chunk.equil))
        {
            .DebugCat(self, 'The chunk is in equilibrium')
        }
        else
        {
            .DebugCat(self, paste('Chunk values for', 
                                  paste(which(!chunk.equil), collapse=','), 
                                  'are not in equilibrium'))
        }

        return (all(chunk.equil))
    }
}

