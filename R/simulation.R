RunChunk <- function(self, time, current.state)
{
    # Returns a chunk
    # Chunk is a matrix
    # First column is times. Time should be monotonically increasing
    # Other columns are systems states
    UseMethod("RunChunk")
}

LSODASimulation <- function(use.atol, ...)
{
    # use.atol was stupid - lsoda always uses atol
    # Swallow use.atol
    warning("LSODASimulation is deprecated. Please use 'ODESimulation' instead")
    return (ODESimulation(method='lsoda', ...))
}

ODESimulation <- function(model, params, sampling.interval=1, 
                          chunk.time=100, extinction.threshold=NULL, 
                          print.debug.msgs=FALSE, 
                          method='lsoda', atol, ...)
{
    # extinction.threshold should be in biomass units
    # extinction.threshold should be either single numbers or vectors - one 
    # number for every species
    # TODO Allow some extinction.threshold to be NA?
    stopifnot(is.null(extinction.threshold) || all(extinction.threshold>0))

    self <- new.env(hash=TRUE, parent=emptyenv())
    ode.params <- list(func=model, parms=params)
    ode.params <- c(ode.params, list(method=method), list(...))

    # Set lsoda's atol to be 1/10 of the extinction threshold(s) if appropriate
    if('lsoda'==method && !is.null(extinction.threshold) && missing(atol))
    {
        # Set lsoda's atol to be 1/10 of the extinction threshold(s)
        ode.params[['atol']] <- extinction.threshold
        if(print.debug.msgs)
        {
            print(paste("ODESimulation set lsoda's atol param to 1/10 of", 
                        "extinction.threshold"))
                  
        }
    }
    else if(!missing(atol))
    {
        ode.params[['atol']] <- atol
    }

    assign('ode.params', ode.params, self)

    assign('print.debug.msgs', print.debug.msgs, self)

    stopifnot(sampling.interval>0)
    stopifnot(chunk.time>0)
    simulation.settings <- list(sampling.interval=sampling.interval, 
                                chunk.time=chunk.time, 
                                extinction.threshold=extinction.threshold)
    assign('simulation.settings', simulation.settings, self)

    class(self) <- c('ODESimulation', class(self))
    return(self)
}

RunChunk.ODESimulation <- function(self, time, current.state)
{
    simulation.settings <- get('simulation.settings', self)
    times <- with(simulation.settings, 
                  seq(from=time, 
                      to=(time+chunk.time),
                      by=sampling.interval))

    args <- c(get('ode.params', self), list(times=times, y=current.state))
    chunk <- do.call('ode', args)

    chunk <- .ChunkEvents(chunk, simulation.settings$extinction.threshold, 
                          get('print.debug.msgs', self))

    stopifnot(nrow(chunk)>1)
    return (chunk)
}

.ChunkEvents <- function(chunk, extinction.threshold, print.debug.msgs=FALSE)
{
    # A private helper that handles extinctions and resurrections.
    # Returns a subset of chunk, up until the row at which the first 
    # extinction/resurrection occurred.

    # Keep the first n rows where n is the first row where:
    #    an extant species went lower that its extinction threshold
    #    an extinct species abundance was not 0
    # For each species that had an event on row n, report:
    #    extinction threshold
    #    abundance at n
    #    abundance at n+1
    # Some cases
    #    extinct at start of sim
    #    extinct at end of chunk

    # What if n is 1 or nrow(chunk)?

    # Get a matrix of bools; if TRUE, species was extinct at that time
    if(is.null(extinction.threshold))
    {
        # If no thresholds, check for -ve species
        went.extinct <- chunk[,-1,drop=FALSE] < 0
    }
    else if(1==length(extinction.threshold))
    {
        went.extinct <- chunk[,-1,drop=FALSE] < extinction.threshold
    }
    else if(length(extinction.threshold)==ncol(chunk)-1)
    {
        went.extinct <- t(apply(chunk[,-1], 1, '<', extinction.threshold))
    }
    else
    {
        stop(paste('extinction.threshold is of length', 
                   length(extinction.threshold), 'but should be either 1', 
                   'or', ncol(chunk)-1))
    }

    # Ignore species that are known to be extinct
    previously.extinct <- 0==chunk[1,-1]

    went.extinct[,previously.extinct] <- FALSE

    # Matrix of bools. TRUE if a species came back from the dead.
    resurrected <- apply(chunk[,-1,drop=FALSE], 
                         2, 
                         function(col) c(FALSE, 0==col[1] & col[-1]!=0))

    event <- went.extinct | resurrected

    if(any(event))
    {
        # A vector containing the row index of the first event for each 
        # species or NA if the species did not experience an event
        rows <- apply(event, 
                      2, 
                      function(col)
                      {
                          e <- which(col)
                          if(length(e))  return (min(e))
                          else           return (NA)
                      })

        # The first row on which an event occurred
        first.row <- min(rows, na.rm=TRUE)

        # Which species had an event on this first row
        species.with.events.on.row <- which(first.row == rows)

        if(print.debug.msgs)
        {
            # Print a matrix containing the time, state and extinction 
            # threshold of each relevant species
            print(paste('Event on row', first.row, 'time', 
                        chunk[first.row, 1], ':'))
            if(1==first.row)
            {
                fb <- chunk[first.row, c(1,1+species.with.events.on.row), 
                            drop=FALSE]
            }
            else
            {
                fb <- chunk[c(first.row-1, first.row), 
                            c(1,1+species.with.events.on.row), drop=FALSE]
            }

            if(!is.null(extinction.threshold))
            {
                if(length(extinction.threshold)>1)
                {
                    et <- extinction.threshold[species.with.events.on.row]
                    fb <- rbind(fb, et=c(NA, et))
                }
                else
                {
                    et <- rep(extinction.threshold, 
                              times=length(species.with.events.on.row))
                    fb <- rbind(fb, et=c(NA, et))
                }
            }
                
            print(fb)
        }

        # Do we allow resurrections?
        # I have seen extinct species resurrected with miniscule biomass 
        # abundances if atol is not set. See Dropbox/melt/R/resurrections.R for 
        # an example.
        resurrections.are.errors <- FALSE
        if(resurrections.are.errors && any(resurrected))
        {
            msg <- 'Species came back from the dead'
            save(msg, chunk, event, resurrected, extinction.threshold, 
                 file=paste(.UTCTimeString(), 'simulation.error.out'))
            stop(msg)
        }

        # Discard everything in this chunk after the event
        chunk <- head(chunk, first.row)

        # Set all species for which an event occurred on this row to 0. 
        # 1+ is required to offset time column
        chunk[first.row, 1+species.with.events.on.row] <- 0
    }

    return(chunk)
}

