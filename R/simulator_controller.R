ControlOnChunk <- function(self, chunk)
{
    # Should return a list containing the following:
    #   terminate - TRUE if the simulation should terminate
    #   result    - dependent on the controller
    UseMethod("ControlOnChunk")
}


MaxTimeController <- function(max.time)
{
    # Stops when time reaches or exceeds a maximum
    stopifnot(max.time>=0)
    self <- max.time
    class(self) <- c('MaxTimeController', class(self))
    return(self)
}

# The ControlOnChunk override for MaxTimeController
ControlOnChunk.MaxTimeController <- function(self, chunk)
{
    terminate <- tail(chunk, 1)[1]>=self
    if(terminate)
    {
        return (list(time=tail(chunk, 1)[,1], 
                     final.state=tail(chunk, 1)[,-1], 
                     terminate=tail(chunk, 1)[1]>=self))
    }
    else
    {
        return (list(time=tail(chunk, 1)[,1], 
                     current.state=tail(chunk, 1)[,-1], 
                     terminate=terminate))
    }
}


EquilibriumController <- function(max.time, equilibrium.fraction, 
                                  min.chunk.equil.rows=100, 
                                  print.debug.msgs=FALSE)
{
    # Stops either when species' numerical abundances change by less than a 
    # fraction or when a maximum time limit is reached.
    # Chunks containing less than min.chunk.equil.rows rows are not considered.
    stopifnot(max.time>=0)
    stopifnot(all(equilibrium.fraction>0, equilibrium.fraction<1))
    stopifnot(min.chunk.equil.rows>=10)

    self <- list(max.time=max.time, 
                 equilibrium.fraction=equilibrium.fraction, 
                 min.chunk.equil.rows=min.chunk.equil.rows, 
                 print.debug.msgs=print.debug.msgs)
    class(self) <- c('EquilibriumController', class(self))
    return(self)
}

# The ControlOnChunk override for EquilibriumController
ControlOnChunk.EquilibriumController <- function(self, chunk)
{
    stopifnot(1==length(self$equilibrium.fraction) || 
              length(self$equilibrium.fraction)==ncol(chunk)-1)

    if(nrow(chunk)>=self$min.chunk.equil.rows)
    {
        chunk.min <- apply(chunk[,-1], 2, min)
        chunk.max <- apply(chunk[,-1], 2, max)
        equil <- (chunk.max-chunk.min) <= self$equilibrium.fraction*chunk.max

        if(self$print.debug.msgs)
        {
            if(all(equil))
            {
                cat('All species are in equilibrium\n')
            }
            else
            {
                cat(paste(paste(which(!equil), collapse=','), 
                          'not in equilibrium\n'))
            }
        }
    }
    else
    {
        equil <- FALSE

        if(self$print.debug.msgs)
        {
            cat(paste('Not examining equilibrium for chunk as only', 
                      nrow(chunk), 'rows\n'))
        }
    }

    equil <- all(equil)

    time.limit <- tail(chunk, 1)[1] >= self$max.time
    terminate <- equil || time.limit
    state <- apply(chunk[,-1], 2, mean)

    if(terminate)
    {
        return (list(time=tail(chunk, 1)[,1], 
                     final.state=state, 
                     terminate=terminate, equilibrium=equil))
    }
    else
    {
        return (list(time=tail(chunk, 1)[,1], 
                     current.state=state, 
                     terminate=terminate, equilibrium=equil))
    }
}

