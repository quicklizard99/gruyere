# SimStart is called by RunSimulation
SimStart <- function(self, initial.state)
{
    UseMethod("SimStart")
}

# SimChunk is called by RunSimulation
SimChunk <- function(self, chunk)
{
    UseMethod("SimChunk")
}

# SimEnd is called by RunSimulation
SimEnd <- function(self, time, final.state)
{
    UseMethod("SimEnd")
}


CurrentTimeObserver <- function()
{
    # Prints current time
    self <- new.env(hash=TRUE, parent=emptyenv())
    class(self) <- c('CurrentTimeObserver', class(self))
    return(self)
}

SimStart.CurrentTimeObserver <- function(self, initial.state)
{
    cat('Simulation started\n')
}

SimChunk.CurrentTimeObserver <- function(self, chunk)
{
    cat(paste('Simulation time', tail(chunk, 1)[1], '\n'))
}

SimEnd.CurrentTimeObserver <- function(self, time, final.state)
{
    cat(paste('Simulation finished at time', time, '\n'))
}


ExtinctionsFeedbackObserver <- function()
{
    # Prints information about extinctions
    self <- new.env(hash=TRUE, parent=emptyenv())
    assign('extinct', NULL, self)
    class(self) <- c('ExtinctionsFeedbackObserver', class(self))
    return(self)
}

.ExtinctionsFeedbackWorker <- function(self, state)
{
    # Does the job of recording and reporting new extinctions
    previous <- get('extinct', self)

    current <- 0==state
    new <- xor(previous, current)

    if(any(new))
    {
        cat(paste('Extinctions:', 
                    paste(which(new), 
                          names(new)[new], 
                          collapse=', '), '\n'))
        assign('extinct', current, self)
    }
}

SimStart.ExtinctionsFeedbackObserver <- function(self, initial.state)
{
    assign('extinct', rep(FALSE, length(initial.state)), self)
    .ExtinctionsFeedbackWorker(self, initial.state)
}

SimChunk.ExtinctionsFeedbackObserver <- function(self, chunk)
{
    current.state <- tail(chunk, 1)[,-1]
    .ExtinctionsFeedbackWorker(self, current.state)
}

SimEnd.ExtinctionsFeedbackObserver <- function(self, time, final.state)
{
    .ExtinctionsFeedbackWorker(self, final.state)
    cat(paste(sum(get('extinct', self)), 'extinctions\n'))
}


CollectChunksObserver <- function()
{
    # Collects all simulation chunks into one matrix
    self <- new.env(hash=TRUE, parent=emptyenv())
    assign('tseries', NULL, self)
    class(self) <- c('CollectChunksObserver', class(self))
    return(self)
}

SimStart.CollectChunksObserver <- function(self, initial.state)
{
    # Colnames will be set in SimChunk.CollectChunksObserver?
    t <- matrix(c(0, initial.state), nrow=1)
    colnames(t) <- c('time', names(initial.state))
    assign('tseries', t, self)
}

SimChunk.CollectChunksObserver <- function(self, chunk)
{
    # Ignore the first row - it was recorded last time
    chunk <- tail(chunk, -1)
    assign('tseries', rbind(get('tseries', self), chunk), self)
}

SimEnd.CollectChunksObserver <- function(self, time, final.state)
{
    # Clear rownames
    m <- get('tseries', self)
    rownames(m) <- NULL
    assign('tseries', m, self)
}

GetTimeSeries <- function(self)
{
    return (get('tseries', self))
}


WriteChunksObserver <- function(fname)
{
    # Writes simulation chunks to disk
    self <- fname
    class(self) <- c('WriteChunksObserver', class(self))
    return(self)
}

SimStart.WriteChunksObserver <- function(self, initial.state)
{
    x <- matrix(c(0, initial.state), nrow=1)
    colnames(x) <- c('time', names(initial.state))
    write.table(x=x, file=self, append=FALSE, row.names=FALSE, col.names=TRUE,
                sep=',')
}

SimChunk.WriteChunksObserver <- function(self, chunk)
{
    # Ignore the first row - it was written last time
    write.table(x=tail(chunk, -1), file=self, append=TRUE, row.names=FALSE, 
                col.names=FALSE, sep=',')
}

SimEnd.WriteChunksObserver <- function(self, time, final.state)
{
}


TimeBetweenChunksObserver <- function(fname)
{
    # Prints the time between chunks
    self <- new.env(hash=TRUE, parent=emptyenv())
    assign('time', NULL, self)
    class(self) <- c('TimeBetweenChunksObserver', class(self))
    return(self)
}

SimStart.TimeBetweenChunksObserver <- function(self, initial.state)
{
    assign('time', proc.time()['elapsed'], self)
}

SimChunk.TimeBetweenChunksObserver <- function(self, chunk)
{
    was <- get('time', self)
    now <- proc.time()['elapsed']
    cat(paste(now-was, 'seconds since last chunk\n'))
    assign('time', now, self)
}

SimEnd.TimeBetweenChunksObserver <- function(self, time, final.state)
{
}


PlotBDeviationsObserver <- function(community, ...)
{
    self <- list(community=community, plot.params=list(...))
    class(self) <- c('PlotBDeviationsObserver', class(self))
    return(self)
}

SimStart.PlotBDeviationsObserver <- function(self, initial.state)
{
    do.call('PlotBDeviations', c(list(community=self$community, 
                                       B.now=initial.state,
                                       main="time=0"), 
                                       self$plot.params))
}

SimChunk.PlotBDeviationsObserver <- function(self, chunk)
{
    time <- tail(chunk, 1)[1]
    current.state <- tail(chunk, 1)[2:ncol(chunk)]
    do.call('PlotBDeviations', c(list(community=self$community, 
                                      B.now=current.state, 
                                      main=paste('time=', time, sep='')), 
                                      self$plot.params))
}

SimEnd.PlotBDeviationsObserver <- function(self, time, final.state)
{
    do.call('PlotBDeviations', c(list(community=self$community, 
                                      B.now=final.state, 
                                      main=paste('final time=', time, sep='')), 
                                      self$plot.params))
}


PlotNDeviationsObserver <- function(community, ...)
{
    self <- list(community=community, plot.params=list(...))
    class(self) <- c('PlotNDeviationsObserver', class(self))
    return(self)
}

SimStart.PlotNDeviationsObserver <- function(self, initial.state)
{
    M <- NP(self$community,'M')
    do.call('PlotNDeviations', c(list(community=self$community, 
                                       N.now=initial.state/M,
                                       main="time=0"), 
                                       self$plot.params))
}

SimChunk.PlotNDeviationsObserver <- function(self, chunk)
{
    time <- tail(chunk, 1)[1]
    current.state <- tail(chunk, 1)[2:ncol(chunk)]
    M <- NP(self$community,'M')
    do.call('PlotNDeviations', c(list(community=self$community, 
                                      N.now=current.state/M, 
                                      main=paste('time=', time, sep='')), 
                                      self$plot.params))
}

SimEnd.PlotNDeviationsObserver <- function(self, time, final.state)
{
    M <- NP(self$community,'M')
    do.call('PlotNDeviations', c(list(community=self$community, 
                                      N.now=final.state/M, 
                                      main=paste('final time=', time, sep='')), 
                                      self$plot.params))
}


PlotBvTObserver <- function(community, ...)
{
    # Inherits from CollectChunksObserver
    self <- CollectChunksObserver()

    assign('community', community, self)
    assign('plot.params', list(...), self)

    class(self) <- c('PlotBvTObserver', class(self))
    return(self)
}

SimStart.PlotBvTObserver <- function(self, initial.state)
{
    NextMethod()
}

SimChunk.PlotBvTObserver <- function(self, chunk)
{
    NextMethod()
    time <- tail(chunk, 1)[1]
    do.call('PlotBvT', c(list(community=get('community', self), 
                              tseries=get('tseries', self), 
                              main=paste('time=', time, sep='')), 
                              get('plot.params', self)))
}

SimEnd.PlotBvTObserver <- function(self, time, final.state)
{
    NextMethod()
    do.call('PlotBvT', c(list(community=get('community', self), 
                              tseries=get('tseries', self), 
                              main=paste('final time=', time, sep='')), 
                              get('plot.params', self)))
}


PlotNvTObserver <- function(community, ...)
{
    # Inherits from CollectChunksObserver
    self <- CollectChunksObserver()

    assign('community', community, self)
    assign('plot.params', list(...), self)

    class(self) <- c('PlotNvTObserver', class(self))
    return(self)
}

SimStart.PlotNvTObserver <- function(self, initial.state)
{
    NextMethod()
}

SimChunk.PlotNvTObserver <- function(self, chunk)
{
    NextMethod()
    time <- tail(chunk, 1)[1]
    do.call('PlotNvT', c(list(community=get('community', self), 
                              tseries=get('tseries', self), 
                              main=paste('time=', time, sep='')), 
                              get('plot.params', self)))
}

SimEnd.PlotNvTObserver <- function(self, time, final.state)
{
    NextMethod()
    do.call('PlotNvT', c(list(community=get('community', self), 
                              tseries=get('tseries', self), 
                              main=paste('final time=', time, sep='')), 
                              get('plot.params', self)))
}


ElapsedTimeObserver <- function()
{
    # Inherits from CollectChunksObserver
    self <- new.env(hash=TRUE, parent=emptyenv())
    class(self) <- c('ElapsedTimeObserver', class(self))
    return(self)
}

SimStart.ElapsedTimeObserver <- function(self, initial.state)
{
    assign('start', proc.time(), self)
}

SimChunk.ElapsedTimeObserver <- function(self, chunk)
{
}

SimEnd.ElapsedTimeObserver <- function(self, time, final.state)
{
    end <- proc.time()
    print('Simulation time:')
    print(end-get('start', self))
}

