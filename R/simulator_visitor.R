# VisitStart is called by RunSimulation
VisitStart <- function(self, initial.state)
{
    UseMethod("VisitStart")
}

# VisitChunk is called by RunSimulation
VisitChunk <- function(self, chunk)
{
    UseMethod("VisitChunk")
}

# VisitEnd is called by RunSimulation
VisitEnd <- function(self, time, final.state)
{
    UseMethod("VisitEnd")
}


CurrentTimeVisitor <- function()
{
    # Prints current time
    self <- new.env(hash=TRUE, parent=emptyenv())
    class(self) <- c('CurrentTimeVisitor', class(self))
    return(self)
}

VisitStart.CurrentTimeVisitor <- function(self, initial.state)
{
    cat('Simulation started\n')
}

VisitChunk.CurrentTimeVisitor <- function(self, chunk)
{
    cat(paste('Simulation time', tail(chunk, 1)[1], '\n'))
}

VisitEnd.CurrentTimeVisitor <- function(self, time, final.state)
{
    cat(paste('Simulation finished at time', time, '\n'))
}


ExtinctionsFeedbackVisitor <- function()
{
    # Prints information about extinctions
    self <- new.env(hash=TRUE, parent=emptyenv())
    assign('extinct', NULL, self)
    class(self) <- c('ExtinctionsFeedbackVisitor', class(self))
    return(self)
}

.extinctions.feedback.worker <- function(self, state)
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

VisitStart.ExtinctionsFeedbackVisitor <- function(self, initial.state)
{
    assign('extinct', rep(FALSE, length(initial.state)), self)
    .extinctions.feedback.worker(self, initial.state)
}

VisitChunk.ExtinctionsFeedbackVisitor <- function(self, chunk)
{
    current.state <- tail(chunk, 1)[,-1]
    .extinctions.feedback.worker(self, current.state)
}

VisitEnd.ExtinctionsFeedbackVisitor <- function(self, time, final.state)
{
    .extinctions.feedback.worker(self, final.state)
    cat(paste(sum(get('extinct', self)), 'extinctions\n'))
}


CollectChunksVisitor <- function()
{
    # Collects all simulation chunks into one matrix
    self <- new.env(hash=TRUE, parent=emptyenv())
    assign('tseries', NULL, self)
    class(self) <- c('CollectChunksVisitor', class(self))
    return(self)
}

VisitStart.CollectChunksVisitor <- function(self, initial.state)
{
    # Colnames will be set in VisitChunk.CollectChunksVisitor?
    t <- matrix(c(0, initial.state), nrow=1)
    colnames(t) <- c('time', names(initial.state))
    assign('tseries', t, self)
}

VisitChunk.CollectChunksVisitor <- function(self, chunk)
{
    # Ignore the first row - it was recorded last time
    chunk <- tail(chunk, -1)
    assign('tseries', rbind(get('tseries', self), chunk), self)
}

VisitEnd.CollectChunksVisitor <- function(self, time, final.state)
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


WriteChunksVisitor <- function(fname)
{
    # Writes simulation chunks to disk
    self <- fname
    class(self) <- c('WriteChunksVisitor', class(self))
    return(self)
}

VisitStart.WriteChunksVisitor <- function(self, initial.state)
{
    x <- matrix(c(0, initial.state), nrow=1)
    colnames(x) <- c('time', names(initial.state))
    write.table(x=x, file=self, append=FALSE, row.names=FALSE, col.names=TRUE,
                sep=',')
}

VisitChunk.WriteChunksVisitor <- function(self, chunk)
{
    # Ignore the first row - it was written last time
    write.table(x=tail(chunk, -1), file=self, append=TRUE, row.names=FALSE, 
                col.names=FALSE, sep=',')
}

VisitEnd.WriteChunksVisitor <- function(self, time, final.state)
{
}


TimeBetweenChunksVisitor <- function(fname)
{
    # Prints the time between chunks
    self <- new.env(hash=TRUE, parent=emptyenv())
    assign('time', NULL, self)
    class(self) <- c('TimeBetweenChunksVisitor', class(self))
    return(self)
}

VisitStart.TimeBetweenChunksVisitor <- function(self, initial.state)
{
    assign('time', proc.time()['elapsed'], self)
}

VisitChunk.TimeBetweenChunksVisitor <- function(self, chunk)
{
    was <- get('time', self)
    now <- proc.time()['elapsed']
    cat(paste(now-was, 'seconds since last chunk\n'))
    assign('time', now, self)
}

VisitEnd.TimeBetweenChunksVisitor <- function(self, time, final.state)
{
}


PlotNDeviationsVisitor <- function(community, ...)
{
    self <- list(community=community, plot.params=list(...))
    class(self) <- c('PlotNDeviationsVisitor', class(self))
    return(self)
}

VisitStart.PlotNDeviationsVisitor <- function(self, initial.state)
{
    M <- NP(self$community,'M')
    do.call('PlotNDeviationsVisitor', c(list(community=self$community, 
                                             N.now=initial.state/M,
                                             main="time=0"), 
                                             self$plot.params))
}

VisitChunk.PlotNDeviationsVisitor <- function(self, chunk)
{
    time <- tail(chunk, 1)[1]
    current.state <- tail(chunk, 1)[2:ncol(chunk)]
    M <- NP(self$community,'M')
    do.call('PlotNDeviations', c(list(community=self$community, 
                                      N.now=current.state/M, 
                                      main=paste('time=', time, sep='')), 
                                      self$plot.params))
}

VisitEnd.PlotNDeviationsVisitor <- function(self, time, final.state)
{
    M <- NP(self$community,'M')
    do.call('PlotNDeviations', c(list(community=self$community, 
                                      N=final.state/M, 
                                      main=paste('final time=', time, sep='')), 
                                      self$plot.params))
}


PlotBvTVisitor <- function(community, ...)
{
    # Inherits from CollectChunksVisitor
    self <- CollectChunksVisitor()

    assign('community', community, self)
    assign('plot.params', list(...), self)

    class(self) <- c('PlotBvTVisitor', class(self))
    return(self)
}

VisitStart.PlotBvTVisitor <- function(self, initial.state)
{
    NextMethod()
}

VisitChunk.PlotBvTVisitor <- function(self, chunk)
{
    NextMethod()
    time <- tail(chunk, 1)[1]
    do.call('PlotBvT', c(list(community=get('community', self), 
                              tseries=get('tseries', self), 
                              main=paste('time=', time, sep='')), 
                              get('plot.params', self)))
}

VisitEnd.PlotBvTVisitor <- function(self, time, final.state)
{
    NextMethod()
    do.call('PlotBvT', c(list(community=get('community', self), 
                              tseries=get('tseries', self), 
                              main=paste('final time=', time, sep='')), 
                              get('plot.params', self)))
}


ElapsedTimeVisitor <- function()
{
    # Inherits from CollectChunksVisitor
    self <- new.env(hash=TRUE, parent=emptyenv())
    class(self) <- c('ElapsedTimeVisitor', class(self))
    return(self)
}

VisitStart.ElapsedTimeVisitor <- function(self, initial.state)
{
    assign('start', proc.time(), self)
}

VisitChunk.ElapsedTimeVisitor <- function(self, chunk)
{
}

VisitEnd.ElapsedTimeVisitor <- function(self, time, final.state)
{
    end <- proc.time()
    print('Simulation time:')
    print(end-get('start', self))
}

