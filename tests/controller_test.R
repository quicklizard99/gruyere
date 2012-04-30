TestRunningAverageController <- function()
{
    com <- Community(nodes=data.frame(node=paste('Species', 1:10), 
                                      category=rep('producer', 10), 
                                      M=rep(0.1, 10), 
                                      N=rep(10, 10)), 
                     properties=list(title='Testing', M.units='kg', 
                                     N.units='m^-2'))

    means.to.keep <- 10
    ct <- RunningAverageController(community=com, 
                                     equilibrium.fraction=0.01, 
                                     time.limit=10000, 
                                     means.to.keep=means.to.keep, 
                                     burn.in.time=0, 
                                     print.debug.msgs=FALSE)

    time.series <- matrix(NA, nrow=0, ncol=11)
    colnames(time.series) <- c('time', paste('Species', 1:10))
    means <- matrix(NA, ncol=10, nrow=0)

    # Check that running averages are being calculated correctly
    for(i in 0:20)
    {
        # Time series data
        chunk <- matrix(i+(1:10), ncol=10, nrow=2+i*2, byrow=TRUE)
        chunk <- cbind(seq(1+nrow(time.series), length.out=nrow(chunk)), chunk)
        colnames(chunk) <- colnames(time.series)

        result <- ControlOnChunk.RunningAverageController(ct, chunk)

        time.series <- rbind(time.series, chunk)

        stopifnot(all.equal(apply(time.series[,-1], 2, mean), 
                            unlist(data.frame(tail(get('means', ct), 1), 
                                              check.names=FALSE))))
        stopifnot(nrow(means)<=means.to.keep)
    }

    # Check that controller terminates
    for(i in 0:100)
    {
        chunk[,'time'] <- chunk[,'time'] + nrow(chunk)

        result <- ControlOnChunk.RunningAverageController(ct, chunk)
        if(result$terminate)
        {
            break
        }

        time.series <- rbind(time.series, chunk)

        stopifnot(all.equal(apply(time.series[,-1], 2, mean), 
                            unlist(data.frame(tail(get('means', ct), 1), 
                                              check.names=FALSE))))
        stopifnot(nrow(means)<=means.to.keep)
    }

    stopifnot(result$terminate)

    # Check that controller terminates if chunk is in equilibrium
    ct <- RunningAverageController(community=com, 
                                     equilibrium.fraction=0.01, 
                                     time.limit=10000, 
                                     means.to.keep=10, 
                                     burn.in.time=0, 
                                     print.debug.msgs=FALSE)

    chunk <- matrix(1:10, ncol=10, nrow=100, byrow=TRUE)
    chunk <- cbind(seq(1, length.out=nrow(chunk)), chunk)
    colnames(chunk) <- c('time', paste('Species', 1:10))

    result <- ControlOnChunk.RunningAverageController(ct, chunk)

    stopifnot(result$terminate)
    stopifnot(!result$means.equilibrium)
    stopifnot(result$chunk.equilibrium)


    # Check that the controller terminates with a short chunk in equilibrium
    ct <- RunningAverageController(community=com, 
                                     equilibrium.fraction=0.01, 
                                     time.limit=10000, 
                                     means.to.keep=10, 
                                     burn.in.time=0, 
                                     print.debug.msgs=FALSE)

    chunk <- matrix(1:10, ncol=10, nrow=10, byrow=TRUE)
    chunk <- cbind(seq(1, length.out=nrow(chunk)), chunk)
    colnames(chunk) <- c('time', paste('Species', 1:10))

    result <- ControlOnChunk.RunningAverageController(ct, chunk)

    stopifnot(result$terminate)
    stopifnot(!result$means.equilibrium)
    stopifnot(result$chunk.equilibrium)
}

