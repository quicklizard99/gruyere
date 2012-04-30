TestChunkEvents <- function()
{
    # A test chunk
    chunk <- matrix(10, ncol=10, nrow=10)
    chunk[ 1, 1] <- 2     # Species 1, time 1
    chunk[ 1,10] <- 3     # Species 10, time 1
    chunk[10, 7] <- 8     # Species 7, time 10
    chunk[ 8, 5] <- 1     # Species 5, time 8
    chunk <- cbind(1:10, chunk) # time column
    colnames(chunk) <- c('time', paste('Species', 1:10))

    # No extinctions
    new.chunk <- gruyere:::.ChunkEvents(chunk, 0)
    stopifnot(all.equal(chunk, new.chunk))

    # Two extinctions in row 2
    new.chunk <- gruyere:::.ChunkEvents(chunk, 5)
    check.chunk <- head(chunk, 1)
    check.chunk[1, c('Species 1', 'Species 10')] <- 0
    stopifnot(all.equal(check.chunk, new.chunk))

    # Species 10 extinct in row 1
    new.chunk <- gruyere:::.ChunkEvents(chunk, c(0,0,0,0,0,0,0,0,0,4))
    check.chunk <- head(chunk, 1)
    check.chunk[1, c('Species 10')] <- 0
    stopifnot(all.equal(check.chunk, new.chunk))

    # Species 1 extinct in row 1
    new.chunk <- gruyere:::.ChunkEvents(chunk, c(4,0,0,0,0,0,0,0,0,0))
    check.chunk <- head(chunk, 1)
    check.chunk[1, c('Species 1')] <- 0
    stopifnot(all.equal(check.chunk, new.chunk))

    # Species 5 extinct in row 8
    new.chunk <- gruyere:::.ChunkEvents(chunk, c(0,0,0,0,2,0,0,0,0,0))
    check.chunk <- head(chunk,8)
    check.chunk[8, c('Species 5')] <- 0
    stopifnot(all.equal(check.chunk, new.chunk))

    # Species 7 extinct in row 10
    new.chunk <- gruyere:::.ChunkEvents(chunk, c(0,0,0,0,0,0,9,0,0,0))
    check.chunk <- head(chunk,10)
    check.chunk[10, c('Species 7')] <- 0
    stopifnot(all.equal(check.chunk, new.chunk))

    # Extinction in subsequent rows
    chunk <- matrix(10, ncol=4, nrow=5)
    chunk <- cbind(1:nrow(chunk), chunk) # time column
    colnames(chunk) <- c('time', paste('Species', 2:ncol(chunk) - 1))
    chunk[,2] <- 0
    chunk[2:nrow(chunk), 3] <- 0
    chunk[3:nrow(chunk), 4] <- 0
    chunk[4:nrow(chunk), 5] <- 0
    new.chunk <- gruyere:::.ChunkEvents(chunk, 1)
    stopifnot(all.equal(chunk[1:2,], new.chunk))

    new.chunk <- gruyere:::.ChunkEvents(chunk[2:nrow(chunk),], 1)
    stopifnot(all.equal(chunk[2:3,], new.chunk))

    new.chunk <- gruyere:::.ChunkEvents(chunk[3:nrow(chunk),], 1)
    stopifnot(all.equal(chunk[3:4,], new.chunk))

    # These tests disabled - resurrections are not errors
#    # Resurrection in second row
#    chunk <- matrix(10, ncol=4, nrow=5)
#    chunk <- cbind(1:nrow(chunk), chunk) # time column
#    colnames(chunk) <- c('time', paste('Species', 2:ncol(chunk) - 1))
#    chunk[,2] <- 0
#    chunk[2,2] <- 1
#    res <- tryCatch(gruyere:::.ChunkEvents(chunk, 1), error=function(e) return(TRUE))
#    stopifnot(res)

#    # Resurrection in last row
#    chunk <- matrix(10, ncol=4, nrow=5)
#    chunk <- cbind(1:nrow(chunk), chunk) # time column
#    colnames(chunk) <- c('time', paste('Species', 2:ncol(chunk) - 1))
#    chunk[,2] <- 0
#    chunk[5,2] <- 1
#    res <- tryCatch(gruyere:::.ChunkEvents(chunk, 10), error=function(e) return(TRUE))
#    stopifnot(res)
}

