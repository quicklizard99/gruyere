# Stuff that is useful and doesn't belong anywhere else
LoadMatrix <- function(path)
{
    # Returns a matrix 
    # check.names is FALSE in order to preserve spaces in column names
    names <- colnames(read.csv(path, nrows=1, check.names=FALSE))

    # scan() much faster than read.csv()
    ts <- matrix(scan(path, quiet=TRUE, sep=',', skip=1), 
                 ncol=length(names),
                 byrow=TRUE)
    colnames(ts) <- names
    return(ts)
}

.UTCTimeString <- function()
{
    # Textual representation of current time
    return(format(Sys.time(), '%Y-%m-%dT%H.%M.%S', tz="UTC", usetz=FALSE))
}

