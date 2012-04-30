SSEModel <- function(N0, N.final)
{
    # Returns the sum of squares of differences between log10(N0) and 
    # log10(N.final).

    # N0 and N.final should be equal length vectors

    stopifnot(length(N0)==length(N.final))

    y <- N0         # Observed values
    f <- N.final    # Fitted values

    which.alive <- which(f>0)
    log10y <- log10(y[which.alive])
    log10f <- log10(f[which.alive])

    return (sum((log10y - log10f)^2))
}

ModelScore <- function(N0, N.final)
{
    # A very simple scoring system

    stopifnot(length(N0)==length(N.final))

    # Score = 
    #   + 1 for every species that persists
    #   + 1/(1+SSEmodel)

    sse.model <- SSEModel(N0, N.final)
    names(sse.model) <- NULL
    return (sum(N.final>0) + 1/(1+sse.model))
}

SSEModelOfScore <- function(score)
{
    # Returns the SSEmodel of a ModelScore
    fractional.part <- score %% 1
    return (1/fractional.part - 1)
}

SSERegression <- function(community)
{
    # Returns the sum of squares of residuals from a linear regression 
    # through log10(N) vs log10(M)
    m <- lm(log10(NP(community,'N')) ~ log10(NP(community,'M')))
    return (sum(m$residuals^2))
}

