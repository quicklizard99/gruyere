.Beta <- function(E)
{
    # Beta distribution given α=1 has the form (β*(1-x) ^ (β-1))
    # See Williams and Martinez. 2000. Simple rules yield complex food 
    # webs, caption for figure 1.
    b <- (1/E) - 1
    y <- runif(E)
    return (1 - (1-y)^(1/b))
}

NicheModelTrophicLinks <- function(C, S, nwebs=1, nodes=paste('Species', 1:S))
{
    # Connectance, C
    # Number of species, S
    stopifnot(0<C && C<1)
    stopifnot(1<S)
    stopifnot(0<nwebs)

    F <- function()
    {
        n <- sort(runif(S))
        E <- rep(2 * C, S)      # Expected values
        r <- .Beta(E) * n        # Niche ranges
        c <- (n-r/2) * runif(S) + r/2   # Niche centres - Dan's
        c <- runif(S, min=r/2, max=pmin(n, 1-r/2))    # Niche centres - mine

        stopifnot(all(c<n))

        # The upper and lower bounds of each species niche
        niche.lower <- c-r/2
        niche.upper <- c+r/2

        trophic.links <- lapply(1:S, function(i)
        {
            diet <- which(n>=niche.lower[i] & n<=niche.upper[i])
            if(length(diet)>0)
            {
                return (cbind(resource=nodes[diet], consumer=nodes[i]))
            }
        })

        return (data.frame(do.call('rbind', c(trophic.links)), 
                           stringsAsFactors=FALSE))
    }

    return (replicate(nwebs, F(), simplify=FALSE))
}

NicheModelPredationMatrix <- function(C, S, nwebs=1, 
                                      nodes=paste('Species', 1:S), 
                                      feasible.webs=TRUE)
{
    # Connectance, C
    # Number of species, S
    # Number of matrices, nwebs
    # feasible.webs - if TRUE, returns only predation matrices for which 
    #                 TrophicLevels() succeeds, indicating that the food web 
    #                 is energetically feasible. See TrophicLevels() for more 
    #                 information.

    stopifnot(0<C && C<1)
    stopifnot(1<S)
    stopifnot(0<nwebs)

    F <- function()
    {
        n <- sort(runif(S))
        E <- rep(2 * C, S)       # Expected values
        r <- .Beta(E) * n        # Niche ranges
        c <- (n-r/2) * runif(S) + r/2   # Niche centres - Dan's
        c <- runif(S, min=r/2, max=pmin(n, 1-r/2))    # Niche centres - mine

        stopifnot(all(c<n))

        # The upper and lower bounds of each species niche
        niche.lower <- c-r/2
        niche.upper <- c+r/2

        pm <- matrix(0, nrow=S, ncol=S)
        for(i in 1:S)
        {
            diet <- which(n>=niche.lower[i] & n<=niche.upper[i])
            pm[diet,i] <- 1
        }
        colnames(pm) <- rownames(pm) <- nodes
        return (pm)
    }

    webs <- NULL
    while(length(webs)<nwebs)
    {
        new.webs <- replicate(nwebs-length(webs), F(), simplify=FALSE)
        if(feasible.webs)
        {
            # Remove non-feasible webs
            new.webs <- new.webs[sapply(new.webs, PMIsEnergeticallyFeasable)]
        }
        webs <- c(webs, new.webs)
    }

    return (webs)
}

PMIsEnergeticallyFeasable <- function(pm)
{
    return (!all(is.na(MatrixInversionTrophicLevel(pm))))
}

MatrixInversionTrophicLevel <- function(pm)
{
    W <- t(pm)
    rs <- rowSums(W)
    W <- W/matrix(rs, ncol = ncol(W), nrow = nrow(W))
    W[0 == rs, ] <- 0
    I <- diag(ncol(W))
    result <- tryCatch(solve(I - W), error = function(e) e)
    if ("error" %in% class(result))
    {
        tl <- rep(NA, ncol(pm))
        names(tl) <- colnames(pm)
    }
    else
    {
        tl <- rowSums(result)
    }
    return (tl)
}

