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
        r <- .Beta(E) * n       # Niche ranges
        c <- runif(S, min=r/2, max=pmin(n, 1-r/2))    # Niche centres

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

