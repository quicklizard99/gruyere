TestModelScore <- function()
{
    N0 <- 1:10
    Nfinal <- N0 + runif(N0, min=-0.1, max=0.1)

    sse.model <- as.numeric(SSEModel(N0, Nfinal))
    score <- ModelScore(N0, Nfinal)
    stopifnot(all.equal(SSEModelOfScore(score), sse.model))
}

