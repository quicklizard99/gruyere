# Lotka-Volterra
# Demonstrates LSODASimulation(), RunSimulation() and visitors
LVDyDt <- function(time, y, params)
{
    R <- y[1]
    C <- y[2]

    with(params, 
    {
        return (list(c(R=r * R * (K-R)/K - a * R * C, 
                       C=e * e * R * C - d * C), 
                     globals=NULL))
    })
}

params <- list(r = 1.5, a=0.1, e=0.1, d=0.05, K=100)

simulation <- LSODASimulation(model=LVDyDt, 
                              params=params, 
                              sampling.interval=0.1,
                              use.atol=FALSE)

# Collect simulation results in memory
collector <- CollectChunksVisitor()

res <- RunSimulation(initial.state=c(R=100, C=5), 
                     simulation=simulation,
                     controller=MaxTimeController(200), 
                     visitors=list(collector, TimeSimulationVisitor()))

# Equilibria derived by hand (See Mick's aie herbivory slides 89 and 94
Re <- with(params, d / (e * a) )
Ce <- with(params, r * (K-Re) / (a * K))

# Equilibria calculated by Mathematica 8
stopifnot(isTRUE(all.equal(Re, with(params, d / (a*e)))))
stopifnot(isTRUE(all.equal(Ce, with(params, r * (a*e*K-d) / (e*K*a^2)))))

tseries <- get('tseries', collector)

# Plot the results
split.screen(c(1,2))
screen(1)

matplot(tseries[,1], log10(tseries[,-1]), type="l", lty=1, col=1:2)
abline(h=log10(Re), lty=2, col=1)
mtext(~R[e], side=4, at=log10(Re), las=1)
abline(h=log10(Ce), lty=2, col=2)
mtext(~C[e], side=4, at=log10(Ce), las=1)


screen(2, new=TRUE)
plot(log10(tseries[,2]), log10(tseries[,3]), xlab=~log[10](R), 
     ylab=~log[10](C), type="l", main="Consumer vs resource")
abline(v=log10(Re), lty=2, col=1)
mtext(~R[e], side=3, at=log10(Re), las=1, line=0.15)
abline(h=log10(Ce), lty=2, col=2)
mtext(~C[e], side=4, at=log10(Ce), las=1)

