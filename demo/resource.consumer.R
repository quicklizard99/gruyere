library(gruyere)
# Resource-consumer
community <- Community(nodes=data.frame(node=c('R','C'), 
                                        category=c('producer', 'invertebrate'), 
                                        M=c(0.1, 1), 
                                        N=c(100, 1)), 
                       trophic.links=data.frame(resource='R', consumer='C'), 
                       properties=c(title='Resource-consumer', 
                                    M.units='kg', 
                                    N.units='m^-2'))

# Model parameters
spec <- ModelParamsSpec(f.constants=AllFConstantsEqual())
params <- IntermediateModelParams(community, spec)
params <- BuildModelParams(community, params) # containing rho,x,z etc

simulation <- LSODASimulation(model=YodzisInnesDyDt, 
                              params=params, 
                              sampling.interval=0.1,
                              use.atol=FALSE)

# Collect simulation results in memory
collector <- CollectChunksVisitor()

res <- RunSimulation(initial.state=Biomass(community), 
                     simulation=simulation,
                     controller=MaxTimeController(200), 
                     visitors=list(collector, TimeSimulationVisitor()))

# Equilibria: eqns 12 and 13 of Yodzis and Innes (1992) on p.1160 using x and 
# y given in eqns 10 and 11, p 1156.
Re <- with(params, B0[1,2] / ( (y[1,2]-1)^ (1/(q+1))))
Ce <- as.numeric(with(params, (fe[1,2]*e[1,2] / x[2]) * Re * (1-Re/K)))

# Equilibria calculated by Mathematica 8. 
# q and d must be 0. 
Re.m <- with(params, B0[1,2] / (y[1,2]-1))
Ce.m <- with(params, -(B0[1,2]*e[1,2]*fe[1,2]*rho[1]*(K+B0[1,2]-K*y[1,2])) / 
                       ( K*x[2]*(y[1,2]-1)^2 ) )
Ce.m <- as.numeric(Ce.m)
stopifnot(isTRUE(all.equal(Re, Re.m)))
stopifnot(isTRUE(all.equal(Ce, Ce.m)))

tseries <- get('tseries', collector)

# Plot the results
par(mfrow=c(1,2))
PlotBvT(community, tseries, col=c(1,2))
abline(h=log10(Re), lty=2)
mtext(~R[e], side=4, at=log10(Re), las=1)
abline(h=log10(Ce), lty=2, col=2)
mtext(~C[e], side=4, at=log10(Ce), las=1)

plot(log10(tseries[,2]), log10(tseries[,3]), 
     xlab=Log10BLabel(community, name="italic(B[R])"), 
     ylab=Log10BLabel(community, name="italic(B[C])"), 
     type="l", main="Consumer vs resource")
abline(v=log10(Re))
mtext(~R[e], side=3, at=log10(Re), las=1, line=0.15)
abline(h=log10(Ce), col=2)
mtext(~C[e], side=4, at=log10(Ce), las=1)

