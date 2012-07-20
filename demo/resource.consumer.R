# Runs a simulation of a simple resource-consumer community
community <- Community(nodes=data.frame(node=c('R','C'), 
                                        category=c('producer', 'invertebrate'), 
                                        M=c(0.1, 1), 
                                        N=c(100, 1)), 
                       trophic.links=data.frame(resource='R', consumer='C'), 
                       properties=c(title='Resource-consumer', 
                                    M.units='kg', 
                                    N.units='m^-2'))

# Model
model <- YodzisInnesDyDt

# Model parameters
spec <- ModelParamsSpec(f.constants=AllFConstantsEqual())
params <- IntermediateModelParams(community, spec)
params <- BuildModelParams(community, params) # containing rho,x,z etc

simulation <- ODESimulation(model=YodzisInnesDyDt, 
                            params=params, 
                            sampling.interval=0.1)

controller <- MaxTimeController(200)

collector <- CollectChunksObserver() # Collect simulation results in memory
observers <- list(collector, ElapsedTimeObserver())

res <- RunSimulation(initial.state=Biomass(community), 
                     simulation=simulation,
                     controller=controller, 
                     observers=observers)
res

tseries <- GetTimeSeries(collector)
head(tseries)

# Plot the results
par(mfrow=c(1,2))
PlotBvT(community, tseries, col=c(1,2))

# Equilibria: eqns 12 and 13 of Yodzis and Innes (1992) on p.1160 using x and 
# y given in eqns 10 and 11, p 1156.
Re <- with(params, W[1,2] / ( (y[1,2]-1)^ (1/(q+1))))
Ce <- as.numeric(with(params, (fe[1,2]*e[1,2] / x[2]) * Re * (1-Re/K)))
abline(h=log10(Re), lty=2)
mtext(~R[e], side=4, at=log10(Re), las=1, line=0.5)
abline(h=log10(Ce), lty=2, col=2)
mtext(~C[e], side=4, at=log10(Ce), las=1, line=0.5)

plot(log10(tseries[,'R']), log10(tseries[,'C']), 
     xlab=Log10BLabel(community, name="italic(B[R])"), 
     ylab=Log10BLabel(community, name="italic(B[C])"), 
     type="l", main="Consumer vs resource")
axis(side=3, labels=FALSE)
axis(side=4, labels=FALSE)
abline(v=log10(Re))
mtext(~R[e], side=3, at=log10(Re), las=1, line=0.5)
abline(h=log10(Ce))
mtext(~C[e], side=4, at=log10(Ce), las=1, line=0.5)

