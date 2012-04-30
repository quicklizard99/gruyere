# 3 stages:
#   1) Parameters specification
#      Per metabolic category, growth model and functional response 
#      ModelParamsSpec()
#
#   2) Intermediate params
#      Allows you to make per-species changes to f* values
#      Williams et al 2007, eqs 2.8 - 2.12
#      IntermediateModelParams()
#
#   3) Model params
#      What the model functions use - x, y, rho etc
#      Williams et al 2007, eqs 2.14 - 2.18
#      BuildModelParams()

YodzisInnes92AValues <- function()
{
    warning("The new name for this function is 'YodzisInnes92AConstants'")
    return (YodzisInnes92AConstants())
}

YodzisInnes92AConstants <- function()
{
    # per-metabolic category values proposed by Yodzis and Innes 1992 Am Nat
    return (list(ar=c(producer=0.386),
                 aT=c(invertebrate=0.5, vert.ecto=2.3, vert.endo=54.9),
                 aJ=c(invertebrate=9.7, vert.ecto=8.9, vert.endo=89.2)))
}

BroseEtAl06AValues <- function()
{
    warning("The new name for this function is 'BroseEtAl06AConstants'")
    return (BroseEtAl06AConstants())
}

BroseEtAl06AConstants <- function()
{
    # per-metabolic category values proposed by Brose Et Al 2006 Ecol Lett
    # No constants for vertebrate endotherms
    # ar and aJ are not maxima but realised rates, so should be used with 
    # all f*=1, I think.
    return (list(ar=c(producer=1),
                 aT=c(invertebrate=0.314, vert.ecto=0.88),
                 aJ=c(invertebrate=2.512, vert.ecto=3.52)))
}

YodzisInnes92FValues <- function()
{
    warning("The new name for this function is 'YodzisInnes92FConstants'")
    return (YodzisInnes92FConstants())
}

YodzisInnes92FConstants <- function()
{
    # per-metabolic category values proposed by Yodzis and Innes, 1992
    return (list(fr=c(producer=0.1),
                 fT=c(invertebrate=1, vert.ecto=1, vert.endo=1),
                 fJ=c(invertebrate=0.3, vert.ecto=0.2, vert.endo=1)))
}

AllFValuesEqual <- function()
{
    warning("The new name for this function is 'AllFConstantsEqual'")
    return (AllFConstantsEqual())
}

AllFConstantsEqual <- function(v=1)
{
    # All set to the same value
    return (list(fr=c(producer=v),
                 fT=c(invertebrate=v, vert.ecto=v, vert.endo=v),
                 fJ=c(invertebrate=v, vert.ecto=v, vert.endo=v)))
}

.CheckModelParamsSpec <- function(spec)
{
    # A private helper that checks spec
    with(as.list(spec), 
    {
        stopifnot(1==length(e.producer))
        stopifnot(0<=e.producer && e.producer<=1)
        stopifnot(1==length(e.consumer))
        stopifnot(0<=e.consumer && e.consumer<=1)
        stopifnot(1==length(fe))
        stopifnot(0<=fe && fe<=1)
        stopifnot(1==length(B0))
        stopifnot(0<B0 && B0<Inf)
        stopifnot(1==length(d))
        stopifnot(0<=d && d<=Inf)
        stopifnot(1==length(q))
        stopifnot(0<=q && q<=Inf)
        stopifnot(1==length(K))
        stopifnot(0<K && K<Inf)
        stopifnot(1==length(a))
        stopifnot(0<=a && a<Inf)
    })
}

ModelParamsSpec <- function(a.constants=YodzisInnes92AConstants(),
                            f.constants=YodzisInnes92FConstants(), 
                            e.producer=0.45, 
                            e.consumer=0.85, 
                            fe=1, 
                            B0=1000, # Half-saturation biomass
                            d=0,     # Predator interference
                            q=0,     # Shape of response
                            K=500,   # Carrying capacity
                            a=1)
{
    p <- c(unlist(a.constants), unlist(f.constants), e.producer=e.producer, 
           e.consumer=e.consumer, fe=fe, B0=B0, d=d, q=q, K=K, a=a)
    .CheckModelParamsSpec(p)
    return (p)
}

SaveModelParamsSpec <- function(params.spec, dir)
{
    write.csv(do.call('cbind.data.frame', params.spec), 
              file.path(dir, 'params_spec.csv'), row.names=FALSE)
}

LoadModelParamsSpec <- function(dir)
{
    params.spec <- read.csv(file.path(dir, 'params_spec.csv'))
    stopifnot(1==nrow(params.spec))
    params.spec <- do.call(list, params.spec[1,,drop=FALSE])
    return(params.spec)
}

IsProducer <- function(community)
{
    return ('producer'==NP(community, 'category'))
}

Producers <- function(community)
{
    return (names(which(IsProducer(community))))
}

IsConsumer <- function(community)
{
    return ('producer'!=NP(community, 'category'))
}

Consumers <- function(community)
{
    return (names(which(IsConsumer(community))))
}

IntermediateModelParams <- function(community, spec)
{
    # Generates model parameters. Returns a list that contains the vectors:
    #   ar
    #   aT
    #   aJ
    #   fr
    #   fT
    #   fJ

    # the matrices:
    #   fJ
    #   fe
    #   e
    #   a   (growth model)
    #   B0  (growth model)
    #   d   (functional response)
    #   q   (functional response)

    # and the integer value:
    #   s

    # Is this community suitable for running this simulation?
    stopifnot('M' %in% NodePropertyNames(community))
    stopifnot('N' %in% NodePropertyNames(community))
    stopifnot(all(NP(community, 'category') %in% c('producer', 'invertebrate', 
                                                   'vert.ecto', 'vert.endo')))
    stopifnot('kg'==CP(community, 'M.units'))

    # The model requires at least one producer. Time is normalised to the 
    # producer with the smallest mass.
    stopifnot(length(Producers(community))>0)

    # Some basic checks on the parameters spec
    .CheckModelParamsSpec(spec)

    ParamVector <- function(constant)
    {
        # Returns a vector of the given constant (a character - 'ar', 'fr' etc)
        # One entry for each population in the community. Values are set using 
        # the metabolic categories. NA if constant is not appropriate for the 
        # metabolic category.
        indices <- paste(constant, '.', NP(community,'category'), sep='')
        v <- spec[indices]
        names(v) <- NP(community,'node')
        return (v)
    }

    with(as.list(spec), 
    {
        n <- NumberOfNodes(community)
        pm <- PredationMatrix(community)
        node <- NP(community, 'node')
        fJ <- matrix(ParamVector('fJ'), nrow=n, ncol=n, byrow=TRUE)
        fJ[0==pm] <- NA
        rownames(fJ) <- colnames(fJ) <- node

        e <- matrix(NA, nrow=n, ncol=n)
        colnames(e) <- rownames(e) <- node
        e[Producers(community),] <- e.producer
        e[Consumers(community),] <- e.consumer
        e[0==pm] <- NA

        fe <- matrix(fe, nrow=n, ncol=n)
        fe[0==pm] <- NA
        rownames(fe) <- colnames(fe) <- node

        # The smallest producer is the reference for normalizing time
        # Could have several populations with the same smallest mass - take 
        # the first of these.
        M <- NP(community, 'M')
        s <- order(M)[1]
        s <- Producers(community)[order(M[IsProducer(community)])[1]]
        s <- NodeNameIndices(community, s)

        # Functional response
        B0 <- matrix(B0, nrow=n, ncol=n)
        B0[pm==0] <- NA
        rownames(B0) <- colnames(B0) <- node

        d <- rep(d, n)
        names(d) <- node
        d[Producers(community)] <- NA

        # Growth model
        a <- matrix(a, nrow=n, ncol=n)
        rownames(a) <- colnames(a) <- node
        a[Consumers(community),] <- NA
        a[,Consumers(community)] <- NA

        return (list(ar=ParamVector('ar'),
                     aT=ParamVector('aT'),
                     aJ=ParamVector('aJ'),
                     fr=ParamVector('fr'),
                     fT=ParamVector('fT'),
                     fJ=fJ,
                     s=s,
                     e=e,
                     fe=fe, 
                     B0=B0, d=d, q=q,   # Functional response
                     K=K, a=a))         # Growth model
    })
}

BuildModelParams <- function(community, params, exponent=0.25)
{
    # Returns a list containing the vectors of model parameters:
    #   rho (would rather use ρ rather than rho but R recommends using ASCII)
    #   x
    #   e
    #   B0
    #   d
    #   q
    #   K
    #   a

    # and the matrices of model parameters:
    #   y
    #   fe

    # Also the following vectors, cached here for performance:
    #   n.species
    #   producers
    #   producers.c
    #   n.producers
    #   consumers
    #   consumers.c
    #   n.consumers

    # Checks on model parameters
    with(params, 
    {
        n <- NumberOfNodes(community)
        producers <- Producers(community)
        consumers <- Consumers(community)
        pm <- PredationMatrix(community)

        stopifnot(all(is.na(ar[consumers])))
        stopifnot(all(0<ar[producers] & ar[producers]<Inf))
        stopifnot(n==length(ar))

        stopifnot(all(is.na(aT[producers])))
        stopifnot(all(0<aT[consumers] & aT[consumers]<Inf))
        stopifnot(n==length(aT))

        stopifnot(all(is.na(aJ[producers])))
        stopifnot(all(0<aJ[consumers] & aJ[consumers]<Inf))
        stopifnot(n==length(aJ))

        stopifnot(all(is.na(fr[consumers])))
        stopifnot(all(0<=fr[producers] & fr[producers]<Inf))
        stopifnot(n==length(fr))

        stopifnot(all(is.na(fT[producers])))
        stopifnot(all(0<=fT[consumers] & fT[consumers]<Inf))
        stopifnot(n==length(fT))

        stopifnot(all(is.na(fJ[0==pm])))
        stopifnot(all(0<=fJ[1==pm] & 
                         fJ[1==pm]<=Inf))
        stopifnot(all(c(n,n) == dim(fJ)))

        stopifnot(1==length(s)) 
        stopifnot(0<s & s<=n) 

        stopifnot(all(is.na(e[0==pm])))
        stopifnot(all(0<=e[1==pm] & 
                         e[1==pm]<=1))
        stopifnot(all(c(n,n) == dim(e)))

        stopifnot(all(is.na(fe[0==pm])))
        stopifnot(all(0<=fe[1==pm] & 
                         fe[1==pm]<=1))
        stopifnot(all(c(n,n) == dim(fe)))

        stopifnot(all(is.na(B0[0==pm])))
        stopifnot(all(0<B0[1==pm] & 
                        B0[1==pm]<Inf))
        stopifnot(all(c(n,n) == dim(B0)))

        stopifnot(all(is.na(d[producers]))) 
        stopifnot(all(0<=d[consumers] & d[consumers]<=Inf))
        stopifnot(all(1==length(q)))
        stopifnot(0<=q & q<=Inf)
        stopifnot(all(0<K & K<Inf))

        stopifnot(all(is.na(a[consumers,consumers])))
        stopifnot(all(0<a[producers,producers] & 
                        a[producers,producers]<Inf))
        stopifnot(all(c(n,n) == dim(a)))
    })

    with(c(community, params), 
    {
        M <- NP(community, 'M')

        m.term <- ((M[s]/M)^exponent)
        rho <- ((fr * ar)/(fr[s] * ar[s])) * m.term
        x <- ((fT * aT)/(fr[s] * ar[s])) * m.term
        y <- fJ * 
             matrix(aJ, nrow=length(M), ncol=length(M), byrow=TRUE) / 
             matrix(fT*aT, nrow=length(M), ncol=length(M), byrow=TRUE)

        # Good to cache vectors of producers and consumers for performance.

        # R is 1-indexed
        producers <- NodeNameIndices(community, Producers(community))
        consumers <- NodeNameIndices(community, Consumers(community))

        # C is 0-indexed.
        producers.c <- as.integer(producers-1)
        if(length(consumers)>0)
        {
            consumers.c <- as.integer(consumers-1)
        }
        else
        {
            consumers.c <- 0    # Can't pass value of NULL as arg to .C()
        }

        # A multiplier to convert units of t' (model equations) to years
        # Homage to Yodzis and Innes, 2.13 (p 42):
        one.t.prime <- 1/as.numeric(fr[s]*ar[s]*(M[s]^-exponent))

        params <- list(community=community, 
                       n.species=NumberOfNodes(community), 
                       producers=producers, 
                       producers.c=producers.c, 
                       n.producers=length(producers),
                       consumers=consumers, 
                       consumers.c=consumers.c, 
                       n.consumers=length(consumers),
                       rho=rho, x=x, y=y, e=e, fe=fe, 
                       B0=B0, d=d, q=q, # Functional response
                       K=K, a=a,        # Growth model
                       # The following are not required by model equations but 
                       # are useful for subsequent analyses.
                       s=s, one.t.prime=one.t.prime)

        .CheckModelParams(params)

        return(params)
    })
}

.CheckModelParams <- function(params)
{
    # A private helper that checks model parameters
    with(params, 
    {
         pm <- PredationMatrix(community)

         stopifnot(all(0<=rho[producers] & rho[producers]<Inf))
         stopifnot(all(is.na(rho[consumers])))

         stopifnot(all(is.na(x[producers])))
         stopifnot(all(0<x[consumers] & x[consumers]<Inf))

         stopifnot(all(is.na(y[0==pm])))
         stopifnot(all(0<y[1==pm] & 
                         y[1==pm]<Inf))
         stopifnot(all(c(n.species,n.species) == dim(y)))

         stopifnot(all(is.na(e[0==pm])))
         stopifnot(all(0<=e[1==pm] & 
                          e[1==pm]<=1))

         stopifnot(all(is.na(fe[0==pm])))
         stopifnot(all(0<=fe[1==pm] & 
                          fe[1==pm]<=1))

         stopifnot(all(is.na(B0[0==pm])))
         stopifnot(all(0<B0[1==pm] & 
                         B0[1==pm]<Inf))

         stopifnot(all(is.na(d[producers])))
         stopifnot(all(0<=d[consumers] & d[consumers]<=Inf))
         stopifnot(all(0<=q & q<=Inf))
         stopifnot(all(0<K & K<Inf))

         stopifnot(all(is.na(e[0==pm])))
         stopifnot(all(0<=e[1==pm] & 
                          e[1==pm]<=1))

         stopifnot(all(0<a[producers,producers] & 
                         a[producers,producers]<Inf))
    })
}

SimulationDays <- function(params, t.prime)
{
    # Returns the number of days represented by the dimensionless t' units
    return (365.25 * t.prime * params$one.t.prime)
}

