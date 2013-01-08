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
    # Per-metabolic category values proposed by Brose Et Al 2006 Ecol Lett
    # No constants for vertebrate endotherms
    # ar and aJ are not maxima but realised rates, so should be used with 
    # all f*=1, I think.
    return (list(ar=c(producer=1),
                 aT=c(invertebrate=0.314, vert.ecto=0.88),
                 aJ=c(invertebrate=2.512, vert.ecto=3.52)))
}

OttoEtAl07AConstants <- function()
{
    # Per-metabolic category values proposed by Otto et al 2006 Ecol Lett
    # No constants for vertebrate ectotherms or endotherms
    # ar and aJ are not maxima but realised rates, so should be used with 
    # all f*=1, I think.
    return(list(ar = c(producer = 1), 
                aT = c(invertebrate = 0.2227), 
                aJ = c(invertebrate = 1.7816)))
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

AllFValuesEqual <- function(v=1)
{
    warning("The new name for this function is 'AllFConstantsEqual'")
    return (AllFConstantsEqual(v=v))
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
        stopifnot(1==length(W))
        stopifnot(0<W && W<Inf)
        stopifnot(1==length(d))
        stopifnot(0<=d && d<=Inf)
        stopifnot(1==length(q))
        stopifnot(0<=q && q<=Inf)
        stopifnot(1==length(K))
        stopifnot(0<K && K<Inf)
        stopifnot(1==length(a))
        stopifnot(0<=a && a<Inf)
        stopifnot(1==length(V))
        stopifnot(0<=V)
        stopifnot(1==length(Z))
        stopifnot(0<=Z)
    })
}

ModelParamsSpec <- function(a.constants=YodzisInnes92AConstants(),
                            f.constants=YodzisInnes92FConstants(), 
                            e.external=0.85, 
                            e.producer=0.45, 
                            e.consumer=0.85, 
                            fe=1, 
                            W=1000,  # Half-saturation biomass
                            B0,      # Old symbol for half-saturation biomass
                            d=0,     # Predator interference
                            q=0,     # Shape of response
                            K=500,   # Carrying capacity
                            a=1, 
                            V=0,
                            Z=0)
{
    if(!missing(B0))
    {
        warning(paste('B0 is the old symbol for functional response', 
                      'half-saturation biomass. Please use W in place of B0.'))
        B0 <- W
    }

    p <- c(unlist(a.constants), unlist(f.constants), e.external=e.external, 
           e.producer=e.producer, e.consumer=e.consumer, 
           fe=fe, W=W, d=d, q=q, K=K, a=a, V=V, Z=Z)
    .CheckModelParamsSpec(p)
    return (p)
}

.IsProducer <- function(community)
{
    return ('producer'==NP(community, 'category'))
}

.Producers <- function(community)
{
    return (names(which(.IsProducer(community))))
}

.IsConsumer <- function(community)
{
    category <- NP(community, 'category')
    return (category!='producer' & category!='')
}

.Consumers <- function(community)
{
    return (names(which(.IsConsumer(community))))
}

.IsExternal <- function(community)
{
    return (''==NP(community, 'category'))
}

.Externals <- function(community)
{
    return (names(which(.IsExternal(community))))
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
    #   W   (growth model)
    #   d   (functional response)
    #   q   (functional response)

    # and the integer value:
    #   m

    # Is this community suitable for running this simulation?
    stopifnot('M' %in% NodePropertyNames(community))
    stopifnot(all(!is.na(NP(community, 'M')[!.IsExternal(community)])))
    stopifnot(all(is.na(NP(community, 'M')[.IsExternal(community)])))
    stopifnot(all(NP(community, 'category') %in% c('producer', 'invertebrate', 
                                                   'vert.ecto', 'vert.endo', 
                                                   '')))
    stopifnot('kg'==CP(community, 'M.units'))

    # The model requires at least one producer. Time is normalised to the 
    # producer with the smallest mass.
    stopifnot(length(.Producers(community))>0)

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
        names(v) <- unname(NP(community,'node'))
        return (v)
    }

    with(as.list(spec), 
    {
        n <- NumberOfNodes(community)
        pm <- PredationMatrix(community)
        node <- unname(NP(community, 'node'))
        fJ <- matrix(ParamVector('fJ'), nrow=n, ncol=n, byrow=TRUE)
        fJ[0==pm] <- NA
        rownames(fJ) <- colnames(fJ) <- node

        e <- matrix(NA, nrow=n, ncol=n)
        colnames(e) <- rownames(e) <- node
        e[.Externals(community),] <- e.external
        e[.Producers(community),] <- e.producer
        e[.Consumers(community),] <- e.consumer
        e[0==pm] <- NA

        fe <- matrix(fe, nrow=n, ncol=n)
        fe[0==pm] <- NA
        rownames(fe) <- colnames(fe) <- node

        # The smallest producer is the reference for normalizing time
        # Could have several populations with the same smallest mass - take 
        # the first of these.
        M <- NP(community, 'M')
        m <- order(M)[1]
        m <- .Producers(community)[order(M[.IsProducer(community)])[1]]
        m <- NodeNameIndices(community, m)

        # Functional response
        W <- matrix(W, nrow=n, ncol=n)
        W[pm==0] <- NA
        rownames(W) <- colnames(W) <- node

        d <- rep(d, n)
        names(d) <- node
        d[!.IsConsumer(community)] <- NA

        # Growth model
        a <- matrix(a, nrow=n, ncol=n)
        rownames(a) <- colnames(a) <- node
        a[!.IsProducer(community),] <- NA
        a[,!.IsProducer(community)] <- NA

        V <- rep(V, n)
        names(V) <- node
        V[!.IsExternal(community)] <- NA

        Z <- rep(Z, n)
        names(Z) <- node
        Z[!.IsExternal(community)] <- NA

        return (list(ar=ParamVector('ar'),
                     aT=ParamVector('aT'),
                     aJ=ParamVector('aJ'),
                     fr=ParamVector('fr'),
                     fT=ParamVector('fT'),
                     fJ=fJ,
                     m=m,
                     e=e,
                     fe=fe, 
                     W=W, d=d, q=q,   # Functional response
                     K=K, a=a,        # Growth model
                     V=V, Z=Z))
    })
}

BuildModelParams <- function(community, params, exponent=1/4)
{
    # Returns a list containing the vectors of model parameters:
    #   rho (would rather use ρ rather than rho but R recommends using ASCII)
    #   x
    #   e
    #   W
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

    stopifnot(is.list(params))
    stopifnot(!is.null(names(params)))

    with(params, 
    {
        n <- NumberOfNodes(community)
        producers <- .Producers(community)
        consumers <- .Consumers(community)
        externals <- .Externals(community)
        pm <- PredationMatrix(community)

        stopifnot(all(is.na(ar[consumers])))
        stopifnot(all(is.na(ar[externals])))
        stopifnot(all(0<ar[producers] & ar[producers]<Inf))
        stopifnot(n==length(ar))

        stopifnot(all(is.na(aT[producers])))
        stopifnot(all(is.na(aT[externals])))
        stopifnot(all(0<aT[consumers] & aT[consumers]<Inf))
        stopifnot(n==length(aT))

        stopifnot(all(is.na(aJ[producers])))
        stopifnot(all(is.na(aJ[externals])))
        stopifnot(all(0<aJ[consumers] & aJ[consumers]<Inf))
        stopifnot(n==length(aJ))

        stopifnot(all(is.na(fr[consumers])))
        stopifnot(all(is.na(fr[externals])))
        stopifnot(all(0<=fr[producers] & fr[producers]<Inf))
        stopifnot(n==length(fr))

        stopifnot(all(is.na(fT[producers])))
        stopifnot(all(is.na(fT[externals])))
        stopifnot(all(0<=fT[consumers] & fT[consumers]<Inf))
        stopifnot(n==length(fT))

        stopifnot(all(is.na(fJ[0==pm])))
        stopifnot(all(0<=fJ[1==pm] & 
                         fJ[1==pm]<=Inf))
        stopifnot(all(c(n,n) == dim(fJ)))

        stopifnot(1==length(m)) 
        stopifnot(0<m & m<=n) 

        stopifnot(all(is.na(e[0==pm])))
        stopifnot(all(0<=e[1==pm] & 
                         e[1==pm]<=1))
        stopifnot(all(c(n,n) == dim(e)))

        stopifnot(all(is.na(fe[0==pm])))
        stopifnot(all(0<=fe[1==pm] & 
                         fe[1==pm]<=1))
        stopifnot(all(c(n,n) == dim(fe)))

        stopifnot(all(is.na(W[0==pm])))
        stopifnot(all(0<W[1==pm] & 
                        W[1==pm]<Inf))
        stopifnot(all(c(n,n) == dim(W)))

        stopifnot(all(is.na(d[producers]))) 
        stopifnot(all(is.na(d[externals]))) 
        stopifnot(all(0<=d[consumers] & d[consumers]<=Inf))
        stopifnot(all(1==length(q)))
        stopifnot(0<=q & q<=Inf)
        stopifnot(all(0<K & K<Inf))

        stopifnot(all(is.na(a[consumers,consumers])))
        stopifnot(all(is.na(a[externals,externals])))
        stopifnot(all(0<a[producers,producers] & 
                        a[producers,producers]<Inf))
        stopifnot(all(c(n,n) == dim(a)))

        stopifnot(all(is.na(V[producers])))
        stopifnot(all(is.na(V[consumers])))
        stopifnot(all(0<=V[externals] & V[externals]<Inf))

        stopifnot(all(is.na(Z[producers])))
        stopifnot(all(is.na(Z[consumers])))
        stopifnot(all(0<=Z[externals] & Z[externals]<Inf))
    })

    with(c(community, params), 
    {
        M <- NP(community, 'M')

        m.term <- ((M[m]/M)^exponent)
        rho <- ((fr * ar)/(fr[m] * ar[m])) * m.term
        x <- ((fT * aT)/(fr[m] * ar[m])) * m.term
        y <- fJ * 
             matrix(aJ, nrow=length(M), ncol=length(M), byrow=TRUE) / 
             matrix(fT*aT, nrow=length(M), ncol=length(M), byrow=TRUE)

        # Good to cache vectors of producers and consumers for performance.

        # R is 1-indexed
        producers <- NodeNameIndices(community, .Producers(community))
        consumers <- NodeNameIndices(community, .Consumers(community))
        externals <- NodeNameIndices(community, .Externals(community))

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

        if(length(externals)>0)
        {
            externals.c <- as.integer(externals-1)
        }
        else
        {
            externals.c <- 0    # Can't pass value of NULL as arg to .C()
        }

        # A multiplier to convert units of t' (model equations) to years
        # Homage to Yodzis and Innes, 2.13 (p 42):
        one.t.prime <- 1/as.numeric(fr[m]*ar[m]*(M[m]^-exponent))

        v <- V * one.t.prime
        z <- Z * one.t.prime

        params <- list(community=community, 
                       n.species=NumberOfNodes(community), 
                       producers=producers, 
                       producers.c=producers.c, 
                       n.producers=length(producers),
                       consumers=consumers, 
                       consumers.c=consumers.c, 
                       n.consumers=length(consumers),
                       externals=externals, 
                       externals.c=externals.c, 
                       n.externals=length(externals),
                       rho=rho, x=x, y=y, e=e, fe=fe, 
                       W=W, d=d, q=q, # Functional response
                       K=K, a=a,      # Growth model
                       v=v, z=z,      # Inputs and decays
                       # The following are not required by model equations but 
                       # are useful for subsequent analyses.
                       m=m, one.t.prime=one.t.prime)

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

         stopifnot(all(is.na(W[0==pm])))
         stopifnot(all(0<W[1==pm] & 
                         W[1==pm]<Inf))

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

