.PlotDeviations <- function(community, 
                            y.0, 
                            y.now, 
                            xlim=NULL, 
                            ylim=NULL, 
                            colour.by,
                            colour.spec,
                            col=NULL, 
                            symbol.by,
                            symbol.spec,
                            pch.0=NULL, 
                            bg.by,
                            bg.spec,
                            bg=NULL, 
                            cex.by=NULL,
                            cex.spec=NULL,
                            cex=NULL, 
                            label.colour.by=NULL,
                            label.colour.spec=NULL,
                            label.colour=NULL, 
                            highlight.links=NULL,
                            highlight.nodes=Cannibals,
                            show.nodes.as='points', 
                            node.labels=NULL, 
                            label.cex=0.6, 
                            line.col=DefaultLinkColour(), 
                            lwd=1,
                            pch.now=4, 
                            ...)
{
    stopifnot(length(y.0)==NumberOfNodes(community))
    stopifnot(length(y.now)==NumberOfNodes(community))

    n <- NumberOfNodes(community)

    log10M <- Log10M(community)
    log10y.0 <- log10(y.0)
    log10y.0[is.infinite(log10y.0)] <- NA
    log10y.now <- log10(y.now)
    log10y.now[is.infinite(log10y.now)] <- NA

    points <- PlaceMissingPoints(rep(log10M, 2), xlim, c(log10y.0, log10y.now), 
                                 ylim)

    plot(points[,1], points[,2], type='n', xlim=xlim, ylim=ylim, ...)
    cheddar:::.AddAxisTicks(...)

    # Draw lines between log10y.0 and log10y.now
    px <- points[,1]
    py <- points[,2]

    lx <- rep(NA, 3*n)
    lx[seq(1, length(lx), by=3)] <- px[1:n]
    lx[seq(2, length(lx), by=3)] <- px[1:n]

    ly <- rep(NA, 3*n)
    ly[seq(1, length(ly), by=3)] <- py[1:n]
    ly[seq(2, length(ly), by=3)] <- py[(1+n):length(py)]

    # No lines for extinct species
    extinct <- which(is.na(log10y.now))
    ly[1+3*(extinct-1)] <- NA
    lines(lx, ly, col=line.col, lwd=lwd)

    # Plot y.0
    cheddar:::.PlotNodes(community=community, 
                         x=points[1:n,1], 
                         y=points[1:n,2], 
                         colour.by=colour.by, colour.spec=colour.spec, col=col, 
                         symbol.by=symbol.by, symbol.spec=symbol.spec,pch=pch.0,
                         bg.by=bg.by, bg.spec=bg.spec, bg=bg, 
                         cex.by=cex.by, cex.spec=cex.spec, cex=cex, 
                         label.colour.by=label.colour.by, 
                         label.colour.spec=label.colour.spec, 
                         label.colour=label.colour, 
                         highlight.nodes=highlight.nodes, 
                         lowlight.nodes=0==y.0, 
                         show.nodes.as=show.nodes.as, 
                         node.labels=node.labels, label.cex=label.cex, 
                         ...)

    # Plot y.now

    cheddar:::.PlotNodes(community=community, 
                         x=points[(1:n)+n,1], 
                         y=points[(1:n)+n,2], 
                         colour.by=colour.by, colour.spec=colour.spec, col=col, 
                         symbol.by=NULL, symbol.spec=NULL, pch=pch.now,
                         bg.by=bg.by, bg.spec=bg.spec, bg=bg, 
                         cex.by=cex.by, cex.spec=cex.spec, cex=cex, 
                         label.colour.by=NULL, 
                         label.colour.spec=NULL, 
                         label.colour=NULL, 
                         highlight.nodes=NULL, 
                         lowlight.nodes=0==y.now, 
                         show.nodes.as='points', 
                         node.labels=NULL, label.cex=NULL, 
                         ...)
}

PlotNDeviations <- function(community, N.now, 
                            xlab=Log10MLabel(community), 
                            main=CPS(community)$title, 
                            ylab=Log10NLabel(community), 
                            ...)
{
    # Plot log10(N) and log10(N.now) against log10(M)
    cheddar:::.RequireM(community)
    cheddar:::.RequireN(community)

    .PlotDeviations(community=community, y.0=NP(community, 'N'), N.now, 
                    xlab=xlab, ylab=ylab, main=main, ...)
}

PlotBDeviations <- function(community, B.now, 
                            xlab=Log10MLabel(community), 
                            main=CPS(community)$title, 
                            ylab=Log10BLabel(community), 
                            ...)
{
    # Plot log10(N) and log10(N.now) against log10(M)
    cheddar:::.RequireM(community)
    cheddar:::.RequireN(community)

    .PlotDeviations(community=community, y.0=Biomass(community), B.now, 
                    xlab=xlab, ylab=ylab, main=main, ...)
}

.PlotYvT <- function(community, tseries, category, from, from.time, 
                     to, to.time, xlab, xlim, ylab, ylim, show.legend, 
                     col, lty, divide.Y.by.M, plot.species.labels, 
                     cex=1, ...)
{
    # TODO colour.spec scheme as used by cheddar
    # A private helper that plot timeseries values against time
    # If divide.Y.by.M is TRUE, the tseries data are divided by community$M 
    # before plotting. Only the subset of tseries to be plotted is divided 
    # by M, giving a performance gain when plotting a small subset of a large 
    # tseries.

    if(!is.null(from.time))
    {
        from <- which(tseries[,1,drop=FALSE]==from.time)
    }

    if(!is.null(to.time))
    {
        to <- which(tseries[,1,drop=FALSE]==to.time)
    }

    stopifnot(all.equal(1, length(to), length(from)))

    time <- tseries[from:to,1,drop=FALSE]
    tseries <- tseries[from:to,-1,drop=FALSE]

    if(divide.Y.by.M)
    {
        tseries <- t(apply(tseries,1,'/',NP(community, 'M')))
    }
    y <- sapply(unique(category), 
                function(c) apply(tseries[,category==c,drop=FALSE], 1, sum))
    colnames(y) <- unique(category)

    if(missing(col))
    {
        default <- DefaultCategoryColours()
        if(all(colnames(y) %in% names(default)))
        {
            col <- default[colnames(y)]
        }
        else
        {
            col <- 1
        }
    }

    if(FALSE && !is.null(names(col)))
    {
        col <- col[colnames(y)]
    }

    matplot(time, log10(y), type="l", xlab=xlab, xlim=xlim, ylab=ylab, 
            ylim=ylim, col=col, lty=lty, cex=cex, ...)

    if(FALSE)
    {
        if(plot.species.labels)
        {
            # Plot species labels where extinctions occur
            extinct <- apply(y, 2, 
                             function(col)
                             {
                                 x <- which(0!=col)
                                 if(length(col)==length(x)) return(NA) 
                                 else return(tail(x,1))
                             })

            # Intersection of extinct species and what the user wants to see
            # Perform the lookup by name
            extinct <- extinct[!is.na(extinct)]
            extinct <- extinct[names(extinct) %in% species]

            sapply(names(extinct), 
                   function(e) text(tseries[extinct[e],'time'], 
                                    log10(y[extinct[e],e]), 
                                    species.labels[e], 
                                    cex=cex))
        }

        if(plot.species.labels)
        {
            # Plot species labels for extant species
            extant <- NP(community, 'node')[which(tseries[to, -1]>0)]

            # Intersection of extinct species and what the user wants to see
            # Perform the lookup by name
            extant <- extant[extant %in% species]

            if(length(extant)>0)
            {
                text(x=tseries[to, 'time'], y=log10(tail(y, 1)[,extant]), 
                     labels=species.labels[extant], pos=4, cex=cex, ...)
            }
        }
    }

    if(FALSE && show.legend)
    {
        legend("topright", legend=colnames(y), lty=lty, col=col, cex=cex)
    }

    # Tick marks at top and bottom of plot
    axis(3, labels=FALSE)
    axis(4, labels=FALSE)
}

PlotNvT <- function(community, tseries, 
          category=relevel(factor(NP(community, 'category')), 'producer'),
          from=1, from.time=NULL, 
          to=nrow(tseries), to.time=NULL, 
          main=CPS(community)$title, 
          xlab="time (t')", xlim=NULL, 
          ylab=Log10NLabel(community), ylim=NULL, 
          col, 
          show.legend=FALSE, 
          lty=1, ...)
{
    .PlotYvT(community, 
             tseries, 
             category, 
             from=from, from.time=from.time, 
             to=to, to.time=to.time, 
             xlab=xlab, xlim=xlim, ylab=ylab, ylim=ylim, 
             show.legend=show.legend, col=col, lty=lty, 
             divide.Y.by.M=TRUE, main=main, ...)
}

PlotBvT <- function(community, tseries, 
          category=relevel(factor(NP(community, 'category')), 'producer'),
          from=1, from.time=NULL, 
          to=nrow(tseries), to.time=NULL, 
          main=CPS(community)$title, 
          xlab="time (t')", xlim=NULL, 
          ylab=Log10BLabel(community), ylim=NULL, 
          col, 
          show.legend=FALSE, 
          lty=1, ...)
{
    .PlotYvT(community, 
             tseries, 
             category, 
             from=from, from.time=from.time, 
             to=to, to.time=to.time, 
             xlab=xlab, xlim=xlim, ylab=ylab, ylim=ylim, 
             show.legend=show.legend, col=col, lty=lty, 
             divide.Y.by.M=FALSE, main=main, ...)
}

FormatSSEModel <- function(sse.model)
{
    # A string representation for SSE model
    return (substitute(italic(SSE[model]) == s, 
                       list(s=sprintf("%.2f",sse.model))))
}

FormatSSERegression <- function(community)
{
    return (substitute(italic(SSE[regression]) == sse, 
                       list(sse=sprintf("%.2f", SSERegression(community)))))
}

PlotFluxBars <- function(community, flux, time, ylim=NULL, 
                         ylab=as.formula(paste('log[10](italic(B)) ~ (', 
                                               CP(community, 'M.units'), '~', 
                                               CP(community, 'N.units'), 
                                               '~ year^-1)',
                                               sep='')),
                         col=c(growth='green3', respiration='purple3', 
                               consumption='blue3', assimilation='red3'))
{
    # Plots log10-transformed flux per species at the given time index
    growth <- log10(flux$growth[time,])
    respiration <- log10(flux$respiration[time,])
    total.consumption <- log10(flux$total.consumption[time,])
    total.assimilation <- log10(flux$total.assimilation[time,])

    growth[is.infinite(growth)] <- NA
    respiration[is.infinite(respiration)] <- NA
    total.consumption[is.infinite(total.consumption)] <- NA
    total.assimilation[is.infinite(total.assimilation)] <- NA

    if(is.null(ylim))
    {
        ylim <- range(c(0, growth, respiration, 
                        total.consumption, 
                        total.assimilation), na.rm=TRUE)
    }

    plot(0, 0, xlim=c(0, NumberOfNodes(community)), ylim=ylim, ylab=ylab, 
         type='n', main=paste('Flux at', flux$time[time]))

    producers <- Producers(community)
    segments(producers-1, growth[producers], producers, 
             growth[producers], col=col['growth'])

    segments( (1:NumberOfNodes(community))-1, total.consumption, 
             1:NumberOfNodes(community), total.consumption, 
             col=col['consumption'])

    consumers <- Consumers(community)
    segments(consumers-1, respiration[consumers], consumers, 
             respiration[consumers], col=col['respiration'])

    segments(consumers-1, total.assimilation[consumers], consumers, 
             total.assimilation[consumers], col=col['assimilation'])
}

PlotFlux <- function(community, flux, species, xlab="time (t')", lty=1, 
                     col=c(growth='green3', respiration='purple3', 
                                  consumption='blue3', assimilation='red3'), 
                     plot.log10=FALSE)
{
    # Plots absolute flux for a single species
    species <- cheddar:::.ResolveToNodeIndices(community, species)

    stopifnot(1==length(species))

    if(FALSE)
    {
        flux$respiration <- -flux$respiration
        flux$consumption <- -flux$consumption
        flux$total.consumption <- -flux$total.consumption
    }

    if(plot.log10)
    {
        ylab <- as.formula(paste('log[10](italic(B)) ~ (', 
                                 CP(community,'M.units'), '~', 
                                 CP(community,'N.units'), '~ year^-1)',
                                 sep=''))

        growth <- log10(flux$growth[,species])
        respiration <- log10(flux$respiration[,species])
        consumption <- log10(flux$consumption[,species,])
        total.consumption <- log10(flux$total.consumption[,species])
        assimilation <- log10(flux$assimilation[,,species])
        total.assimilation <- log10(flux$total.assimilation[,species])

        growth[is.infinite(growth)] <- NA
        respiration[is.infinite(respiration)] <- NA
        consumption[is.infinite(consumption)] <- NA
        total.consumption[is.infinite(total.consumption)] <- NA
        assimilation[is.infinite(assimilation)] <- NA
        total.assimilation[is.infinite(total.assimilation)] <- NA
        total <- rep(NA, length(growth))
    }
    else
    {
        ylab <- as.formula(paste('B ~ (', 
                                 CP(community, 'M.units'), '~', 
                                 CP(community, 'N.units'), '~ year^-1)',
                                 sep=''))

        growth <- flux$growth[,species]
        respiration <- flux$respiration[,species]
        consumption <- flux$consumption[,species,]
        total.consumption <- flux$total.consumption[,species]
        assimilation <- flux$assimilation[,,species]
        total.assimilation <- flux$total.assimilation[,species]
        total <- flux$total[,species]
    }

    .PlotFluxLines(community, species, flux$time, growth, respiration, 
                   consumption, total.consumption, assimilation, 
                   total.assimilation, total, col, lty, ylab)
    .PlotFluxLabels(community, species, flux$time, growth, respiration, 
                    consumption, total.consumption, assimilation, 
                    total.assimilation, col, lty)
}

.PlotFluxLines <- function(community, species, time, growth, respiration, 
                           consumption, total.consumption, assimilation, 
                           total.assimilation, total, col, lty, ylab)
{
    # na.rm=TRUE because some entries will be NA, e.g. respiration for a 
    # producer
    ylim <- range(c(growth, respiration, consumption, total.consumption, 
                    assimilation, total.assimilation, total), na.rm=TRUE)

    plot(0, 0, xlim=c(time[1], tail(time,1)), ylim=ylim, type="n", 
         xlab="time (t')", ylab=ylab, 
         main=paste("Absolute biomass fluxes for ", 
                    NP(community, 'node')[species], ' [', species, ']', sep=''))

    # A zero flux line
    abline(h=0, col="grey")

    # Need to be careful here because species might not have consumers
    # and consumers might not have any resource species
    if(IsProducer(community)[species])
    {
        lines(time, growth, col=col['growth'], lty=lty)
    }
    else if(IsConsumer(community)[species])
    {
        lines(time, respiration, col=col['respiration'], lty=lty)
        resources <- ResourcesByNode(community)[[species]]
        if(length(resources)>0)
        {
            matlines(time, assimilation[,resources], col=col['assimilation'], 
                     lty=lty)
            if(length(resources)>1)
            {
                lines(time, total.assimilation, col=col['assimilation'], 
                      lty=lty)
            }
        }
    }

    consumers <- ConsumersByNode(community)[[species]]
    if(length(consumers)>0)
    {
        matlines(time, consumption[,consumers], 
                 col=col['consumption'], lty=lty)
        if(length(consumers)>1)
        {
            lines(time, total.consumption, col=col['consumption'], 
                  lty=lty)
        }
    }

    if(!is.null(total))
    {
        lines(time, total)
    }
}

.PlotFluxLabels <- function(community, species, time, growth, respiration, 
                            consumption, total.consumption, assimilation, 
                            total.assimilation, col, lty)
{
    # Place labels next to each flux line
    # Labels are positioned either at left or right of each trace 
    # depending on where the space is
    first <- 1
    last <- length(time)
    range.start <- range(c(growth[first], 
                           respiration[first], 
                           total.consumption[first], 
                           total.assimilation[first]), na.rm=TRUE)
    range.end <- range(c(growth[last], 
                         respiration[last], 
                         total.consumption[last], 
                         total.assimilation[last]), na.rm=TRUE)

    if(diff(range.start) > diff(range.end))
    {
        text.i <- first
        text.pos <- 2
    }
    else
    {
        text.i <- last
        text.pos <- 4
    }

    if(IsProducer(community)[species])
    {
        text(x=time[text.i], y=growth[text.i], labels="G+", 
             cex=0.8, col=col['growth'], pos=text.pos)
    }
    else if(IsConsumer(community)[species])
    {
        text(x=time[text.i], y=respiration[text.i], 
             labels="R-", cex=0.8, col=col['respiration'], pos=text.pos)
        resources <- ResourcesByNode(community)[[species]]
        text(x=rep(time[text.i],length(resources)), 
             y=assimilation[text.i, resources], 
             labels=paste(ResourcesByNode(community)[[species]], '+'), 
             cex=0.8, col=col['assimilation'], pos=text.pos)
        if(length(resources)>1)
        {
            text(x=time[text.i], 
                 y=total.assimilation[text.i], 
                 labels="T+", cex=0.8, pos=text.pos, col=col['assimilation'])
        }
    }

    consumers <- ConsumersByNode(community)[[species]]
    if(length(consumers)>0)
    {
        text(x=rep(time[text.i],length(consumers)), 
             y=consumption[text.i, consumers], 
             labels=paste(ConsumersByNode(community)[[species]], '-'), 
             cex=0.8, col=col['consumption'], pos=text.pos)

        if(length(consumers)>1)
        {
            text(x=time[text.i], 
                 y=total.consumption[text.i], 
                 labels="T-", cex=0.8, pos=text.pos, col=col['consumption'])
        }
    }
}

