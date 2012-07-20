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
                            main=CPS(community)$title, 
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

    plot(points[,1], points[,2], type='n', xlim=xlim, ylim=ylim, main=main, ...)
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
                            ylab=Log10NLabel(community), 
                            ...)
{
    # Plot log10(N) and log10(N.now) against log10(M)
    cheddar:::.RequireM(community)
    cheddar:::.RequireN(community)

    .PlotDeviations(community=community, y.0=NP(community, 'N'), N.now, 
                    xlab=xlab, ylab=ylab, ...)
}

PlotBDeviations <- function(community, B.now, 
                            xlab=Log10MLabel(community), 
                            ylab=Log10BLabel(community), 
                            ...)
{
    # Plot log10(N) and log10(N.now) against log10(M)
    cheddar:::.RequireM(community)
    cheddar:::.RequireN(community)

    .PlotDeviations(community=community, y.0=Biomass(community), B.now, 
                    xlab=xlab, ylab=ylab, ...)
}

.PlotYvT <- function(community, tseries, 
                     category, col, 
                     from=1, from.time=NULL, 
                     to=nrow(tseries), to.time=NULL, 
                     main=CPS(community)$title, 
                     xlab="time (t')", xlim=NULL, 
                     ylab=NULL, ylim=NULL, 
                     show.legend=FALSE, 
                     lty=1, 
                     divide.Y.by.M, 
                     ...)
{
    # A private helper that plot timeseries values against time

    # If divide.Y.by.M is TRUE, the tseries data are divided by community$M 
    # before plotting. Only the subset of tseries to be plotted is divided 
    # by M, giving a performance gain when plotting a small subset of a large 
    # tseries.

    if(missing(category))
    {
        category <- NP(community, 'node')
    }

    if(missing(col))
    {
        col <- DefaultCategoryColours()[NP(community, 'category')]
    }

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

    if(all(category==colnames(tseries)))
    {
        y <- tseries
    }
    else
    {
        y <- sapply(unique(category), 
                    function(c) apply(tseries[,category==c,drop=FALSE], 1, sum))
        colnames(y) <- unique(category)
    }

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

    matplot(time, log10(y), type="l", xlab=xlab, xlim=xlim, ylab=ylab, 
            ylim=ylim, col=col, lty=lty, main=main, ...)

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
                                    species.labels[e]))
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
                     labels=species.labels[extant], pos=4, ...)
            }
        }
    }

    if(FALSE && show.legend)
    {
        legend("topright", legend=colnames(y), lty=lty, col=col)
    }

    # Tick marks at top and bottom of plot
    axis(3, labels=FALSE)
    axis(4, labels=FALSE)
}

PlotNvT <- function(community, tseries, ylab=Log10NLabel(community), ...)
{
    .PlotYvT(community, tseries, divide.Y.by.M=TRUE, ylab=ylab, ...)
}

PlotBvT <- function(community, tseries, ylab=Log10BLabel(community), ...)
{
    .PlotYvT(community, tseries, divide.Y.by.M=FALSE, ylab=ylab, ...)
}

