######################################################################
### Script to plot covid pileups

### 0. Some functions:
plotMutation <- function(mutation=mutation, ylim=ylim1){
    polygon(x=rep(mutation,each=2),
            y=c(-ylim[2]/100,-ylim[2]/30,-ylim[2]/30,-ylim[2]/100),
            col="salmon")
}

plotSequence <- function(sequence, ylim=ylim, cex){
    text(x=1:length(sequence), y=-ylim[2]/50, labels=toupper(sequence),
         col=c(rep('black', 607), rep('red', 9), rep('black', 600)),
         cex=cex)
}

plotCoverage <- function(xAcov, xPileup, readLength, main, xlim,
                         sequence=FALSE, scex=1, ylim=NULL){
    if(length(readLength) > 1){
        xTot <- rowSums(xPileup[,as.character(readLength)])
        xA <- rowSums(xAcov[,as.character(readLength)])
    }else{
        xTot <- xPileup[,as.character(readLength)]
        xA <- xAcov[,as.character(readLength)]
    }
    if(length(ylim) != 2)
        ylim <- c(0, max(xTot))
    cmax <- dim(xPileup)[1]
    par(mar=c(5,5,2,2))
    plot(1,0, pch=".", xlim=xlim, ylim=ylim, xlab="", frame.plot=F,
         main=main, ylab="Read coverage", cex.axis=2,
         cex.main=3, cex.lab=2)
    polygon(x=sort(c(0.5, seq(0.5, cmax, 1), seq(1.5, cmax, 1), cmax)),
            y=c(0, xTot[sort(c(1:cmax, 1:(cmax-1)))], 0),
            col="lightgray")
    polygon(x=sort(c(seq(0.5, cmax, 1), seq(1.5, cmax, 1))),
            y=c(xA[sort(c(1:cmax, 1:(cmax-1)))]),
            col="skyblue", border="black")
    if (length(sequence) > 1)
        plotSequence(sequence=sequence, ylim=ylim, cex=scex)
}

getFrame <- function(x, readLength, cds=cds, offset, region=c(0, length(x))){
    readFrame <- matrix(ncol=3)
    x.length <- x[,3] - x[,2]
    for(I in readLength){
        rData <- x[x.length == I, 2] + 1 + offset[as.character(I)]
        ind <- which(rData >= region[1] & rData <= region[2])
        if(length(ind) >=1){
            rData <- rData[which(rData >= region[1] & rData <= region[2])]
            d <- cbind(table((rData - cds) %% 3), c(0,1,-1), I)
        }else{
            d <- cbind(c(NA,NA,NA), c(0,1,-1), I)
        }
        readFrame <- rbind(readFrame, d)
    }
    return(readFrame[-1,])
}

getCoverage <- function(x, readLength, offset, cLength, ref=5){
    x.depth <- matrix(0, ncol=length(readLength),
                      nrow=cLength)
    colnames(x.depth) <- readLength
    x.length <- x[,3] - x[,2]
    for(I in readLength){
        ind <- which(x.length == I)
        x.reg <- x[ind,]
        if(ref == 5)
            wt.depth <- table(x.reg[,2]+ 1 + offset[as.character(I)])
        else if(ref == 3)
            wt.depth <- table(x.reg[,3] - offset[as.character(I)])
        ind1 <- as.numeric(names(wt.depth))
        x.depth[ind1, as.character(I)] <- wt.depth
    }
    return(x.depth)
}

extendCoverage <- function(x, readLength){
    for(I in 1:length(readLength)){
        idx <- which(x[,I] > 0)
        idx.ext <- vector()
        for(K in idx){
            idx.ext <- seq(K, K + readLength[I] - 1, by=1)
            x[idx.ext, I] <- x[idx.ext, I] + x[K,I]
        }
    }
    return(x)
}

extendRead <- function(x, l){
    return(seq(x, x+l, 1))
}

getReadCoverage <- function(x, readLength, cLength){
    x.depth <- matrix(0, ncol=length(readLength),
                      nrow=cLength)
    colnames(x.depth) <- readLength
    x.length <- x[,3] - x[,2]
    for(I in readLength){
        ind <- which(x.length == I)
        x.reg <- x[ind,]
        ind1 <- unlist(tapply(x.reg[,2], ind, extendRead, I))
        wt.depth <- rep(0, cLength)
        for(K in ind1){wt.depth[K] <- wt.depth[K] + 1}
        x.depth[, as.character(I)] <- wt.depth
    }
    return(x.depth)
}

readCumSum <- function(bed, cLength){
    cSum <- rep(0, cLength)
    for(I in 2:cLength){
        reads <- length(which(bed[,2] == I))
        cSum[I] <- sum(cSum[(I-1)], reads)
    }
    return(cSum)
}
