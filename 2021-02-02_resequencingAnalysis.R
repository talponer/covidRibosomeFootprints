######################################################################
## This script contains instruction to analyse the newly sequenced
## PF14 and PD14, the samples already analysed on 2021-01-28 but with
## higher depth in the hope that they show a clearer signal in the
## fram-shifting region. There are two run for each library (1 and 2)

######################################################################
### Script to plot covid pileups
library(ggplot2)
library(seqinr)
source('covidFunctions.R')
library(RColorBrewer)
library(zoo)
library(scales)

redCols <- gradient_n_pal(brewer.pal(9, 'Reds'))(seq(0, 1, length.out = 50))
blueCols <- gradient_n_pal(brewer.pal(9, 'Blues'))(seq(0, 1, length.out = 50))

######################################################################
######################################################################
## Some stats:

## Trimming stats
trimStatA <- read.table('../2021-01-27_trimStat', header=T, sep="\t")
trimStat <- matrix(0, ncol=dim(trimStatA)[2]-1,
                   nrow=length(unique(trimStatA$Sample)))
rownames(trimStat) <- unique(trimStatA$Sample)
colnames(trimStat) <- colnames(trimStatA)[-1]
for (I in colnames(trimStat)) trimStat[,I] <- tapply(trimStatA[,I],
                                                     trimStatA$Sample,
                                                     sum)
trimStat <- as.data.frame(trimStat)
noAdapter <- trimStat$Total - trimStat$With.adapter
cutadaptFiltered <- trimStat$With.adapter - trimStat$written
trimStatNew <- data.frame(no.adapter = noAdapter,
                          cutadapt.filtered = cutadaptFiltered,
                          trimStat[,4:6])

pdf("figures/2021-01-27/trimmingStats.pdf", width=5)
par(mar=c(5,4,4,8))
barplot(t(as.matrix(trimStatNew)),
        col=rainbow(5, start=0.1, end=0.9, alpha=0.6), main='Filtering')
par(xpd=TRUE, mar=c(5,4,4,1))
legend(dim(newMapStats)[1]*1.3, max(trimStat[,1]), legend=colnames(trimStatNew), col=rainbow(5, start=0.1, end=0.9, alpha=0.6), pch=0, pt.lwd=8, bty="n")
dev.off()

## Mapping stats
newMapStatsA <- read.table('../2021-01-27_mappingStats', header=T, sep="\t")
newMapStats <- matrix(0, ncol=dim(newMapStatsA)[2]-1,
                   nrow=length(unique(newMapStatsA$Sample)))
rownames(newMapStats) <- unique(newMapStatsA$Sample)
colnames(newMapStats) <- colnames(newMapStatsA)[-1]
for (I in colnames(newMapStats)) newMapStats[,I] <- tapply(newMapStatsA[,I],
                                                     newMapStatsA$Sample,
                                                     sum)
newMapStats <- as.data.frame(newMapStats)
newMapStatPer <- newMapStats/newMapStats[,1]

pdf("figures/2021-01-27/mappingStatsNew.pdf", width=5)
par(mar=c(5,4,4,8))
barplot(t(as.matrix(newMapStats[,-1])),
        col=rainbow(4, start=0.1, end=0.9, alpha=0.6), main='Mapping')
par(xpd=TRUE, mar=c(5,4,4,1))
legend(dim(newMapStats)[1]*1.3, max(newMapStats[,1]), legend=colnames(newMapStats)[-1], col=rainbow(4, start=0.1, end=0.9, alpha=0.6), pch=0, pt.lwd=8, bty="n")
dev.off()

pdf("figures/2021-01-27/mappingStatsNewPercentage.pdf", width=5)
par(mar=c(5,4,4,8))
barplot(t(as.matrix(newMapStatPer[,-1])),
        col=rainbow(4, start=0.1, end=0.9, alpha=0.6), main='Mapping')
par(xpd=TRUE, mar=c(5,4,4,1))
legend(dim(newMapStats)[1]*1.3, 1, legend=colnames(newMapStats)[-1], col=rainbow(4, start=0.1, end=0.9, alpha=0.6), pch=0, pt.lwd=8, bty="n")
dev.off()


######################################################################
### Get the sample annotation:
samplesAnno <- read.table('../2021-01-27_samples.anno', header=T,
                          sep='\t')
libraries <- c('PF14', 'PD14')
samplesAnno <- samplesAnno[which(samplesAnno$ID %in% libraries),]

dSamples <- c('PD14.AGCTA', 'PD14.TCTAG')
mSamples <- c('PF14.AGCTA', 'PF14.TCTAG')
names(dSamples) <- names(mSamples) <- c('WildType', 'Mutant')

dSamplesCX <- c('PD14.ATCGT', 'PD14.CGTAA')
mSamplesCX <- c('PF14.ATCGT', 'PF14.CGTAA')
names(dSamplesCX) <- names(mSamplesCX) <- c('WildType', 'Mutant')

######################################################################
### Read length distribution:
readLength <- list()
for(I in libraries){
    for(K in unique(samplesAnno$Adapter)){
        fname <- paste("../trimmed_data/", I, "_354_1_001.", K,
                       ".readLengthDist", sep='')
        name <- paste(I, K, sep='.')
        readLength[[name]] <- read.table(fname, header=F)
    }
}


pdf('figures/2021-02-02/disomeSamplesReadLength.pdf', height=15, width=10)
par(mfrow=c(4,1))
idx <- grep('PD14', samplesAnno$ID)
for(I in unique(samplesAnno$ID[idx])){
    for(K in unique(samplesAnno$Adapter)){
        name <- paste(I, K, sep='.')
        barplot(readLength[[name]][,2],
                names.arg=readLength[[name]][,1], main=name)
    }
}
dev.off()

pdf('figures/2021-02-02/monosomesSamplesReadLength.pdf', height=20, width=5)
par(mfrow=c(4,1))
idx <- grep('PF14', samplesAnno$ID)
for(I in unique(samplesAnno$ID[idx])){
    for(K in unique(samplesAnno$Adapter)){
        name <- paste(I, K, sep='.')
        barplot(readLength[[name]][,2],
                names.arg=readLength[[name]][,1], main=name)
    }
}
dev.off()

######################################################################
######################################################################
## Calculate the best offset for each read length (to get majority of
## reads at frame 0 in a region downsteam the CDS start and before
## the frameshifting site)

covidSeq <- read.fasta('../annotations/sarsCov2.fa')

cds <- 116
mutation <- c(608, 617)

cLength <- read.table('../annotations/sarsCov2.lenght', header=F,
                      row.names=1)

## get the bed files
bed <- list()
for(I in unique(samplesAnno$ID)){
    for(K in unique(samplesAnno$Adapter)){
        fname <- paste("../mapping_data/sam/", I, ".", K, ".bed",
                       sep='')
        name <- paste(I, K, sep='.')
        tmp <- read.table(fname, header=F)
        type <- samplesAnno$Genotype[which(samplesAnno$ID == I & samplesAnno$Adapter == K)]
        bed[[name]] <- tmp[tmp[,1] == type,]
    }
}

str(bed)

lDepth <- lapply(bed, dim)

######################################################################
######################################################################
### Plot read coverage stratified by read length
bestOffset <- rep(14, length(26:100))
names(bestOffset) <- 26:100

## Monosomes
wtMCoverage <- getCoverage(bed[[mSamples['WildType']]],
                           readLength=26:35, offset=bestOffset,
                           cLength=cLength[1,1])
mtMCoverage <- getCoverage(bed[[mSamples['Mutant']]],
                           readLength=26:35, offset=bestOffset,
                           cLength=cLength[1,1])

wt.m.read.coverage <- getReadCoverage(bed[[mSamples['WildType']]],
                                      readLength=26:35,
                                      cLength=cLength[1,1])

mt.m.read.coverage <- getReadCoverage(bed[[mSamples['Mutant']]],
                                      readLength=26:35,
                                      cLength=cLength[1,1])

### Write out the wig fies for submission
write.table(rowSums(wt.m.read.coverage),
            file="GEOsubmission/WildType_monosome.wig", row.names=F,
            col.names=F)

write.table(rowSums(mt.m.read.coverage),
            file="GEOsubmission/Mutant_monosome.wig", row.names=F,
            col.names=F)


## Calculate read per 10000:
mil <- 10000
wtMCoverage <- mil * wtMCoverage / lDepth[[mSamples['WildType']]][1]
mtMCoverage <- mil * mtMCoverage / lDepth[[mSamples['Mutant']]][1]
wt.m.read.coverage <- mil * wt.m.read.coverage / lDepth[[mSamples['WildType']]][1]
mt.m.read.coverage <- mil * mt.m.read.coverage / lDepth[[mSamples['Mutant']]][1]

pdf("figures/2021-02-02/WT-Mutant_Monosomes_readCoverage.pdf", width=30,
    height=15)
par(mfrow=c(2,1))
plotCoverage(xAcov=wtMCoverage, xPileup=wt.m.read.coverage,
             readLength=26:35, main='Wild Type', xlim=c(525,675),
             sequence=covidSeq$WildType, scex=1)
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
plotCoverage(xAcov=mtMCoverage, xPileup=mt.m.read.coverage,
             readLength=26:35, main='Mutant', xlim=c(525,675),
             sequence=covidSeq$Mutant, scex=1)
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
dev.off()

## Disomes:
wtDCoverage <- getCoverage(bed[[dSamples['WildType']]],
                           readLength=36:100, offset=bestOffset,
                           cLength=cLength[1,1])
mtDCoverage <- getCoverage(bed[[dSamples['Mutant']]],
                           readLength=36:100, offset=bestOffset,
                           cLength=cLength[1,1])


wt.d.read.coverage <- getReadCoverage(bed[[dSamples['WildType']]],
                                      readLength=36:100,
                                      cLength=cLength[1,1])

mt.d.read.coverage <- getReadCoverage(bed[[dSamples['Mutant']]],
                                      readLength=36:100,
                                      cLength=cLength[1,1])

### Write out the wig fies for submission
write.table(rowSums(wt.d.read.coverage),
            file="GEOsubmission/WildType_disomes.wig", row.names=F,
            col.names=F)

write.table(rowSums(mt.d.read.coverage),
            file="GEOsubmission/Mutant_disomes.wig", row.names=F,
            col.names=F)

## Calculate read per 1000:
wtDCoverage <- mil * wtDCoverage / lDepth[[mSamples['WildType']]][1]
mtDCoverage <- mil * mtDCoverage / lDepth[[mSamples['Mutant']]][1]
wt.d.read.coverage <- mil * wt.d.read.coverage / lDepth[[mSamples['WildType']]][1]
mt.d.read.coverage <- mil * mt.d.read.coverage / lDepth[[mSamples['Mutant']]][1]

pdf("figures/2021-02-02/WT-Mutant_disome_readCoverage.pdf", width=30,
    height=15)
par(mfrow=c(2,1))
plotCoverage(xAcov=wtDCoverage, xPileup=wt.d.read.coverage,
             readLength=36:100, main='Wild Type', xlim=c(525,675),
             sequence=covidSeq$WildType, scex=1)
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
plotCoverage(xAcov=mtDCoverage, xPileup=mt.d.read.coverage,
             readLength=36:100, main='Mutant', xlim=c(525,675),
             sequence=covidSeq$Mutant, scex=1)
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
dev.off()

## Restrict the disome to fragments from 61 to 70
wtDCoverage1 <- getCoverage(bed[[dSamples['WildType']]],
                           readLength=61:70, offset=bestOffset,
                           cLength=cLength[1,1])
mtDCoverage1 <- getCoverage(bed[[dSamples['Mutant']]],
                           readLength=61:70, offset=bestOffset,
                           cLength=cLength[1,1])


wt.d.read.coverage1 <- getReadCoverage(bed[[dSamples['WildType']]],
                                      readLength=61:70,
                                      cLength=cLength[1,1])

mt.d.read.coverage1 <- getReadCoverage(bed[[dSamples['Mutant']]],
                                      readLength=61:70,
                                      cLength=cLength[1,1])

## Calculate read per 1000:
wtDCoverage1 <- mil * wtDCoverage1 / lDepth[[mSamples['WildType']]][1]
mtDCoverage1 <- mil * mtDCoverage1 / lDepth[[mSamples['Mutant']]][1]
wt.d.read.coverage1 <- mil * wt.d.read.coverage1 / lDepth[[mSamples['WildType']]][1]
mt.d.read.coverage1 <- mil * mt.d.read.coverage1 / lDepth[[mSamples['Mutant']]][1]


pdf("figures/2021-02-02/WT-Mutant_disome_readCoverage_61-70.pdf", width=30,
    height=15)
par(mfrow=c(2,1))
plotCoverage(xAcov=wtDCoverage1, xPileup=wt.d.read.coverage1,
             readLength=61:70, main='Wild Type', xlim=c(525,675),
             sequence=covidSeq$WildType, scex=1)
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
plotCoverage(xAcov=mtDCoverage1, xPileup=mt.d.read.coverage1,
             readLength=61:70, main='Mutant', xlim=c(525,675),
             sequence=covidSeq$Mutant, scex=1)
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
dev.off()



######################################################################
######################################################################
### Make a heat map for the read 5'end
zeroOff <- rep(0, length(26:100))
names(zeroOff) <- 26:100

mt5end.d <- getCoverage(bed[[dSamples['Mutant']]], readLength=36:100,
                        offset=zeroOff, cLength=cLength[1,1])
wt5end.d <- getCoverage(bed[[dSamples['WildType']]],
                        readLength=36:100, offset=zeroOff,
                        cLength=cLength[1,1])
## Calculate read per 1000:
mt5end.d <- mil * mt5end.d / lDepth[[mSamples['WildType']]][1]
wt5end.d <- mil * wt5end.d / lDepth[[mSamples['Mutant']]][1]
mt5end.d[mt5end.d<1] <- 0
wt5end.d[wt5end.d<1] <- 0


png('figures/2021-02-02/Mutant_disome_heatmap.png', width=30000,
    height=1000, res=200)
image(x=seq(0.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
      z=log(mt5end.d), xlab='Read 5-end position', ylab='Read length')
abline(v=c(551,581,611,1147))
dev.off()

png('figures/2021-02-02/WildType_disome_heatmap.png', width=30000,
    height=1000, res=200)
image(x=seq(0.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
      z=log(wt5end.d), xlab='Read 5-end position', ylab='Read length')
abline(v=c(551,581,611,1147))
dev.off()

mt5end.d.ext <- extendCoverage(mt5end.d, readLength=36:100)
wt5end.d.ext <- extendCoverage(wt5end.d, readLength=36:100)

### Do teh same for the 3'-end

mt3end.d <- getCoverage(bed[[dSamples['Mutant']]], readLength=36:100,
                        offset=zeroOff, cLength=cLength[1,1], ref=3)
wt3end.d <- getCoverage(bed[[dSamples['WildType']]],
                        readLength=36:100, offset=zeroOff,
                        cLength=cLength[1,1], ref=3)
## Calculate read per 1000:
mt3end.d <- mil * mt3end.d / lDepth[[mSamples['WildType']]][1]
wt3end.d <- mil * wt3end.d / lDepth[[mSamples['Mutant']]][1]
mt3end.d[mt3end.d<1] <- 0
wt3end.d[wt3end.d<1] <- 0

png('figures/2021-02-02/WildType_3end_disome_heatmap.png',
    width=30000, height=1800, res=200)
image(x=seq(0.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
      z=log(wt3end.d), xlab='Read 3-end position',
      ylab='Read length')
abline(v=c(551,581,611,1147))
dev.off()


png('figures/2021-02-02/Mutant_3end_disome_heatmap.png',
    width=30000, height=1800, res=200)
image(x=seq(0.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
      z=log(mt3end.d), xlab='Read 3-end position',
      ylab='Read length')
abline(v=c(551,581,611,1147))
dev.off()


####################
### Zooming in:

png('figures/2021-02-02/WildType_5end_disomes_zoom_heatmap.png',
    width=6000, height=1800, res=200)
image(x=seq(0.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
      z=log(wt5end.d), xlab='5-end position',
      ylab='Read length', xlim=c(500,700), xaxt='n',
      yaxt='n', col=brewer.pal(9, 'Reds'))
axis(1, seq(500, 700, 10))
axis(2, seq(36, 100, 2))
plotSequence(covidSeq$WildType, ylim=c(0,-1250), cex=1)
abline(v=c(551,581,611,1147))
dev.off()


png('figures/2021-02-02/Mutant_5end_disomes_zoom_heatmap.png',
    width=6000, height=1800, res=200)
image(x=seq(0.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
      z=log(mt5end.d), xlab='5-end position',
      ylab='Read length', xlim=c(500,700), xaxt='n',
      yaxt='n', col=brewer.pal(9, 'Reds'))
axis(1, seq(500, 700, 10))
axis(2, seq(36, 100, 2))
plotSequence(covidSeq$Mutant, ylim=c(0,-1250), cex=1)
abline(v=c(551,581,611,1147))
dev.off()

## extended reads
png('figures/2021-02-02/Mutant_coverageByLength_disomes_zoom_heatmap.png',
    width=6000, height=1800, res=200)
image(x=seq(0.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
      z=log(mt5end.d.ext), xlab='Read coverage',
      ylab='Read length', xlim=c(500,700), xaxt='n',
      yaxt='n', col=redCols)
axis(1, seq(500, 700, 10))
axis(2, seq(36, 100, 2))
plotSequence(covidSeq$Mutant, ylim=c(0,-1250), cex=1)
abline(v=c(551,581,611,1147))
dev.off()


####################
### 3'-end
png('figures/2021-02-02/WildType_3end_disomes_zoom_heatmap.png',
    width=6000, height=1800, res=200)
image(x=seq(0.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
      z=log(wt3end.d), xlab='3-end position',
      ylab='Read length', xlim=c(500,700), xaxt='n',
      yaxt='n', col=brewer.pal(9, 'Blues'))
axis(1, seq(500, 700, 10))
axis(2, seq(36, 100, 2))
plotSequence(covidSeq$WildType, ylim=c(0,-1250), cex=1)
abline(v=c(551,581,611,1147))
dev.off()


png('figures/2021-02-02/Mutant_3end_disomes_zoom_heatmap.png',
    width=6000, height=1800, res=200)
image(x=seq(0.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
      z=log(mt3end.d), xlab='3-end position',
      ylab='Read length', xlim=c(500,700), xaxt='n',
      yaxt='n', col=brewer.pal(9, 'Blues'))
axis(1, seq(500, 700, 10))
axis(2, seq(36, 100, 2))
plotSequence(covidSeq$Mutant, ylim=c(0,-1250), cex=1)
abline(v=c(551,581,611,1147))
dev.off()


## Combine 5'-end and 3'-end
png('figures/2021-02-02/WildType_5-3end_disomes_zoom_heatmap.png',
    width=6000, height=2200, res=200)
image(x=seq(0.5, 630+1, 1), y=seq(35.5, 101, 1),
      z=log(wt5end.d[1:630,]), xlab='5 and 3-end position',
      ylab='Read length', xlim=c(550,700), xaxt='n', yaxt='n',
      col=brewer.pal(9, 'Reds'), ylim=c(34.5, 100.5),
      main='Wild type')
par(new=TRUE)
image(x=seq(630.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
      z=log(wt3end.d[631:cLength[1,1],]), xlab='5 and 3-end position',
      ylab='Read length', xlim=c(550,700), xaxt='n',
      yaxt='n', col=brewer.pal(9, 'Blues'), ylim=c(34.5, 100.5))
axis(1, seq(500, 700, 10))
axis(2, seq(36, 100, 2))
plotSequence(covidSeq$WildType, ylim=c(0,-(35*50)), cex=0.8)
abline(v=c(551,581,611,1147))
dev.off()


## png('figures/2021-02-02/WildType_5-3end_zoom2_heatmap.png',
##     width=4000, height=1200, res=200)
## image(x=seq(0.5, 630+1, 1), y=seq(35.5, 101, 1),
##       z=log(wt5end[1:630,]), xlab='5 and 3-end position',
##       ylab='Read length', xlim=c(570,700), xaxt='n',
##       yaxt='n', col=brewer.pal(9, 'Reds'), ylim=c(64.5, 100.5))
## par(new=TRUE)
## image(x=seq(630.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
##       z=log(wt3end[631:cLength[1,1],]), xlab='5 and 3-end position',
##       ylab='Read length', xlim=c(570,700), xaxt='n',
##       yaxt='n', col=brewer.pal(9, 'Blues'), ylim=c(64.5, 100.5))
## axis(1, seq(500, 700, 10))
## axis(2, seq(36, 100, 2))
## plotSequence(covidSeq$WildType, ylim=c(0,-(65*50)), cex=0.8)
## abline(v=c(551,581,611,1147))
## dev.off()


png('figures/2021-02-02/Mutant_5-3end_disomes_zoom_heatmap.png',
    width=6000, height=2200, res=200)
image(x=seq(0.5, 630+1, 1), y=seq(35.5, 101, 1),
      z=log(mt5end.d[1:630,]), xlab='5 and 3-end position',
      ylab='Read length', xlim=c(550,700), xaxt='n', yaxt='n',
      col=brewer.pal(9, 'Reds'), ylim=c(34.5, 100.5), main='Mutant')
par(new=TRUE)
image(x=seq(630.5, cLength[1,1]+1, 1), y=seq(35.5, 101, 1),
      z=log(mt3end.d[631:cLength[1,1],]), xlab='5 and 3-end position',
      ylab='Read length', xlim=c(550,700), xaxt='n',
      yaxt='n', col=brewer.pal(9, 'Blues'), ylim=c(34.5, 100.5))
axis(1, seq(500, 700, 10))
axis(2, seq(36, 100, 2))
plotSequence(covidSeq$Mutant, ylim=c(0,-(35*50)), cex=0.8)
abline(v=c(551,581,611,1147))
dev.off()

######################################################################
## Do exactly the same for the monosome data:
mt5end.m <- getCoverage(bed[[mSamples['Mutant']]], readLength=26:35,
                        offset=zeroOff, cLength=cLength[1,1])
wt5end.m <- getCoverage(bed[[mSamples['WildType']]],
                        readLength=26:35, offset=zeroOff,
                        cLength=cLength[1,1])
mt3end.m <- getCoverage(bed[[mSamples['Mutant']]], readLength=26:35,
                        offset=zeroOff, cLength=cLength[1,1], ref=3)
wt3end.m <- getCoverage(bed[[mSamples['WildType']]],
                        readLength=26:35, offset=zeroOff,
                        cLength=cLength[1,1], ref=3)
mt5end.m <- mil * mt5end.m / lDepth[[mSamples['WildType']]][1]
wt5end.m <- mil * wt5end.m / lDepth[[mSamples['Mutant']]][1]
mt5end.m[mt5end.m<1] <- 0
wt5end.m[wt5end.m<1] <- 0
mt3end.m <- mil * mt3end.m / lDepth[[mSamples['WildType']]][1]
wt3end.m <- mil * wt3end.m / lDepth[[mSamples['Mutant']]][1]
mt3end.m[mt3end.m<1] <- 0
wt3end.m[wt3end.m<1] <- 0

mt5end.m.ext <- extendCoverage(mt5end.m, readLength=26:35)
wt5end.m.ext <- extendCoverage(wt5end.m, readLength=26:35)


png('figures/2021-02-02/WildType_monosome_readCoverageLength_zoom_heatmap.png',
    width=6000, height=800, res=200)
image(x=seq(1, cLength[1,1], 1), y=seq(25.5, 36, 1),
      z=log(wt5end.m.ext), xlab='Read coverage',
      ylab='Read length', xlim=c(550,700), xaxt='n', yaxt='n',
      col=blueCols, ylim=c(24.5, 35.5),
      main='Wild type')
axis(1, seq(500, 700, 10))
axis(2, seq(26, 35, 2))
plotSequence(covidSeq$WildType, ylim=c(0,-(25*50)), cex=0.8)
abline(v=c(551,581,611,1147))
dev.off()

png('figures/2021-02-02/WildType_monosome_5-3end_zoom_heatmap.png',
    width=6000, height=800, res=200)
image(x=seq(0.5, 630+1, 1), y=seq(25.5, 36, 1),
      z=log(wt5end.m[1:630,]), xlab='5 and 3-end position',
      ylab='Read length', xlim=c(550,700), xaxt='n', yaxt='n',
      col=redCols, ylim=c(24.5, 35.5),
      main='Wild type')
## par(new=TRUE)
## image(x=seq(609.5, 630.5, 1), y=seq(25.5, 36, 1),
##       z=log(wt3end.m[610:630,]),
##       xlab='', ylab='', xlim=c(550,700),
##       ylim=c(24.5, 36), xaxt='n', yaxt='n',
##       col=brewer.pal(9, 'Blues'))
par(new=TRUE)
image(x=seq(609.5, cLength[1,1]+1, 1), y=seq(25.5, 36, 1),
      z=log(wt3end.m[610:cLength[1,1],]), xlab='', ylab='',
      xlim=c(550,700), xaxt='n', yaxt='n', col=blueCols,
      ylim=c(24.5, 35.5))
axis(1, seq(500, 700, 10))
axis(2, seq(26, 35, 2))
plotSequence(covidSeq$WildType, ylim=c(0,-(25*50)), cex=0.8)
abline(v=c(551,581,611,1147))
dev.off()

png('figures/2021-02-02/Mutant_monosome_readCoverageLength_zoom_heatmap.png',
    width=6000, height=800, res=200)
image(x=seq(0.5, cLength[1,1], 1), y=seq(25.5, 36, 1),
      z=log(mt5end.m.ext), xlab='Read coverage', ylab='Read length',
      xlim=c(550,700), xaxt='n', yaxt='n',
      col=blueCols,
      ylim=c(24.5, 35.5), main='Mutant')
axis(1, seq(500, 700, 10))
axis(2, seq(26, 35, 2))
plotSequence(covidSeq$WildType, ylim=c(0,-(25*50)), cex=0.8)
abline(v=c(551,581,611,1147))
dev.off()

png('figures/2021-02-02/Mutant_monosome_5-3end_zoom_heatmap.png',
    width=6000, height=800, res=200)
image(x=seq(0.5, 630+1, 1), y=seq(25.5, 36, 1),
      z=log(mt5end.m[1:630,]), xlab='5 and 3-end position',
      ylab='Read length', xlim=c(550,700), xaxt='n', yaxt='n',
      col=redCols, ylim=c(24.5, 35.5),
      main='Mutant')
par(new=TRUE)
image(x=seq(609.5, 630.5, 1), y=seq(25.5, 36, 1),
      z=log(wt3end.m[610:630,]),
      xlab='', ylab='', xlim=c(550,700),
      ylim=c(24.5, 36), xaxt='n', yaxt='n',
      col=blueCols)
par(new=TRUE)
image(x=seq(630.5, cLength[1,1]+1, 1), y=seq(25.5, 36, 1),
      z=log(mt3end.m[631:cLength[1,1],]), xlab='', ylab='',
      xlim=c(550,700), xaxt='n', yaxt='n', col=blueCols,
      ylim=c(24.5, 35.5))
axis(1, seq(500, 700, 10))
axis(2, seq(26, 35, 2))
plotSequence(covidSeq$Mutant, ylim=c(0,-(25*50)), cex=0.8)
abline(v=c(551,581,611,1147))
dev.off()




######################################################################
######################################################################
### Plot read coverage for monosome with best offset
bos <- 13:87
names(bos) <- 26:100
bos['29'] <- 15
## bos['32'] <- 13

## Monosomes
wtMCoverage <- getCoverage(bed[[mSamples['WildType']]],
                           readLength=26:35, offset=bos,
                           cLength=cLength[1,1])
mtMCoverage <- getCoverage(bed[[mSamples['Mutant']]],
                           readLength=26:35, offset=bos,
                           cLength=cLength[1,1])

wt.m.read.coverage <- getReadCoverage(bed[[mSamples['WildType']]],
                                      readLength=26:35,
                                      cLength=cLength[1,1])

mt.m.read.coverage <- getReadCoverage(bed[[mSamples['Mutant']]],
                                      readLength=26:35,
                                      cLength=cLength[1,1])

## Calculate read per 10000:
mil <- 10000
wtMCoverage <- mil * wtMCoverage / lDepth[[mSamples['WildType']]][1]
mtMCoverage <- mil * mtMCoverage / lDepth[[mSamples['Mutant']]][1]
wt.m.read.coverage <- mil * wt.m.read.coverage / lDepth[[mSamples['WildType']]][1]
mt.m.read.coverage <- mil * mt.m.read.coverage / lDepth[[mSamples['Mutant']]][1]

pdf("figures/2021-02-02/WT-Mutant_Monosomes_readCoverage_bestOffset.pdf", width=30,
    height=15)
par(mfrow=c(2,1))
plotCoverage(xAcov=wtMCoverage, xPileup=wt.m.read.coverage,
             readLength=26:35, main='Wild Type', xlim=c(525,675),
             sequence=covidSeq$WildType, scex=1)
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
plotCoverage(xAcov=mtMCoverage, xPileup=mt.m.read.coverage,
             readLength=26:35, main='Mutant', xlim=c(525,675),
             sequence=covidSeq$Mutant, scex=1)
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
dev.off()





######################################################################
######################################################################
### Divide reads that map the the target region from reads that map
### outside of it and plot them together as pileup. Concentrate on
### reds with length 28 and 29 only
targetRegion <- c(590, 630)
outsideRegion <- c(500, 700)

wtMono <- bed[[mSamples['WildType']]]
mtMono <- bed[[mSamples['Mutant']]]

mtMonoT <- mtMono[which(mtMono[,2] >= targetRegion[1] & mtMono[,3] <= targetRegion[2]),]
wtMonoT <- wtMono[which(wtMono[,2] >= targetRegion[1] & wtMono[,3] <= targetRegion[2]),]

mtMonoO <- rbind(mtMono[which(mtMono[,2] > outsideRegion[1] & mtMono[,2] <= targetRegion[1]),],
                 mtMono[which(mtMono[,2] > targetRegion[2] & mtMono[,2] <= outsideRegion[2]),])
wtMonoO <- rbind(wtMono[which(wtMono[,2] > outsideRegion[1] & wtMono[,2] <= targetRegion[1]),],
                 wtMono[which(wtMono[,2] > targetRegion[2] & wtMono[,2] <= outsideRegion[2]),])

## Monosomes in target region
wtMCoverageT <- getCoverage(wtMonoT,
                           readLength=26:35, offset=bos,
                           cLength=cLength[1,1])
mtMCoverageT <- getCoverage(mtMonoT,
                           readLength=26:35, offset=bos,
                           cLength=cLength[1,1])
wt.m.read.coverageT <- getReadCoverage(wtMonoT,
                                      readLength=26:35,
                                      cLength=cLength[1,1])
mt.m.read.coverageT <- getReadCoverage(mtMonoT,
                                      readLength=26:35,
                                      cLength=cLength[1,1])

## Outside region:
wtMCoverageO <- getCoverage(wtMonoO,
                           readLength=26:35, offset=bos,
                           cLength=cLength[1,1])
mtMCoverageO <- getCoverage(mtMonoO,
                           readLength=26:35, offset=bos,
                           cLength=cLength[1,1])
wt.m.read.coverageO <- getReadCoverage(wtMonoO,
                                      readLength=26:35,
                                      cLength=cLength[1,1])
mt.m.read.coverageO <- getReadCoverage(mtMonoO,
                                      readLength=26:35,
                                      cLength=cLength[1,1])

## Calculate read per 10000:
mil <- 10000
wtMCoverageT <- mil * wtMCoverageT / lDepth[[mSamples['WildType']]][1]
mtMCoverageT <- mil * mtMCoverageT / lDepth[[mSamples['Mutant']]][1]
wt.m.read.coverageT <- mil * wt.m.read.coverageT / lDepth[[mSamples['WildType']]][1]
mt.m.read.coverageT <- mil * mt.m.read.coverageT / lDepth[[mSamples['Mutant']]][1]
wtMCoverageO <- mil * wtMCoverageO / lDepth[[mSamples['WildType']]][1]
mtMCoverageO <- mil * mtMCoverageO / lDepth[[mSamples['Mutant']]][1]
wt.m.read.coverageO <- mil * wt.m.read.coverageO / lDepth[[mSamples['WildType']]][1]
mt.m.read.coverageO <- mil * mt.m.read.coverageO / lDepth[[mSamples['Mutant']]][1]



pdf("figures/2021-02-02/WT-Mutant_Monosomes_readCoverage_targetRegion.pdf",
    width=30, height=15)
par(mfrow=c(2,1))
plotCoverage(xAcov=wtMCoverageT, xPileup=wt.m.read.coverageT,
             readLength=26:35, main='Wild Type', xlim=c(525,675),
             sequence=covidSeq$WildType, scex=1, ylim=c(0,400))
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
plotCoverage(xAcov=mtMCoverageT, xPileup=mt.m.read.coverageT,
             readLength=26:35, main='Mutant', xlim=c(525,675),
             sequence=covidSeq$Mutant, scex=1, ylim=c(0,800))
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
dev.off()

pdf("figures/2021-02-02/WT-Mutant_Monosomes_readCoverage_outside.pdf",
    width=30, height=15)
par(mfrow=c(2,1))
plotCoverage(xAcov=wtMCoverageO, xPileup=wt.m.read.coverageO,
             readLength=26:35, main='Wild Type', xlim=c(525,675),
             sequence=covidSeq$WildType, scex=1, ylim=c(0,400))
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
plotCoverage(xAcov=mtMCoverageO, xPileup=mt.m.read.coverageO,
             readLength=26:35, main='Mutant', xlim=c(525,675),
             sequence=covidSeq$Mutant, scex=1, ylim=c(0,800))
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
dev.off()


######################################################################
### Do the same for the disomes (map the A-site of the leading ribosome)

disomeTregion <- c(550, 650)
disomeOregion <- c(400, 700)

wtDisome <- bed[[dSamples['WildType']]]
mtDisome <- bed[[dSamples['Mutant']]]

mtDisomeT <- mtDisome[which(mtDisome[,2] >= disomeTregion[1] & mtDisome[,3] <= disomeTregion[2]),]
wtDisomeT <- wtDisome[which(wtDisome[,2] >= disomeTregion[1] & wtDisome[,3] <= disomeTregion[2]),]

mtDisomeO <- rbind(mtDisome[which(mtDisome[,2] > disomeOregion[1] & mtDisome[,2] <= disomeTregion[1]),],
                   mtDisome[which(mtDisome[,2] > disomeTregion[2] & mtDisome[,2] <= disomeOregion[2]),])
wtDisomeO <- rbind(wtDisome[which(wtDisome[,2] > disomeOregion[1] & wtDisome[,2] <= disomeTregion[1]),],
                 wtDisome[which(wtDisome[,2] > disomeTregion[2] & wtDisome[,2] <= disomeOregion[2]),])

## Disomesomes in target region
wtDCoverageT <- getCoverage(wtDisomeT,
                           readLength=50:70, offset=bestOffset+30,
                           cLength=cLength[1,1])
mtDCoverageT <- getCoverage(mtDisomeT,
                           readLength=50:70, offset=bestOffset+30,
                           cLength=cLength[1,1])
wt.d.read.coverageT <- getReadCoverage(wtDisomeT,
                                      readLength=50:70,
                                      cLength=cLength[1,1])
mt.d.read.coverageT <- getReadCoverage(mtDisomeT,
                                      readLength=50:70,
                                      cLength=cLength[1,1])

## Outside region:
wtDCoverageO <- getCoverage(wtDisomeO,
                           readLength=50:70, offset=bestOffset+30,
                           cLength=cLength[1,1])
mtDCoverageO <- getCoverage(mtDisomeO,
                           readLength=50:70, offset=bestOffset+30,
                           cLength=cLength[1,1])
wt.d.read.coverageO <- getReadCoverage(wtDisomeO,
                                      readLength=50:70,
                                      cLength=cLength[1,1])
mt.d.read.coverageO <- getReadCoverage(mtDisomeO,
                                      readLength=50:70,
                                      cLength=cLength[1,1])

## Calculate read per 10000:
mil <- 10000
wtDCoverageT <- mil * wtDCoverageT / lDepth[[dSamples['WildType']]][1]
mtDCoverageT <- mil * mtDCoverageT / lDepth[[dSamples['Mutant']]][1]
wt.d.read.coverageT <- mil * wt.d.read.coverageT / lDepth[[dSamples['WildType']]][1]
mt.d.read.coverageT <- mil * mt.d.read.coverageT / lDepth[[dSamples['Mutant']]][1]
wtDCoverageO <- mil * wtDCoverageO / lDepth[[dSamples['WildType']]][1]
mtDCoverageO <- mil * mtDCoverageO / lDepth[[dSamples['Mutant']]][1]
wt.d.read.coverageO <- mil * wt.d.read.coverageO / lDepth[[dSamples['WildType']]][1]
mt.d.read.coverageO <- mil * mt.d.read.coverageO / lDepth[[dSamples['Mutant']]][1]



pdf("figures/2021-02-02/WT-Mutant_disome_readCoverage_targetRegion.pdf", width=30,
    height=15)
par(mfrow=c(2,1))
plotCoverage(xAcov=wtDCoverageT, xPileup=wt.d.read.coverageT,
             readLength=50:70, main='Wild Type', xlim=c(450,675),
             sequence=covidSeq$WildType, scex=0.8, ylim=c(0,270))
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
plotCoverage(xAcov=mtDCoverageT, xPileup=mt.d.read.coverageT,
             readLength=50:70, main='Mutant', xlim=c(450,675),
             sequence=covidSeq$Mutant, scex=0.8, ylim=c(0,270))
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
dev.off()

pdf("figures/2021-02-02/WT-Mutant_disome_readCoverage_outsideRegion.pdf", width=30,
    height=15)
par(mfrow=c(2,1))
plotCoverage(xAcov=wtDCoverageO, xPileup=wt.d.read.coverageO,
             readLength=50:70, main='Wild Type', xlim=c(450,675),
             sequence=covidSeq$WildType, scex=0.8, ylim=c(0,270))
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
plotCoverage(xAcov=mtDCoverageO, xPileup=mt.d.read.coverageO,
             readLength=50:70, main='Mutant', xlim=c(450,675),
             sequence=covidSeq$Mutant, scex=0.8, ylim=c(0,270))
segments(x0=c(551, 581, 611), y0=rep(0,4), x1=c(551, 581, 611),
         y1=rep(800,4), col='red', lty=2)
dev.off()


######################################################################
######################################################################
### Study best offset for different reads length (monosome only)
######################################################################
######################################################################
offset14 <- rep(14, 10)
offset15 <- offset14 + 1
offset13 <- offset14 - 1
names(offset14) <- names(offset15) <- names(offset13) <- 26:35


wtFrame13 <- getFrame(bed[[mSamples['WildType']]], readLength=26:35,
                      cds=cds, offset=offset13, region=c(110,500))
wtFrame13sum <- tapply(wtFrame13[,1], wtFrame13[,2], sum)
wtReads13Sum <- tapply(wtFrame13[,1], wtFrame13[,3], sum)

wtFrame15 <- getFrame(bed[[mSamples['WildType']]], readLength=26:35,
                      cds=cds, offset=offset15, region=c(110,500))
wtFrame15sum <- tapply(wtFrame15[,1], wtFrame15[,2], sum)
wtReads15Sum <- tapply(wtFrame15[,1], wtFrame15[,3], sum)

wtFrame14 <- getFrame(bed[[mSamples['WildType']]], readLength=26:35,
                      cds=cds, offset=offset14, region=c(110,500))
wtFrame14sum <- tapply(wtFrame14[,1], wtFrame14[,2], sum)
wtReads14Sum <- tapply(wtFrame14[,1], wtFrame14[,3], sum)



wtFrame15all <- data.frame(Frame = factor(wtFrame15[,2]),
                           Counts = wtFrame15[,1]/wtReads15Sum[as.character(wtFrame15[,3])],
                           Length = factor(wtFrame15[,3]))
fp2 <- ggplot(wtFrame15all, aes(x=Length, y=Counts, fill=Frame)) +
    geom_bar(stat='identity', position=position_dodge()) +
    ggtitle('WT sample - 15nt offset') +
    theme_minimal()
ggsave("figures/2021-02-02/WildType_frames_15.pdf", fp2, width=6)

wtFrame14all <- data.frame(Frame = factor(wtFrame14[,2]),
                           Counts = wtFrame14[,1]/wtReads14Sum[as.character(wtFrame14[,3])],
                           Length = factor(wtFrame14[,3]))
fp2 <- ggplot(wtFrame14all, aes(x=Length, y=Counts, fill=Frame)) +
    geom_bar(stat='identity', position=position_dodge()) +
    ggtitle('WT sample - 14nt offset') +
    theme_minimal()
ggsave("figures/2021-02-02/WildType_frames_14.pdf", fp2, width=6)

wtFrame13all <- data.frame(Frame = factor(wtFrame13[,2]),
                           Counts = wtFrame13[,1]/wtReads13Sum[as.character(wtFrame13[,3])],
                           Length = factor(wtFrame13[,3]))
fp2 <- ggplot(wtFrame13all, aes(x=Length, y=Counts, fill=Frame)) +
    geom_bar(stat='identity', position=position_dodge()) +
    ggtitle('WT sample - 13nt offset') +
    theme_minimal()
ggsave("figures/2021-02-02/WildType_frames_13.pdf", fp2, width=6)


## Now the mutant

mtFrame13 <- getFrame(bed[[mSamples['Mutant']]], readLength=26:35,
                      cds=cds, offset=offset13, region=c(110,500))
mtFrame13sum <- tapply(mtFrame13[,1], mtFrame13[,2], sum)
mtReads13Sum <- tapply(mtFrame13[,1], mtFrame13[,3], sum)
mtFrame15 <- getFrame(bed[[mSamples['Mutant']]], readLength=26:35,
                      cds=cds, offset=offset15, region=c(110,500))
mtFrame15sum <- tapply(mtFrame15[,1], mtFrame15[,2], sum)
mtReads15Sum <- tapply(mtFrame15[,1], mtFrame15[,3], sum)
mtFrame14 <- getFrame(bed[[mSamples['Mutant']]], readLength=26:35,
                      cds=cds, offset=offset14, region=c(110,500))
mtFrame14sum <- tapply(mtFrame14[,1], mtFrame14[,2], sum)
mtReads14Sum <- tapply(mtFrame14[,1], mtFrame14[,3], sum)



mtFrame15all <- data.frame(Frame = factor(mtFrame15[,2]),
                           Counts = mtFrame15[,1]/mtReads15Sum[as.character(mtFrame15[,3])],
                           Length = factor(mtFrame15[,3]))
fp2 <- ggplot(mtFrame15all, aes(x=Length, y=Counts, fill=Frame)) +
    geom_bar(stat='identity', position=position_dodge()) +
    ggtitle('MT sample - 15nt offset') +
    theme_minimal()
ggsave("figures/2021-02-02/Mutant_frames_15.pdf", fp2, width=6)

mtFrame14all <- data.frame(Frame = factor(mtFrame14[,2]),
                           Counts = mtFrame14[,1]/mtReads14Sum[as.character(mtFrame14[,3])],
                           Length = factor(mtFrame14[,3]))
fp2 <- ggplot(mtFrame14all, aes(x=Length, y=Counts, fill=Frame)) +
    geom_bar(stat='identity', position=position_dodge()) +
    ggtitle('MT sample - 14nt offset') +
    theme_minimal()
ggsave("figures/2021-02-02/Mutant_frames_14.pdf", fp2, width=6)

mtFrame13all <- data.frame(Frame = factor(mtFrame13[,2]),
                           Counts = mtFrame13[,1]/mtReads13Sum[as.character(mtFrame13[,3])],
                           Length = factor(mtFrame13[,3]))
fp2 <- ggplot(mtFrame13all, aes(x=Length, y=Counts, fill=Frame)) +
    geom_bar(stat='identity', position=position_dodge()) +
    ggtitle('MT sample - 13nt offset') +
    theme_minimal()
ggsave("figures/2021-02-02/Mutant_frames_13.pdf", fp2, width=6)



######################################################################
### Plot the frame preference before and after the site in all samples
### (monosomes)

## Wild Type
wtFrame <- getFrame(bed[[mSamples['WildType']]], readLength=28:29,
                    cds=cds, offset=bos, region=c(110,550))
wtFramesum <- tapply(wtFrame[,1], wtFrame[,2], sum)
wtReadsSum <- tapply(wtFrame[,1], wtFrame[,3], sum)
wtFrameAfter <- getFrame(bed[[mSamples['WildType']]],
                         readLength=28:29, cds=cds, offset=bos,
                         region=c(630, cLength[1,1]))
wtFrameAftersum <- tapply(wtFrameAfter[,1], wtFrameAfter[,2], sum)
wtReadsSum <- tapply(wtFrameAfter[,1], wtFrameAfter[,3], sum)

## Mutant
mtFrame <- getFrame(bed[[mSamples['Mutant']]], readLength=28:29,
                    cds=cds, offset=bos, region=c(110,550))
mtFramesum <- tapply(mtFrame[,1], mtFrame[,2], sum)
mtReadsSum <- tapply(mtFrame[,1], mtFrame[,3], sum)
mtFrameAfter <- getFrame(bed[[mSamples['Mutant']]], readLength=28:29,
                         cds=cds, offset=bos,
                         region=c(630, cLength[1,1]))
mtFrameAftersum <- tapply(mtFrameAfter[,1], mtFrameAfter[,2], sum)
mtReadsSum <- tapply(mtFrameAfter[,1], mtFrameAfter[,3], sum)

## Wild Type Cyclohex.
wtCxFrame <- getFrame(bed[[mSamplesCX['WildType']]], readLength=28:29,
                    cds=cds, offset=bos, region=c(110,550))
wtCxFramesum <- tapply(wtCxFrame[,1], wtCxFrame[,2], sum)
wtCxReadsSum <- tapply(wtCxFrame[,1], wtCxFrame[,3], sum)
wtCxFrameAfter <- getFrame(bed[[mSamplesCX['WildType']]],
                         readLength=28:29, cds=cds, offset=bos,
                         region=c(630, cLength[1,1]))
wtCxFrameAftersum <- tapply(wtCxFrameAfter[,1], wtCxFrameAfter[,2], sum)
wtCxReadsSum <- tapply(wtCxFrameAfter[,1], wtCxFrameAfter[,3], sum)

## Mutant Cyclohex.
mtCxFrame <- getFrame(bed[[mSamplesCX['Mutant']]], readLength=28:29,
                    cds=cds, offset=bos, region=c(110,550))
mtCxFramesum <- tapply(mtCxFrame[,1], mtCxFrame[,2], sum)
mtCxReadsSum <- tapply(mtCxFrame[,1], mtCxFrame[,3], sum)
mtCxFrameAfter <- getFrame(bed[[mSamplesCX['Mutant']]], readLength=28:29,
                         cds=cds, offset=bos,
                         region=c(630, cLength[1,1]))
mtCxFrameAftersum <- tapply(mtCxFrameAfter[,1], mtCxFrameAfter[,2], sum)
mtCxReadsSum <- tapply(mtCxFrameAfter[,1], mtCxFrameAfter[,3], sum)



frameSumData <- data.frame(Frame  = factor(rep(names(mtFramesum), 8)),
                           Counts = c(mtFramesum/sum(mtFramesum),
                                      mtFrameAftersum/sum(mtFrameAftersum),
                                      wtFramesum/sum(wtFramesum),
                                      wtFrameAftersum/sum(wtFrameAftersum),
                                      mtCxFramesum/sum(mtCxFramesum),
                                      mtCxFrameAftersum/sum(mtCxFrameAftersum),
                                      wtCxFramesum/sum(wtCxFramesum),
                                      wtCxFrameAftersum/sum(wtCxFrameAftersum)
                                      ),
                           Position = factor(rep(rep(c('Before','After'), each=3), 4),
                                             levels=c('Before','After')),
                           Geno = rep(c(rep('Mutant', 6),
                                        rep('WildType', 6)),
                                      2),
                           Treatment = c(rep('None', 12),
                                           rep('Cycloheximide', 12))
                           )

pf3 <- ggplot(frameSumData, aes(x=Position, y=Counts, fill=Frame)) +
    geom_bar(stat='identity', position=position_dodge()) +
    facet_grid(rows=vars(Geno), cols=vars(Treatment)) +
    ggtitle('Frame preference')
ggsave("figures/2021-02-02/frame_preferences_all_samples.pdf", pf3)



######################################################################
######################################################################
### Calculate frame preference along the construct:

dis <- 100
codonFrame <- matrix(ncol=3, nrow=(1002-158)/3)
rownames(codonFrame) <- seq(158, (1002-2), 3)
colnames(codonFrame) <- c(0,1,-1)
head(codonFrame)
for(I in rownames(codonFrame)){
    K <- as.numeric(I)
    r <- c(K-dis,K+dis)
    f <- getFrame(bed[[mSamples['WildType']]], readLength=28:29, cds=cds,
                  offset=offset15, region=r)
    if(sum(f[1:3,1]) > 300)
        codonFrame[I,] <- f[1:3,1] + f[4:6,1]
    else
        codonFrame[I,] <- c(0,0,0)
}

codonFrame[1:110,]
frameMax <- unlist(apply(codonFrame/rowSums(codonFrame), 1, which.max))

codonFrameFreq <- codonFrame/rowSums(codonFrame)

### Now the mutant
codonFrameMt <- matrix(ncol=3, nrow=(1002-158)/3)
rownames(codonFrameMt) <- seq(158, (1002-2), 3)
colnames(codonFrameMt) <- c(0,1,-1)
head(codonFrameMt)
for(I in rownames(codonFrameMt)){
    K <- as.numeric(I)
    r <- c(K-dis,K+dis)
    f <- getFrame(bed[[mSamples['Mutant']]], readLength=28:29, cds=cds,
                  offset=offset15, region=r)
    if(sum(f[1:3,1]) > 300)
        codonFrameMt[I,] <- f[1:3,1] + f[4:6,1]
    else
        codonFrameMt[I,] <- c(0,0,0)
}

codonFrameMt[100:110,]
frameMax <- unlist(apply(codonFrameMt/rowSums(codonFrameMt), 1, which.max))

codonFrameMtFreq <- codonFrameMt/rowSums(codonFrameMt)

####################
### Put them together

frameData1 <- data.frame(Position   = rep(as.numeric(rownames(codonFrameFreq)), 3),
                         Frequency  = c(codonFrameFreq[,1],
                                        codonFrameFreq[,2],
                                        codonFrameFreq[,3]),
                         Frame      = factor(c(rep('0', dim(codonFrameFreq)[1]),
                                               rep('1', dim(codonFrameFreq)[1]),
                                               rep('-1', dim(codonFrameFreq)[1]))),
                         Genotype   = factor(rep('Wild Type', 3*dim(codonFrameFreq)[1]),
                                             levels=c('Wild Type', 'Mutant'))
                        )
frameData2 <- data.frame(Position   = rep(as.numeric(rownames(codonFrameMtFreq)), 3),
                         Frequency  = c(codonFrameMtFreq[,1],
                                        codonFrameMtFreq[,2],
                                        codonFrameMtFreq[,3]),
                         Frame      = factor(c(rep('0', dim(codonFrameMtFreq)[1]),
                                               rep('1', dim(codonFrameMtFreq)[1]),
                                               rep('-1', dim(codonFrameMtFreq)[1]))),
                         Genotype   = factor(rep('Mutant', 3*dim(codonFrameMtFreq)[1]),
                                             levels=c('Wild Type', 'Mutant'))
                         )
frameData <- rbind(frameData1, frameData2)


fp <- ggplot(frameData, aes(x=Position, y=Frequency, col=Frame,
                            linetype=Genotype)) +
    geom_vline(xintercept=610, alpha=0.5) +
    geom_smooth(span=0.5, size=1) +
    ## geom_smooth(method = lm, formula = y ~ splines::bs(x, 3)) +
    ggtitle('Frame preference along the construct') +
    theme_classic()
ggsave('figures/2021-02-02/constructFrame.pdf', height=4, width=10)


######################################################################
######################################################################
### Calculate frame preference along the construct for CX samples:

dis <- 100
codonFrameCx <- matrix(ncol=3, nrow=(1002-158)/3)
rownames(codonFrameCx) <- seq(158, (1002-2), 3)
colnames(codonFrameCx) <- c(0,1,-1)
head(codonFrameCx)
for(I in rownames(codonFrameCx)){
    K <- as.numeric(I)
    r <- c(K-dis,K+dis)
    f <- getFrame(bed[[mSamplesCX['WildType']]], readLength=28:29, cds=cds,
                  offset=offset15, region=r)
    if(sum(f[1:3,1]) > 300)
        codonFrameCx[I,] <- f[1:3,1] + f[4:6,1]
    else
        codonFrameCx[I,] <- c(0,0,0)
}

codonFrameCx[1:110,]
frameMax <- unlist(apply(codonFrameCx/rowSums(codonFrameCx), 1, which.max))

codonFrameCxFreq <- codonFrameCx/rowSums(codonFrameCx)

### Now the mutant
codonFrameCxMt <- matrix(ncol=3, nrow=(1002-158)/3)
rownames(codonFrameCxMt) <- seq(158, (1002-2), 3)
colnames(codonFrameCxMt) <- c(0,1,-1)
head(codonFrameCxMt)
for(I in rownames(codonFrameCxMt)){
    K <- as.numeric(I)
    r <- c(K-dis,K+dis)
    f <- getFrame(bed[[mSamplesCX['Mutant']]], readLength=28:29, cds=cds,
                  offset=offset15, region=r)
    if(sum(f[1:3,1]) > 300)
        codonFrameCxMt[I,] <- f[1:3,1] + f[4:6,1]
    else
        codonFrameCxMt[I,] <- c(0,0,0)
}

codonFrameCxMt[100:110,]
frameMax <- unlist(apply(codonFrameCxMt/rowSums(codonFrameCxMt), 1, which.max))

codonFrameCxMtFreq <- codonFrameCxMt/rowSums(codonFrameCxMt)

frameDataCx1 <- data.frame(Position   = rep(as.numeric(rownames(codonFrameCxFreq)), 3),
                         Frequency  = c(codonFrameCxFreq[,1],
                                        codonFrameCxFreq[,2],
                                        codonFrameCxFreq[,3]),
                         Frame      = factor(c(rep('0', dim(codonFrameCxFreq)[1]),
                                               rep('1', dim(codonFrameCxFreq)[1]),
                                               rep('-1', dim(codonFrameCxFreq)[1]))),
                         Genotype   = factor(rep('Wild Type', 3*dim(codonFrameCxFreq)[1]),
                                             levels=c('Wild Type', 'Mutant'))
                        )
frameDataCx2 <- data.frame(Position   = rep(as.numeric(rownames(codonFrameCxMtFreq)), 3),
                         Frequency  = c(codonFrameCxMtFreq[,1],
                                        codonFrameCxMtFreq[,2],
                                        codonFrameCxMtFreq[,3]),
                         Frame      = factor(c(rep('0', dim(codonFrameCxMtFreq)[1]),
                                               rep('1', dim(codonFrameCxMtFreq)[1]),
                                               rep('-1', dim(codonFrameCxMtFreq)[1]))),
                         Genotype   = factor(rep('Mutant', 3*dim(codonFrameCxMtFreq)[1]),
                                             levels=c('Wild Type', 'Mutant'))
                         )
frameDataCx <- rbind(frameDataCx1, frameDataCx2)


fp <- ggplot(frameDataCx, aes(x=Position, y=Frequency, col=Frame,
                            linetype=Genotype)) +
    geom_vline(xintercept=610, alpha=0.5) +
    geom_smooth(span=0.5, size=1) +
    ## geom_smooth(method = lm, formula = y ~ splines::bs(x, 3)) +
    ggtitle('Frame preference along the construct - Cycloheximide') +
    theme_classic()
ggsave('figures/2021-02-02/constructFrameCyclohex.pdf', height=4, width=10)


######################################################################
######################################################################
### Do the read cumulative sum along the transcrip
wt.cSum <- readCumSum(bed=rbind(bed[[dSamples['WildType']]],
                                bed[[mSamples['WildType']]]),
                      cLength=cLength[1,1])
mt.cSum <- readCumSum(bed=rbind(bed[[dSamples['Mutant']]],
                                bed[[mSamples['Mutant']]]),
                      cLength=cLength[1,1])

### Normalisation factor before teh frame shifting site
nFactor <- wt.cSum[500]/mt.cSum[500]

data.cSum <- data.frame(Position = rep(1:cLength[1,1], 2),
                        cumSum   = c(wt.cSum/nFactor,
                                     mt.cSum),
                        Genotype = factor(rep(c('Wild Type', 'Mutant'),
                                              each=cLength[1,1]),
                                          levels=c('Wild Type', 'Mutant'))
                        )
head(data.cSum)

pcsum <- ggplot(data.cSum, aes(x=Position, y=cumSum, col=Genotype)) +
    geom_vline(xintercept=610, alpha=0.5) +
    geom_line(size=2) +
    xlab('Position') +
    ylab('Reads cumulative sum') +
    theme_light()
ggsave('figures/2021-02-02/constructReadCumSum.pdf', pcsum, width=7,
       height=6)


## Do the CHX samples
wt.cSum <- readCumSum(bed=rbind(bed[[dSamplesCX['WildType']]],
                                bed[[mSamplesCX['WildType']]]),
                      cLength=cLength[1,1])
mt.cSum <- readCumSum(bed=rbind(bed[[dSamplesCX['Mutant']]],
                                bed[[mSamplesCX['Mutant']]]),
                      cLength=cLength[1,1])

### Normalisation factor before teh frame shifting site
nFactor <- wt.cSum[500]/mt.cSum[500]

data.cSum <- data.frame(Position = rep(1:cLength[1,1], 2),
                        cumSum   = c(wt.cSum/nFactor,
                                     mt.cSum),
                        Genotype = factor(rep(c('Wild Type', 'Mutant'),
                                              each=cLength[1,1]),
                                          levels=c('Wild Type', 'Mutant'))
                        )
head(data.cSum)

pcsum <- ggplot(data.cSum, aes(x=Position, y=cumSum, col=Genotype)) +
    geom_vline(xintercept=610, alpha=0.5) +
    geom_line(size=2) +
    xlab('Position') +
    ylab('Reads cumulative sum') +
    theme_light()
ggsave('figures/2021-02-02/constructReadCumSum_CHX.pdf', pcsum, width=7,
       height=6)
