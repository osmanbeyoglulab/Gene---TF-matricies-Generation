rm(list=ls())
options(stringsAsFactors = FALSE)
library(ChIPpeakAnno)
library(stringr)
library(GenomicAlignments)
library(data.table)
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(org.Hs.eg.db)
#######
require(BSgenome.Hsapiens.NCBI.GRCh38)
require(annotables)
# Load peak_motifs
loadFile <- "/Volumes/SeagateBackupPlusDrive/Drake/Summer2020/ATAC-Seq/fimo_post-process/peak-consensus-p-0.00001-PAN-250.csv"
saveFile <- "/Volumes/SeagateBackupPlusDrive/Drake/Summer2020/ATAC-Seq/fimo_post-process/gene-consensus-p-0.00001-PAN-250.csv"

if(i != 1)
{
  loadFile <- sub(cancer.types[i-1], cancer.types[i], loadFile)
  saveFile <- sub(cancer.types[i-1], cancer.types[i], saveFile)
}
  
peak_motif <- read.csv(loadFile, row.names = 1, check.names = T)


grch38 <- as.data.frame(grch38)

###
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(BSgenome.Hsapiens.UCSC.hg38)
# library(phastCons100way.UCSC.hg38)
# phast_hg38 <- phastCons100way.UCSC.hg38
###
# require(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(phastCons100way.UCSC.hg19)
# phast_hg19 <- phastCons100way.UCSC.hg19
###########

xget <- function(x, envir, mode, ifnotfound=NA, inherits, 
                 output=c("all", "first")){
  output <- match.arg(output)
  y <- mget(x=x, envir=envir, mode=mode, 
            ifnotfound=ifnotfound, inherits=inherits)
  switch(output, 
         all=sapply(y, paste, collapse=";"),
         first=sapply(y, `[`, 1)
  )
}

##########################FILTER TFs
annoData <- annoGR(TxDb.Hsapiens.UCSC.hg38.knownGene)
#keepTF <- setdiff(intersect(colnames(peak_motif),rownames(bulkrna)),union(housekeepingTF,hkgenes[,1]))
#length(keepTF)
#peak_motif <- peak_motif[,keepTF]
#peak_motif <- peak_motif[,-which(colSums(peak_motif !=0) > nrow(peak_motif)*0.1 | colSums(peak_motif !=0) < nrow(peak_motif)*0.001)]
####
#meanTFexp <- rowMeans(bulkrna[intersect(colnames(peak_motif),rownames(bulkrna)),])
#sdTFexp <- rowSds(data.matrix(bulkrna[intersect(colnames(peak_motif),rownames(bulkrna)),])); names(sdTFexp) <- names(meanTFexp)
#quantile <- quantile(meanTFexp, probs = c(0.2))
#removeTF <- names(which(meanTFexp < quantile[1]))
#peak_motif <- peak_motif[,setdiff(colnames(peak_motif),union(removeTF,names(which(sdTFexp < quantile(sdTFexp, probs = 0.25)))))]
#dim(peak_motif)
############

require(vegan)
library(corrplot)
library(caret)
a<- t(peak_motif)
vare.dist <- vegdist(data.frame(a), "jaccard", binary=TRUE)
b <- as.matrix(vare.dist)
b <- 1 - b
highlyCor <- findCorrelation(b , 0.7)
sort(colnames(peak_motif)[highlyCor])
peak_motif <- peak_motif[,setdiff(colnames(peak_motif),colnames(peak_motif)[highlyCor])]
dim(peak_motif)
#########################################################################################################
#########################################################################################################
peaks <- rownames(peak_motif)
chr <- unlist(lapply(peaks, function(x) strsplit(x, "_")[[1]][1]))
coordinates <-  unlist(lapply(peaks, function(x) strsplit(x, "_")[[1]][2]))
start <- unlist(lapply(coordinates, function(x) strsplit(x, "-")[[1]][1]))
end <- unlist(lapply(coordinates, function(x) strsplit(x, "-")[[1]][2]))

gr <-     GRanges(seqnames = chr,
                  ranges = IRanges(start = strtoi(start),
                                   end = strtoi(end),
                                   names = peaks))
#gr_phast <- gscores(phast_hg19, GRanges(gr))
############################
annotated = annotatePeakInBatch (gr, AnnotationData = annoData, #featureType = c("TSS"), 
                                 PeakLocForDistance="middle", multiple=FALSE, 
                                 output="shortestDistance", 
                                 FeatureLocForDistance="TSS", select="first")
annotated$symbol <- xget(annotated$feature, org.Hs.egSYMBOL)
annotated <- as.data.frame(annotated)
annotated <- annotated[!is.na(annotated$symbol),]
annotated$peaks <-  paste0(annotated$seqnames, "_",annotated$start, "-", annotated$end)
rownames(annotated) <- annotated$peaks
peak_motif$symbol <- annotated[rownames(peak_motif),"symbol"]
dim(peak_motif)
dist2TSS <- 2000
filter_peaks <- rownames(annotated)[which((annotated$distancetoFeature > -dist2TSS & annotated$distancetoFeature < dist2TSS)
                                          | annotated$insideFeature == "inside"
                                          | annotated$insideFeature == "overlapStart"
                                          | annotated$insideFeature == "overlapEnd"
                                          | annotated$insideFeature == "includeFeature"
)]
length(filter_peaks)
peak_motif2 <- peak_motif[filter_peaks,]
#
a<- t(peak_motif2[,-ncol(peak_motif2)])
vare.dist <- vegdist(data.frame(a), "jaccard", binary=TRUE)
b <- as.matrix(vare.dist)
b <- 1 - b
highlyCor <- findCorrelation(b , 0.7)
sort(colnames(peak_motif2)[highlyCor])
peak_motif2 <- peak_motif2[,setdiff(colnames(peak_motif2),colnames(peak_motif2)[highlyCor])]
dim(peak_motif2)
#####
library(data.table) 
tfDT <- data.table(peak_motif2)
allTFs <- paste0(setdiff(colnames(tfDT), c("symbol")), "=sum(", setdiff(colnames(tfDT), "symbol"), ")")
argsToPaste <- as.list(allTFs)
argsToPaste[["sep"]] <- ","
command <- paste0("data_TF_merge <- tfDT[, list(", do.call(paste, argsToPaste), "), by = symbol]")
eval(parse(text = command))             
D = data.frame(data_TF_merge)
rownames(D) = D[,1]
D = D[,-1]
#D = replace(D, D>1, 1)
sort(apply(D, 2, function(x) (length(x) - length(which(x==0)))/length(x)))
write.csv(D, file=saveFile, quote = F)

#############
D1 <- read.csv("/Volumes/SeagateBackupPlusDrive/Drake/Summer2020/ATAC-Seq/fimo_post-process/peak-consensus-p-0.00001-CHOL-250-gene.csv")
D2 <- read.csv("/Volumes/SeagateBackupPlusDrive/Drake/Summer2020/ATAC-Seq/fimo_post-process/peak-consensus-p-0.00001-CHOL-250.csv")

D1 <- D1[which(!is.na(D1[,1])),]; rownames(D1) <- D1[,1]; D1 <- D1[,-1]
D2 <- D2[which(!is.na(D2[,1])),]; rownames(D2) <- D2[,1]; D2 <- D2[,-1]

cgenes <- intersect(rownames(D1),rownames(D2))
corr <- cor(D1[cgenes,],D2[cgenes,])
sum(D1[cgenes,]-D2[cgenes,])

summary(diag(corr))
