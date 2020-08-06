rm(list=ls())
options(stringsAsFactors = FALSE)

# Using the MEME55 curated Cis-BP56 TF-binding motif reference, we scanned each ATAC-seq tumor type peak atlas and common atlas with FIMO57 to find peaks likely to contain each motif 
# (P < 10−5). We filtered TFs that were not expressed in at least 50% of samples in at least one of the five tumor types. 
library(ChIPpeakAnno)
library(org.Hs.eg.db)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

annoData<-toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene,feature="gene")
cancertpyes <- c("ESCA", "")
myfile <- read.table("/Volumes/SeagateBackupPlusDrive/Drake/Summer2020/ATAC-Seq/FIMO/PAN/fimo.tsv",header = T,stringsAsFactors = F,sep="\t") #[,c(3,4,5,6,10)]

motif2peak <- myfile[which(myfile[,"p.value"] < 10^-5),]
motifs <- sort(unique(motif2peak[,"motif_alt_id"])) ##motifs 731
peaks <- unique(motif2peak[,"sequence_name"])
peak_motif <- matrix(0, ncol=length(motifs), nrow=length(peaks))
rownames(peak_motif) <- peaks
colnames(peak_motif) <- motifs
peak_motif <- as.data.frame(peak_motif)
dim(peak_motif)

for(i in 1:length(peaks)){ #
  motif <- data.frame(table(motif2peak[motif2peak[,"sequence_name"] %in% peaks[i], "motif_alt_id"]))
  tmp <- motif2peak[motif2peak[,"sequence_name"] %in% peaks[i], c("motif_alt_id","score")]
  tmp <- tmp[order(tmp[,"score"], decreasing = T),]
  tmp <- tmp[!duplicated(tmp[,"motif_alt_id"]),]
  peak_motif[peaks[i],as.character(unlist(tmp[,"motif_alt_id"]))] <- tmp[,"score"]
}

write.csv(peak_motif, file="/Volumes/SeagateBackupPlusDrive/Drake/Summer2020/ATAC-Seq/fimo_post-process/peak-consensus-p-0.00001-PAN-250.csv", quote = F)


