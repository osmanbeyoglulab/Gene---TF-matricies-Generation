# to read txt files
library(readr)

# to transform data into GenomicRanges
library(GenomicRanges)

# other ones used to prepare the data
library(tidyr)
library(dplyr)
library(SummarizedExperiment)

# For the t.test loop
library(plyr)

# For easy volcano plot  *Not able to load
library(TCGAbiolinks)

# For heatmap plot
library(ComplexHeatmap)
library(circlize)

# For the bigwig plot
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Other libraries 
library(GenomicAlignments)
library(data.table)
library(Rsamtools)
library(rtracklayer)
library(ChIPpeakAnno)
library(stringr)
library('BSgenome.Hsapiens.UCSC.hg38')

# Load the data in for specific and pancancers #####
setwd("/Volumes/SeagateBackupPlusDrive/Drake/Summer2020/ATAC-Seq/TCGA-ATAC_Cancer_Type-specific_PeakCalls")
atac_acc <- readr::read_tsv("ACC_peakCalls.txt")
atac_blca <- readr::read_tsv("BLCA_peakCalls.txt")
atac_brca <- readr::read_tsv("BRCA_peakCalls.txt")
atac_cesc <- readr::read_tsv("CESC_peakCalls.txt")
atac_chol <- readr::read_tsv("CHOL_peakCalls.txt")
atac_coad <- readr::read_tsv("COAD_peakCalls.txt")
atac_esca <- readr::read_tsv("ESCA_peakCalls.txt")
atac_gbm <- readr::read_tsv("GBM_peakCalls.txt")
atac_hnsc <- readr::read_tsv("HNSC_peakCalls.txt")
atac_kirc <- readr::read_tsv("KIRC_peakCalls.txt")
atac_kirp <- readr::read_tsv("KIRP_peakCalls.txt")
atac_lgg <- readr::read_tsv("LGG_peakCalls.txt")
atac_lihc <- readr::read_tsv("LIHC_peakCalls.txt")
atac_luad <- readr::read_tsv("LUAD_peakCalls.txt")
atac_lusc <- readr::read_tsv("LUSC_peakCalls.txt")
atac_meso <- readr::read_tsv("MESO_peakCalls.txt")
atac_pcpg <- readr::read_tsv("PCPG_peakCalls.txt")
atac_prad <- readr::read_tsv("PRAD_peakCalls.txt")
atac_skcm <- readr::read_tsv("SKCM_peakCalls.txt")
atac_stad <- readr::read_tsv("STAD_peakCalls.txt")
atac_tgct <- readr::read_tsv("TGCT_peakCalls.txt")
atac_thca <- readr::read_tsv("THCA_peakCalls.txt")
atac_ucec <- readr::read_tsv("UCEC_peakCalls.txt")
atac_pan <- readr::read_tsv("/Volumes/SeagateBackupPlusDrive/Drake/Summer2020/ATAC-Seq/TCGA-ATAC_PanCancer_PeakSet.txt")

# Filter peaks from pan
specific.cancers <- c("BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "UCEC")
pan_indexes <- c()
for (i in 1:length(specific.cancers)) {
 pan_indexes <- append(pan_indexes, which(grepl(specific.cancers[i], atac_pan$name, fixed=T) == TRUE))
}
atac_pan <- atac_pan[pan_indexes,]

# Transform Peaks into GRanges #####
make_gr <- function(cancer_peaks) {
  gr <- makeGRangesFromDataFrame(cancer_peaks,keep.extra.columns = T)
  names(gr) <- paste0(as.character(seqnames(gr)), "_", start(gr), "-", end(gr))
  return(gr)
}

atac_acc.gr <- make_gr(atac_acc)
atac_blca.gr <- make_gr(atac_blca)
atac_brca.gr <- make_gr(atac_brca)
atac_cesc.gr <- make_gr(atac_cesc)
atac_chol.gr <- make_gr(atac_chol)
atac_coad.gr <- make_gr(atac_coad)
atac_esca.gr <- make_gr(atac_esca)
atac_gbm.gr <- make_gr(atac_gbm)
atac_hnsc.gr <- make_gr(atac_hnsc)
atac_kirc.gr <- make_gr(atac_kirc)
atac_kirp.gr <- make_gr(atac_kirp)
atac_lgg.gr <- make_gr(atac_lgg)
atac_lihc.gr <- make_gr(atac_lihc)
atac_luad.gr <- make_gr(atac_luad)
atac_lusc.gr <- make_gr(atac_lusc)
atac_meso.gr <- make_gr(atac_meso)
atac_pcpg.gr <- make_gr(atac_pcpg)
atac_prad.gr <- make_gr(atac_prad)
atac_skcm.gr <- make_gr(atac_skcm)
atac_stad.gr <- make_gr(atac_stad)
atac_tgct.gr <- make_gr(atac_tgct)
atac_thca.gr <- make_gr(atac_thca)
atac_ucec.gr <- make_gr(atac_ucec)
atac_pan.gr <- make_gr(atac_pan)

# Quality check #####
peaks_in_pancancer <- function(atac_specific) {
  length(subsetByOverlaps(makeGRangesFromDataFrame(atac_specific,keep.extra.columns = T),atac_pan.gr))
}


# Get sequences for each GRanges ####
GRanges_to_fasta <- function(cancer_type, specific.gr) {
  # nonflanks
  setwd("/Volumes/SeagateBackupPlusDrive/Drake/Summer2020/ATAC-Seq")
  dir.create("nonflanking_fasta/")
  setwd("nonflanking_fasta/")
  nseq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, specific.gr)
  writeXStringSet(nseq, paste0(cancer_type, ".fasta"))
  
  # flanks
  setwd("/Volumes/SeagateBackupPlusDrive/Drake/Summer2020/ATAC-Seq")
  dir.create("250_flanks_fasta/")
  setwd("250_flanks_fasta/")
  flank_left<-flank(specific.gr,width = 250)
  flank_right<-flank(specific.gr,width = 250,start = F)
  flank<-c(flank_left,flank_right)
  flank_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, flank)
  writeXStringSet(flank_seq, paste0(cancer_type,"_flank.fasta"))
}

GRanges_to_fasta("ACC", atac_acc.gr)
GRanges_to_fasta("BLCA", atac_blca.gr)
GRanges_to_fasta("BRCA", atac_brca.gr)
GRanges_to_fasta("CESC", atac_cesc.gr)
GRanges_to_fasta("CHOL", atac_chol.gr)
GRanges_to_fasta("COAD", atac_coad.gr)
GRanges_to_fasta("ESCA", atac_esca.gr)
GRanges_to_fasta("GBM", atac_gbm.gr)
GRanges_to_fasta("HNSC", atac_hnsc.gr)
GRanges_to_fasta("KIRC", atac_kirc.gr)
GRanges_to_fasta("KIRP", atac_kirp.gr)
GRanges_to_fasta("LGG", atac_lgg.gr)
GRanges_to_fasta("LIHC", atac_lihc.gr)
GRanges_to_fasta("LUAD", atac_luad.gr)
GRanges_to_fasta("LUSC", atac_lusc.gr)
GRanges_to_fasta("MESO", atac_meso.gr)
GRanges_to_fasta("PAN", atac_pan.gr)
GRanges_to_fasta("PCPG", atac_pcpg.gr)
GRanges_to_fasta("PRAD", atac_prad.gr)
GRanges_to_fasta("SKCM", atac_skcm.gr)
GRanges_to_fasta("STAD", atac_stad.gr)
GRanges_to_fasta("TGCT", atac_tgct.gr)
GRanges_to_fasta("THCA", atac_thca.gr)
GRanges_to_fasta("UCEC", atac_ucec.gr)

# Run FIMO 
