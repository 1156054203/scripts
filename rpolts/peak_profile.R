source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg38")


afrom <- 46128000
ato <- 46146799
library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
axTrack <- GenomeAxisTrack(col="black",fontcolor="black")
idxTrack <- IdeogramTrack(genome = "hg38", chromosome = "chr6",col="black",fontcolor="black")
knownGenes <- UcscTrack(genome = "hg38", chromosome = "chr6",track = "knownGene", from = afrom, to = ato, trackType = "GeneRegionTrack",rstarts = "exonStarts", rends = "exonEnds", gene = "name",symbol = "name", transcript = "name", strand = "strand",fill = "#8282d2", name = "ENPP4",col="black",fontcolor="black")

WT_1 <- file.path("/Users/rainbow2/Desktop/zhaomm/zhaomm_project/ATAC-seq/已结/韩蒙蒙/resultReport/5.Macs2PeaksCalling/WT_1.bw")
WT_2 <- file.path("/Users/rainbow2/Desktop/zhaomm/zhaomm_project/ATAC-seq/已结/韩蒙蒙/resultReport/5.Macs2PeaksCalling/WT_2.bw")
KO_1 <- file.path("/Users/rainbow2/Desktop/zhaomm/zhaomm_project/ATAC-seq/已结/韩蒙蒙/resultReport/5.Macs2PeaksCalling/KO_1.bw")
KO_2 <- file.path("/Users/rainbow2/Desktop/zhaomm/zhaomm_project/ATAC-seq/已结/韩蒙蒙/resultReport/5.Macs2PeaksCalling/KO_2.bw")

alTrack1 <- DataTrack(range = WT_1, start = afrom, end = ato, chromosome = "chr1", strand = "+",gneome = "hg38",name = "WT_1",stream = FALSE,type="histogram",ylim = c(0, 1.5),col.histogram = "darkblue",fill.histogram = "darkblue")
alTrack2 <- DataTrack(range = WT_2, start = afrom, end = ato, chromosome = "chr1",strand = "+",gneome = "hg38",name = "WT_2",stream = FALSE,type="histogram",ylim = c(0, 1.5),col.histogram = "darkblue",fill.histogram = "darkblue")
alTrack3 <- DataTrack(range = KO_1, start = afrom, end = ato, chromosome = "chr1", strand = "+",gneome = "hg38",name = "KO_1",stream = FALSE,type="histogram",ylim = c(0, 1.5),col.histogram = "darkblue",fill.histogram = "darkblue")
alTrack4 <- DataTrack(range = KO_2, start = afrom, end = ato, chromosome = "chr1",strand = "+",gneome = "hg38",name = "KO_2",stream = FALSE,type="histogram",ylim = c(0, 1.5),col.histogram = "darkblue",fill.histogram = "darkblue")

bedfile <- file.path("/Users/rainbow2/Desktop/zhaomm/zhaomm_project/ATAC-seq/已结/韩蒙蒙/resultReport/13.DiffPeaks/TreatmentVSControl_up.bed")
bed <-AnnotationTrack(range = bedfile,start = afrom, end = ato,chromosome = "chr6", strand = "+",gneome = "hg38",name = "ENPP4",stream = FALSE)

plotTracks(list(idxTrack, axTrack,knownGenes,alTrack1,alTrack2,alTrack3,alTrack4,bed),panel.only = FALSE,from = afrom, to = ato,background.panel = "white", background.title = "white",col.title="black",col.axis="black",title.width=1.5,cex.title=0.8,cex.axis=0.8)

