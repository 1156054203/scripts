#!/online/software/R-3.1.2/bin/Rscript
.libPaths('/online/software/R-3.1.2/library')

args <- commandArgs(trailingOnly = TRUE)
cir_path <- "/online/home/chenyl/circos/scripts"


samplename <- args[1]
assembly <- args[2]

#samplename <- gsub("^\\s+|\\s+$", "", args[1])
#assembly <- ""
#if(length(args) == 3) {
#  assembly <- gsub("^\\s+|\\s+$", "", args[2])
#} else {
#  assembly <- "hg19"
#}

setwd(cir_path)

setEPS()
postscript(paste(samplename, ".circos", ".eps", sep=""))
jpeg(paste(samplename, ".jpg",sep=''),width = 500, height = 500, units = "px")
par(mar = c(1, 1, 1, 1), lwd = 0.5)
library("circlize")
if(assembly == "hg19") {
  circos.initializeWithIdeogram()
} else {
  circos.initializeWithIdeogram(cytoband="hg38_cytoBand.txt")
}

# Nonsynonymous SNV
tvbed <- read.table(paste(samplename, ".tv.bed", sep=""), head=T)
tibed <- read.table(paste(samplename, ".ti.bed", sep=""), head=T)
bed_list <- list(tvbed, tibed)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF0000", "#0000FF"))

# CNV Pvalue < 1e-5 & size > 300K
gainbed <- read.table(paste(samplename, ".cnv.gain.bed", sep=""), head=T)
circos.genomicTrackPlotRegion(gainbed, stack = TRUE, panel.fun = function(region, value, ...) { circos.genomicRect(region, value, col = "red", border = "red") }, bg.border = NA, track.height = 0.05)
lossbed <- read.table(paste(samplename, ".cnv.loss.bed", sep=""), head=T)
circos.genomicTrackPlotRegion(lossbed, stack = TRUE, panel.fun = function(region, value, ...) { circos.genomicRect(region, value, col = "green", border = "green") }, bg.border = NA, track.height = 0.05)

# SV
bed_list = read.table(paste(samplename, ".sv.bed", sep=""), head=T)
bed1 = bed_list[,1:3]
bed2 = bed_list[,4:6]
circos.genomicLink(bed1, bed2, col=sample(10, nrow(bed1), replace=TRUE), border=sample(10, nrow(bed1), replace=TRUE))


# Close output
dev.off()
