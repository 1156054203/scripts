source("library/Intergrate_scRNA_scATAC_lib.R")
#----------------------------------------
# Input
#----------------------------------------
outdir="002_Activity_combined_scRNASeq"
scATAC_rds="scATAC.Activity.rds" # Activity
scRNASeq_rds="sampleCombined_allCellType.rds" #cluster 0-19
markerFile="marker_custom.csv"



#----------------------------------------
# combined
#----------------------------------------
pbmc_atac=readRDS(scATAC_rds)
pbmc_atac$tech <- "atac"
#pbmc.atac <- FindClusters(pbmc.atac)
#DefaultAssay(pbmc_atac) <- "ACTIVITY"

pbmc_rna <- readRDS(scRNASeq_rds)
pbmc_rna$tech <- "rna"
#DefaultAssay(pbmc_rna)='integrated'

Outdir(outdir)

#----------------------------------------
# Identify anchors 
#----------------------------------------

transfer.anchors <- FindTransferAnchors(reference = pbmc_rna, query = pbmc_atac, features = VariableFeatures(object = pbmc_rna), reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype_predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc_rna$CellType, weight.reduction = pbmc_atac[["lsi"]])
pbmc_atac <- AddMetaData(pbmc_atac, metadata = celltype_predictions)

#----------------------------------------
# prediction.score
file1_png=paste(outdir,"/prediction.score.max.png",sep="")
CairoPNG(file=file1_png,width=480*2,height=480)
hist(pbmc_atac$prediction.score.max)
abline(v = 0.5, col = "red")
dev.off()

saveRDS(pbmc_atac, file = paste(outdir,"/scRNA_scATAC.combind.rds",sep=""))
# pbmc_atac=readRDS(paste(outdir,"/scRNA_scATAC.combind.rds",sep=""))

#----------------------------------------
#  choose results whose score > 0.5, and plot
score=0.5
table(pbmc_atac$prediction.score.max > score)
pbmc_atac_filtered <- subset(pbmc_atac, subset = prediction.score.max > score)
pbmc_atac_filtered$predicted.id <- factor(pbmc_atac_filtered$predicted.id, levels = levels(pbmc_rna))  # to make the colors match

TSNE_UMAP_Plot_CellType_tmp(outdir,pbmc_atac_filtered,pbmc_rna)

saveRDS(pbmc_atac_filtered, file = paste(outdir,"/scRNA_scATAC.combind.highscore.rds",sep=""))
# pbmc_atac_filtered=readRDS(paste(outdir,"/scRNA_scATAC.combind.highscore.rds",sep=""))

# save cell count and ratio in each sample
count_cell=table(pbmc_atac_filtered[[]][,c("orig.ident","predicted.id")])
write.csv(t(count_cell),file=paste(outdir,'/cellNum_InCluster.csv',sep=""),row.names=T)
count_cell2 <- read.csv(file=paste(outdir,'/cellNum_InCluster.csv',sep=""),header=T,row.names=1)
ratio=round(100*t(t(count_cell2)/apply(count_cell2,2,sum)),3)
data=cbind(t(count_cell),Total=apply(t(count_cell),1,sum),ratio)
write.csv(data,file=paste(outdir,'/cellNum_InCluster.csv',sep=""),row.names=T)

#----------------------------------------
# 差异分析
#----------------------------------------
# DE genes
DefaultAssay(pbmc_atac_filtered) <- "ACTIVITY"
# markerFind_scRNA_scATACSeq(pbmc_atac_filtered,outdir,,calculate="no")

# plot1, list for known gene with info. of  cell type
markerList <- read.csv(file=markerFile, header=T)
Outdir(paste(outdir,"/known",sep=""))
markerPlot(pbmc_atac_filtered, marker_list=markerList,group="predicted.id",QUICK=TRUE)

# plot2,list for DE gene
markerFile_tmp=paste(outdir,"/top10_list.csv",sep="")
markerList <- read.csv(file=markerFile_tmp, header=T)
Outdir(paste(outdir,"/DEGs",sep=""))
markerPlot(pbmc_atac_filtered, marker_list=markerList,group="predicted.id",QUICK=TRUE)

# plot3, single genes
Outdir(paste(outdir,"/specific_marker",sep=""))
markerPlot_gene_name("CD3g",pbmc_atac_filtered,"predicted.id")
markerPlot_gene_name("CD4",pbmc_atac_filtered,"predicted.id")
markerPlot_gene_name("CD8a",pbmc_atac_filtered,"predicted.id")

#plot4, coverage info
COVERAGE_PLOT(pbmc.atac.filtered,G_VER=EnsDb.Mmusculus.v75,'Cd8a','chr6_71373427_71379173')


