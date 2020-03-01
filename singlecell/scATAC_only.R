#help: https://www.yuque.com/docs/share/e8b17963-83bc-4f10-a787-6f6d426fe41d#

source("/Data2/wanglp/bin/Seurat/Self/library/scATACSeq_lib.R")
#----------------------------------------
# Input
#----------------------------------------

outdir="/Data2/wanglp/scATACSeq/JYSH1903LHT27_Zhengyx_2scATAC_20190328/signac_seurat/Peak/"
scATAC_h5="/Data2/wanglp/scATACSeq/JYSH1903LHT27_Zhengyx_2scATAC_20190328/WE_WT/outs/filtered_peak_bc_matrix.h5"
fragment_path="/Data2/wanglp/scATACSeq/JYSH1903LHT27_Zhengyx_2scATAC_20190328/WE_WT/outs/fragments.tsv.gz"
meta_file="/Data2/wanglp/scATACSeq/JYSH1903LHT27_Zhengyx_2scATAC_20190328/WE_WT/outs/singlecell.csv"
sample_name_order=c("WE","WT")
#G_VER=EnsDb.Hsapiens.v86
G_VER=EnsDb.Mmusculus.v75  #V79b版本有问题, 该参数为包名，不能加引号,

#----------------------------------------
# initial data
#----------------------------------------
Step1_peak=paste(outdir,"/001_Peak",sep="")
Outdir(Step1_peak)
metadata <- read.csv(file = meta_file, header = TRUE, row.names = 1)
counts <- Read10X_h5(scATAC_h5)
pbmc <- CreateSeuratObject(counts = counts, project = 'ATAC',assay = "ATAC", min.cells = 1, meta.data = metadata)
# [1] 82889 11501
pbmc <- SetFragments( object = pbmc, file = fragment_path)
cat(" [] Get cell Count in rawdata:",dim(pbmc)[2] ,"\n")

#-----------------------------------------------------
# QC
#-----------------------------------------------------
pbmc <- NucleosomeSignal(object = pbmc)

# plot: VlnPlot
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
PLOT_QC(pbmc,"celranger",Step1_peak)

# plot: fragment_size_distribution
NS=10
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > NS, paste("'NS > ",NS,"'",sep=""), paste("'NS < ",NS,"'",sep=""))
PLOT_FRAG(pbmc,Step1_peak)

#filter out cells 
# peak_frag_range=RANGER(pbmc)
# pbmc.atac <- subset(pbmc.atac, subset = nCount_ATAC > 500)
pbmc_filter <- subset(pbmc, subset = peak_region_fragments > 500 & 
         pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < NS)

PLOT_QC(pbmc_filter, "filter", Step1_peak)
cat(" [] Get cell rate:",dim(pbmc_filter)[2]/dim(pbmc)[2], " Count:",dim(pbmc_filter)[2],"\n")

#-----------------------------------------------------
# Normalization and linear dimensional reduction
#-----------------------------------------------------
DefaultAssay(pbmc_norm)

pbmc_norm <- RunTFIDF(pbmc_filter) # save data to slot="data"
pbmc_norm=ScaleData(pbmc_norm)  # save data to slot="scale.data"
pbmc_norm <- FindTopFeatures(pbmc_norm, min.cutoff = 'q70') # 此时用于降维的 Features 自动设置为 VariableFeatures
pbmc_norm <- RunSVD( object = pbmc_norm, assay = 'ATAC', reduction.key = 'LSI_', reduction.name = 'lsi')

#-----------------------------------------------------
# Non-linear dimension reduction and clustering
#-----------------------------------------------------

pbmc_red <- RunUMAP(object = pbmc_norm, reduction = 'lsi', dims = 1:50)
pbmc_red <- RunTSNE(object = pbmc_red, reduction = 'lsi', dims = 1:50)
pbmc_red <- FindNeighbors(object = pbmc_red, reduction = 'lsi', dims = 1:50)
pbmc_red = ADD_SAMPLE(sample_name_order,pbmc_red)

for(resolution in 1:5*0.2){
  pbmc_red <- FindClusters(object = pbmc_red, verbose = FALSE,resolution=resolution)
  table(pbmc_red$seurat_clusters)
  PLOT_red(pbmc_red,Step1_peak,resolution=resolution)
}

saveRDS(pbmc_red, file = paste(Step1_peak,"/scATAC.single.rds",sep=""))
# pbmc_red=readRDS(paste(Step1_peak,"/scATAC.single.rds",sep=""))

#-----------------------------------------------------
# DEAnalysis_Peak
#-----------------------------------------------------
# DE peaks
markerFind_scRNA_scATACSeq(pbmc_red,Step1_peak)
# markerFind_scRNA_scATACSeq(pbmc_red,Step1_peak,calculate="no")

#-----------------------------------------------------
# add activity，并进行标准化和scale，差异分析
#-----------------------------------------------------
# EnsDb.Mmusculus.v75 版本有问题，没有feature 信息
for(region in c("activity","genebody","promoter")){
	pbmc_activity=pbmc_red
	if(region == "activity"){
		Step1=paste(outdir,"/001_ACTIVITY",sep="")
		gene_activities = EXT_COORDS(G_VER,fragment_path,pbmc_activity,region) # 第一个参数时包名，不能加引号, ~20min , 1] 21893  9541
		pbmc_activity[['ACTIVITY']] <- CreateAssayObject(counts = gene_activities)
		
	}
	if(region == "genebody"){
		Step1=paste(outdir,"/001_genebody",sep="")
		gene_activities = EXT_COORDS(G_VER,fragment_path,pbmc_activity,region) # 第一个参数时包名，不能加引号, ~20min
		pbmc_activity[['ACTIVITY']] <- CreateAssayObject(counts = gene_region$BODY)
	}
	if(region == "promoter"){
		Step1=paste(outdir,"/001_promoter",sep="")
		gene_activities = EXT_COORDS(G_VER,fragment_path,pbmc_activity,region) # 第一个参数时包名，不能加引号, ~20min
		pbmc_activity[['ACTIVITY']] <- CreateAssayObject(counts = gene_region$PROMOTER)
	}

	Outdir(Step1)
	DefaultAssay(pbmc_activity)="ACTIVITY"
	vlnPlot_corPlot_tmp(Step1,pbmc_activity)

	pbmc_norm <- RunTFIDF(pbmc_activity) # save data to slot="data"
	pbmc_norm=ScaleData(pbmc_norm)  # save data to slot="scale.data"

	saveRDS(pbmc_norm, file = paste(Step1,"/scATAC.Activity.rds",sep=""))
	# pbmc_norm=readRDS(paste(Step1,"/scATAC.Activity.rds",sep=""))

	markerFind_scRNA_scATACSeq(pbmc_norm,Step1)
}