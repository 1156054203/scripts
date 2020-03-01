# help: https://www.yuque.com/docs/share/46b70be0-0444-4825-8488-e9ee7f35a593#

# 软件安装
## if(!require('DESeq2')){  BiocManager::install('DESeq2',ask = F,update = F,  version = "3.8")  }

options(bitmapType='cairo')
## install UMAP by python
# if(!require('reticulate')){  BiocManager::install('reticulate',ask = F,update = F)  }; library(reticulate)
# reticulate::py_config()
# use_python("/rainbow/software/anaconda3/bin/python3.6")
# use_virtualenv("~/myenv")
# reticulate::py_install(packages = 'umap-learn')
## seurat

{
	if(!require('sctransform')){  BiocManager::install('sctransform',ask = F,update = F)  }; library(sctransform)
	if(!require(Cairo)){install.packages("Cairo")};library(Cairo)
	options(bitmapType='cairo')
	if(!require(stringr)){install.packages("stringr")};library(stringr)
	if(!require(ggplot2)){install.packages("ggplot2")};library(ggplot2)
	if(!require(cowplot)){install.packages("cowplot")};library(cowplot)
	## seurat
	if(!require(Seurat)){install.packages("Seurat")};library(Seurat)
	if(!require(pheatmap)){install.packages("pheatmap")};library(pheatmap)
	if(!require(MAST)){install.packages("BiocManager") ;BiocManager::install("MAST")}# Needed to install all Bioconductor packages
	if(!require(dplyr)){install.packages("dplyr")};library(dplyr)
	if(!require(org.Hs.eg.db)){ BiocManager::install('org.Hs.eg.db',ask = F,update = F) }; library(org.Hs.eg.db)
	if(!require(org.Mm.eg.db)){ BiocManager::install('org.Mm.eg.db',ask = F,update = F) }; library(org.Mm.eg.db)
	if(!require(clusterProfiler)){ BiocManager::install('clusterProfiler',ask = F,update = F) };library(clusterProfiler)
}
Outdir<- function(step2){
        if(!file.exists(step2)) {dir.create(step2,recursive=T)}
        setwd(step2)
}


getInitialPNG <- function(seuratObj){
	require(Cairo)
	## view and evaluation the initial data
	mito.genes1 <- grep(pattern = "^mt-",x = rownames(seuratObj), ignore.case = TRUE, value = TRUE)
	mito.genes2 <- grep(pattern = "^mrp",x = rownames(seuratObj), ignore.case = TRUE ,value = TRUE)
	percent.mito <- (Matrix::colSums( GetAssayData(object = seuratObj, slot = "counts")[mito.genes1, ]) + 
		Matrix::colSums( GetAssayData(object = seuratObj, slot = "counts")[mito.genes2, ])) / Matrix::colSums(GetAssayData(object = seuratObj, slot = "counts"))*100
	#add the percent.mito to metagene matrix
	#seuratObj[["percent.mt"]] <- PercentageFeatureSet(object = seuratObj, pattern = "MT-")
	seuratObj <- AddMetaData(object=seuratObj,metadata=percent.mito,col.name="percent.mt")
	
	#VlnPlot the nGene,nUMI and percent.mito
	CairoPNG(file="01.vinPlot.png",width=800,height=500)
	print(VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
	dev.off()
	
	CairoPNG(file="02.corPlot.png",width=1000,height=500)
	plot1 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA",pt.size=0.1, feature2 = "percent.mt")
	plot2 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", pt.size=0.1,feature2 = "nFeature_RNA")
	print(CombinePlots(plots = list(plot1, plot2)))
	dev.off()
	## 
	return(seuratObj)
	cat("getInitialPNG work has done!!!","\n")
	cat("=============================","\n")
}


getFiltObj <- function(seuratObj) {
	require(stringr)
	path_step1 <- getwd()

  #--------------------------------------------
	# save the raw meta_data
	setwd('../')
	metaData<- data.frame(barcode=rownames(seuratObj@meta.data), sample=seuratObj@meta.data$orig.ident, 
		nGene=seuratObj@meta.data$nFeature_RNA, nUMI=seuratObj@meta.data$nCount_RNA, percent.mt=seuratObj@meta.data$percent.mt)
	metaData$barcode <- str_sub(metaData$barcode,start=-16)
	write.csv(metaData,file='raw_meta_data.csv',row.names=F)	
	## save the raw all_genes, matrix
	all_genes <- rownames(x = seuratObj)
	write.csv(data.frame(all_genes), file = "all_genes.csv", row.names = FALSE)
	raw_data <- t(data.frame(FetchData(seuratObj, vars = all_genes))) #cells.use = NULL
	write.table(raw_data, file = "raw_bc_matrix.xls", sep = "\t")
  
	#--------------------------------------------
	## remove low qulality cells and background genes
	## 1-filter data by UMI,Gene and percent.mito 
	setwd(path_step1)
	# nCount_RNA_ranger=quantile(seuratObj$nCount_RNA,c(0.99))
	seuratObj <- subset(x = seuratObj, subset = nFeature_RNA > 200 & percent.mt < 10)

  # seuratObj=seuratObj_tmp
	CairoPNG(file="03.filt_vinPlot.png",width=800,height=500)
	print(VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
	dev.off()
	
	plot1 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	
	CairoPNG(file="04.filt_corPlot.png",width=1000,height=500)
	print(CombinePlots(plots = list(plot1, plot2)))
	dev.off()
#--------------------------------------------
  ## 2-remove ribosome and mitochondrial genes
	setwd('../')
	filtered_genes <- grep(pattern = "^mt-" , x = rownames(x=seuratObj), ignore.case = TRUE, value = TRUE, invert = TRUE)
	filtered_genes <- grep(pattern = "^mrp" , x = filtered_genes, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filtered_genes <- grep(pattern = "^rpl" , x = filtered_genes, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filtered_genes <- grep(pattern = "^rps" , x = filtered_genes, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filtered_genes <- grep(pattern = "^rpf" , x = filtered_genes, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filtered_genes <- grep(pattern = "^rpn" , x = filtered_genes, ignore.case = TRUE, value = TRUE, invert = TRUE)
	write.csv(filtered_genes, file = "filtered_genes.csv", row.names = FALSE )
	## save the filtered matrix
	filtered_genes <- filtered_genes
	filtered_data <- t(data.frame(FetchData(seuratObj, vars = filtered_genes))) #cells.use = NULL
	write.table(filtered_data, file = "filtered_bc_matrix.xls", sep = "\t")
	#filtered_data <- read.table(file = "filtered_bc_matrix.xls", header=T, row.names=1, sep = "\t")
	# using data that filtered low-quality barcodes and mito-ribosome genes 
	###########################################################################################
	## reload the filtered expression matrix and creat seurat object
	setwd(path_step1)
	seuratObj <- CreateSeuratObject(counts = filtered_data, project =  seuratObj@project.name)
	#seuratObj <- AddMetaData(object=seuratObj,metadata=percent.mito,col.name="percent.mt")
	CairoPNG(file="05.final_vinPlot.png",width=550,height=500)
	print(VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2))
	dev.off()
	CairoPNG(file="06.final_corPlot.png",width=600,height=500)
	#plot1 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	print(plot2)
	dev.off()
	############################################################################################
	return(seuratObj)
	cat("getFiltObj1 work has done!!!","\n")
	cat("=============================","\n")
}


RENAME_CELL_TYPE  <- function(seuratObj,pathHome=pathHome,clusterfile=clusterfile) {
	#-----------------------------------
	# add cell type
	#-----------------------------------
	outdir=paste(pathHome,"/004.CellType",sep="")
	if(!file.exists(outdir)) {dir.create(outdir)}	
	setwd(outdir)
	# clusterfile="/Data2/wanglp/scRNASeq/JYNJ1908CYH02_JYNJ1908CYH03_Xiaf_scRNAseq/Seurat/CD45_plus_CellType.txt"
	celltype=read.table(clusterfile,head=T,row.names=1,sep="\t")
	#celltype=read.table("CellType.txt",head=T,row.names=1,sep="\t")
	ident=seuratObj[[]][,"seurat_clusters"]
	ident_new=factor(ident,levels=rownames(celltype),labels=as.character(celltype[,1]))
	Idents(object = seuratObj) <- ident_new
	seuratObj$CellType <- Idents(object = seuratObj)

	cluster_barcode=data.frame(barcode=names(Idents(seuratObj)),id=Idents(seuratObj))
	file_cluster=("cellTypeForCloup.csv")
	write.table(cluster_barcode,file_cluster,quote=F,col.names=T,row.names=F,sep=",")
	TSNE_UMAP_PLOT(seuratObj,outdir,group="CellType")
	return(seuratObj)
}





















