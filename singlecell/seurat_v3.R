#############################################################################################################
#seurat functions including: getInitialPNG, getFiltObj, getOCA_TSNE, getGraphBased, markerPlot, markerFind
#2019
#Seurat v3.0
#############################################################################################################
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
	plot1 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	print(CombinePlots(plots = list(plot1, plot2)))
	dev.off()
	## 
	return(seuratObj)
	cat("getInitialPNG work has done!!!","\n")
	cat("=============================","\n")
}
#############################################################################################################
getFiltObj <- function(seuratObj){
	require(stringr)
	path_step1 <- getwd()
	#percent.mito <- seuratObj@meta.data$percent.mt
	############################################################
	## save the raw meta_data
	setwd('../')
	metaData<- data.frame(barcode=rownames(seuratObj@meta.data), sample=seuratObj@meta.data$orig.ident, 
		nGene=seuratObj@meta.data$nFeature_RNA, nUMI=seuratObj@meta.data$nCount_RNA, percent.mt=seuratObj@meta.data$percent.mt)
	metaData$barcode <- str_sub(metaData$barcode,start=-3)
	write.csv(metaData,file='raw_meta_data.csv',row.names=F)	
	## save the raw all_genes, matrix
	all_genes <- rownames(x = seuratObj)
	write.csv(all_genes, file = "all_genes.csv", row.names = FALSE, col.names = FALSE)
	raw_data <- t(data.frame(FetchData(seuratObj, vars = all_genes))) #cells.use = NULL
	write.table(raw_data, file = "raw_bc_matrix.xls", sep = "\t")
	############################################################
	## remove low qulality cells and background genes
	## 1-filter data by UMI,Gene and percent.mito 
	setwd(path_step1)
	seuratObj <- subset(x = seuratObj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
	CairoPNG(file="03.filt_vinPlot.png",width=800,height=500)
	print(VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
	dev.off()
	CairoPNG(file="04.filt_corPlot.png",width=1000,height=500)
	plot1 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	print(CombinePlots(plots = list(plot1, plot2)))
	dev.off()
	## 2-remove ribosome and mitochondrial genes
	setwd('../')
	filtered_genes <- grep(pattern = "^mt-" , x = rownames(x=seuratObj), ignore.case = TRUE, value = TRUE, invert = TRUE)
	filtered_genes <- grep(pattern = "^mrp" , x = filtered_genes, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filtered_genes <- grep(pattern = "^rpl" , x = filtered_genes, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filtered_genes <- grep(pattern = "^rps" , x = filtered_genes, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filtered_genes <- grep(pattern = "^rpf" , x = filtered_genes, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filtered_genes <- grep(pattern = "^rpn" , x = filtered_genes, ignore.case = TRUE, value = TRUE, invert = TRUE)
	write.csv(filtered_genes, file = "filtered_genes.csv", row.names = FALSE, col.names = FALSE)
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
#############################################################################################################
getStandardization <- function(seuratObj){
	# # data standardization
	# seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
	# #FindVariableGenes
	# seuratObj <- FindVariableFeatures(object = seuratObj, selection.method='vst', nFeatures=2000, mean.cutoff = c(0.1, 10), dispersion.cutoff = c(1, Inf),)
	# # identify the 10 most highly variable genes
	# top10 <- head(x = VariableFeatures(object = seuratObj), 10)
	# # plot variable features with and without labels
	# plot1 <- VariableFeaturePlot(object = seuratObj)
	# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	# CairoPNG(file="07.VarGenes_top10.png",width=800,height=600)
	# print(plot2)
	# dev.off()
	# # center, scale and regress
	# seuratObj <- ScaleData(object = seuratObj, do.scale = TRUE, do.center = TRUE)
	# run sctransform
	require('sctransform')
	seuratObj <- SCTransform(object = seuratObj, do.scale = TRUE, do.center = TRUE, verbose = FALSE)
	#seuratObj <- SCTransform(object = seuratObj, vars.to.regress = "percent.mt", verbose = FALSE)
	# identify the 10 most highly variable genes
	top10 <- head(x = VariableFeatures(object = seuratObj), 10)
	# plot variable features with and without labels
	plot1 <- VariableFeaturePlot(object = seuratObj)
	plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	CairoPNG(file="07.VarGenes_top10.png",width=800,height=600)
	print(plot2)
	dev.off()
	return(seuratObj)
	cat("getStandardization work have done!!!","\n")
	cat("====================================","\n")
}
#############################################################################################################
getPCA_TSNE <- function(seuratObj,PCA_dims=50){
	# PCA reduce the dimension
	# do PCA
	seuratObj <- RunPCA(object = seuratObj, npcs = PCA_dims, ndims.print = 1:3, nfeatures.print = 3)
	# Visualization: distribution
	CairoPNG(file="08.PCA.png",width=700,height=600)
	print(DimPlot(object = seuratObj, reduction = "pca"))
	dev.off()
	#PCs 2
	CairoPNG(file="09.PC1-2.png",width=700,height=600)
	print(VizDimLoadings(object = seuratObj, dims = 1:2, reduction = "pca"))
	dev.off()
	#Determine the ‘dimensionality’ of the dataset
	# PCs elbow line
	CairoPNG(file="10.Elbow.png",width=600,height=400)
	print(ElbowPlot(object = seuratObj, ndims = PCA_dims, reduction = "pca"))
	dev.off()
	# t-SNE plot
	#seuratObj <- RunTSNE(object = seuratObj, dims.use = 1:PCA_dims, do.fast = TRUE, reduction.use = "pca")
	seuratObj <-RunTSNE(seuratObj,reduction='pca',dims=1:PCA_dims,tsne.method='Rtsne',features=rownames(seuratObj))
	p3 <- DimPlot(object = seuratObj, reduction = "tsne", label = TRUE)
	CairoPNG(file='11.combined_tsne.png',width=700,height=600)
	print(p3)
	dev.off()

	# UMAP plot
	seuratObj <- RunUMAP(object = seuratObj, dims = 1:PCA_dims, verbose = FALSE)
	CairoPNG(file="12.combined_UMAP.png",width=700,height=600)
	print(DimPlot(object = seuratObj, reduction = 'umap'))
	dev.off()
	return(seuratObj)
	cat("getPCA_TSNE work have done!!!","\n")
	cat("==============================","\n")
}
#############################################################################################################
getGraphBased <- function(seuratObj,PCA_dims=50,resolution=0.8){
	## We find that setting this parameter between 0.4-1.2 typically 
	## returns good results for single-cell datasets of around 3K cells. 
	## Optimal resolution often increases for larger datasets.
	seuratObj <- FindNeighbors(object = seuratObj, dims = 1:PCA_dims)
	seuratObj <- FindClusters(object = seuratObj, resolution = resolution)
	## tsne plot
	CairoPNG(file="Graph_Based_Clusters_tsne.png",width=700,height=600)
	print(DimPlot(object = seuratObj,  reduction = "tsne", label = TRUE) + NoLegend())
	dev.off()
	# UMAP plot
	CairoPNG(file="Graph_Based_Clusters_UMAP.png",width=700,height=600)
	print(DimPlot(object = seuratObj,  reduction = "umap", label = TRUE) + NoLegend())
	dev.off()
	
	##########################################################
	# export metaData
	# res <- grep(pattern=paste('SCT_snn_res',resolution,sep='.'),x=names(seuratObj@meta.data))

	if(length(levels(as.factor(seuratObj$orig.ident))) > 1){
		metaData<- data.frame(barcode=rownames(seuratObj@meta.data),sample=seuratObj@meta.data$orig.ident,group=seuratObj@meta.data$stim,
			nGene=seuratObj@meta.data$nFeature_RNA,nUMI=seuratObj@meta.data$nCount_RNA,cluster=seuratObj@meta.data$seurat_clusters)
		metaData$barcode <- str_sub(metaData$barcode,end=16)
		write.csv(metaData,file='filtered_meta_data.csv',row.names=F)
		# count the cell number in subCluster and sample
		count_cell1 <- table(metaData[,c(2,6)])
		write.csv(t(count_cell1),file='cellNum_InCluster.csv',row.names=T)
		count_cell1 <- read.csv(file='cellNum_InCluster.csv',header=T)
		count_cell1$SUM <-rowSums(count_cell1[,2:ncol(count_cell1)])
		# count the cell number in subCluster and group
		count_cell2 <- table(metaData[,c(3,6)])
		write.csv(t(count_cell2),file='cellNum_InCluster.csv',row.names=T)
		count_cell2 <- read.csv(file='cellNum_InCluster.csv',header=T)
		##### combine
		count_all <- cbind(count_cell1,count_cell2[,c(2:ncol(count_cell2))])
		colnames(count_all)[1]='clusterID'
		## save
		write.csv(count_all,file='cellNum_InCluster.csv',row.names=F)
	}else{
		metaData<- data.frame(barcode=rownames(seuratObj@meta.data),sample=seuratObj@meta.data$orig.ident,
		nGene=seuratObj@meta.data$nFeature_RNA,nUMI=seuratObj@meta.data$nCount_RNA,cluster=seuratObj@meta.data$seurat_clusters)
		metaData$barcode <- str_sub(metaData$barcode,end=16)
		write.csv(metaData,file='filtered_meta_data.csv',row.names=F)
		# count the cell number in subCluster and sample
		count_cell <- table(metaData[,c(2,5)])
		write.csv(t(count_cell),file='cellNum_InCluster.csv',row.names=T)
		count_cell <- read.csv(file='cellNum_InCluster.csv',header=T)
		colnames(count_cell)[1]='clusterID'
		## save
		write.csv(count_cell,file='cellNum_InCluster.csv',row.names=F)
	}
	return(seuratObj)
	cat("getGraphBased work has done!!!","\n")
	cat("===============================","\n")
}
#############################################################################################################
markerPlot <- function(seuratObj,marker_list){
	##introduction for parameters
	path_save <- getwd()
	colN <- 1:ncol(marker_list)
	marker_type <- colnames(marker_list)
	#initial the value
	gene_list <- data.frame()
	## start ploting
	#DefaultAssay(object = seuratObj) <- "RNA"
	for (counterCol in colN){
		dir.create(marker_type[counterCol])
		path.loop <- paste(path_save,marker_type[counterCol],sep="/")
		setwd(path.loop)
		rowN <- 1:length(marker_list[,counterCol])
		for(counterRow in rowN){
			gene_name <- as.character(marker_list[counterRow,counterCol])
			if(nchar(gene_name)>0){
				if( length( grep(pattern=paste("^",gene_name,"$",sep=""),x=rownames(x=seuratObj),ignore.case = TRUE) ) ){
					gene_name <- grep(pattern=paste("^",gene_name,"$",sep=""),x=rownames(x=seuratObj),ignore.case = TRUE, value=TRUE)
					vinPlot_name <- paste(gene_name,"vin.png",sep="_")
					#tsnePlot_name <- paste(gene_name,"tsne.png",sep="_")
					umapPlot_name <- paste(gene_name,"umap.png",sep="_")
					CairoPNG(file=umapPlot_name,width=700,height=600)
					print(FeaturePlot(object = seuratObj, features = gene_name, reduction = "umap")	)	
					dev.off()
					CairoPNG(file=vinPlot_name,width=700,height=600)
					print(VlnPlot(object = seuratObj, features = gene_name, slot = "counts", log = TRUE))
					dev.off()
					#CairoPNG(file=tsnePlot_name,width=700,height=600)
					#print(FeaturePlot(object = seuratObj, features = gene_name, reduction = "tsne"))
					#dev.off()
					# get the gene to plot dotmap
					gene_list <- c(gene_list, gene_name)
				}
            	else{
					write.table(gene_name, file = "unFoundMarkers.txt", append = TRUE, col.names = FALSE, row.names=FALSE, sep = "\n")
				}
			}
			else next
		}
		setwd(path_save)
	}
	## summary gene dot plot
	gene_list <- unique(as.character(gene_list))[1:20]
	write.csv(gene_list,file='gene_plot.csv',row.names=FALSE)
	CairoPNG(file = "DotPlot_allMarker.png",width=200+50*length(gene_list),height=300+50*length(levels(seuratObj$seurat_clusters)))
	print(DotPlot(seuratObj, features = rev(gene_list), cols = c("blue","red"), dot.scale = 8))
	dev.off()
	if(length(levels(as.factor(seuratObj$orig.ident)))>1) {
		CairoPNG(file = "DotPlot_allSample_allMarker.png",width=200+50*length(gene_list),height=500+100*( length(levels(seuratObj$seurat_clusters))+length(levels(as.factor(seuratObj$orig.ident))) ) )
		print(DotPlot(seuratObj, features = rev(gene_list), cols = c("blue","red","green","purple","coral","gold2","deepskyblue","darkorange"), dot.scale = 8, split.by = "stim") + RotatedAxis())
		dev.off()
	}
	cat("markerPlot work has done!!!","\n")
	cat("=============================","\n")


}
#############################################################################################################
markerFind <- function(seuratObj,test_method="MAST",min.pct=0.25,logFC=0.25){
	if( min.pct>=0 && min.pct<=1 ){
		# ClusterNum <- levels(seuratObj$seurat_clusters)
		# for (numC in ClusterNum ){
		# 	name.cluster = paste("cluster",numC,test_method,"markers.csv",sep="_")
		#   	cluster.markers <- FindMarkers(object = seuratObj, ident.1 = numC, test.use = test_method, min.pct = min.pct, logfc.threshold=logFC)
		#   	write.csv(cluster.markers, file = name.cluster)
		# }
		# find markers for every cluster compared to all remaining cells
		# only the positive ones
		seurat.markers <- FindAllMarkers(object = seuratObj, only.pos = FALSE, test.use = test_method, min.pct = min.pct, logfc.threshold=logFC)
		write.csv(seurat.markers, file = "all_diff_exp_genes.csv",row.names=F)
		#seurat.markers <- read.csv(file='AllMarkers.csv',header=T,row.names=1)
		# plot top gene heatmap
		seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
		#top10 <- seurat.markers %>% group_by(cluster) %>% arrange((p_val_adj),desc(avg_logFC)) %>% top_n(n = 10, wt = avg_logFC)
		CairoPNG(file="Top10Heatmap.png",width=100*length(levels(seuratObj$seurat_clusters)),height=100*length(levels(seuratObj$seurat_clusters)))
		print(DoHeatmap(object = seuratObj, features = top10$gene) + NoLegend())
		dev.off()
	}
	else{
		print("Error!")
		print("Please check input data. The input data may out of range or with incorrect format!")
	}

	#return(seuratObj)
	cat("markerFind work has done!!!","\n")
	cat("=============================","\n")
}
#############################################################################################################
#############################################################################################################
#monocle:analysis pesudotime using seurat result. 
#2019
#############################################################################################################
newimport <- function(otherCDS, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- otherCDS@assays$RNA@counts

    if(class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }

    pd <- tryCatch( {
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, 
    #warning = function(w) { },
    error = function(e) { 
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)

      message("This Seurat object doesn't provide any meta data");
      pd
    })

    # remove filtered cells from Seurat
    if(length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]  
    }

    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- 0

    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }

    valid_data <- data[, row.names(pd)]

    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)

    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL

        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS

      } else {
        # mist_list <- list(ident = ident) 
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    }

    if(1==1) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@assays$RNA@var.features)

    }
    monocle_cds@auxClusteringData$seurat <- mist_list

  } else if (class(otherCDS)[1] == 'SCESet') {
    requireNamespace("scater")

    message('Converting the exprs data in log scale back to original scale ...')    
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset

    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if("is.expr" %in% slotNames(otherCDS))
      lowerDetectionLimit <- otherCDS@is.expr
    else 
      lowerDetectionLimit <- 1

    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }

    if(import_all) {
      # mist_list <- list(iotherCDS@sc3,
      #                   otherCDS@reducedDimension)
      mist_list <- otherCDS 

    } else {
      mist_list <- list()
    }

    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    # monocle_cds@auxClusteringData$sc3 <- otherCDS@sc3
    # monocle_cds@auxOrderingData$scran <- mist_list

    monocle_cds@auxOrderingData$scran <- mist_list

  } else {
    stop('the object type you want to export to is not supported yet')
  }

  return(monocle_cds)
}
#############################################################################################################
monoclePlot <- function(seuratObj,clusRes,anchor=FALSE){
	######################################
	# do monocle analysis
	# 2018.3.7
	######################################
	if(!require(Cairo)){install.packages("Cairo")};library(Cairo)
	if(!require(stringr)){install.packages("stringr")};library(stringr)
	if(!require(ggplot2)){install.packages("ggplot2")};library(ggplot2)
	if(!require(Seurat)){install.packages("Seurat")};library(Seurat)
	if(!require(dplyr)){install.packages("dplyr")};library(dplyr)
	if(!require(monocle)){source("http://bioconductor.org/biocLite.R");biocLite("monocle")};library(monocle)
	if(!require(DDRTree)){source("http://bioconductor.org/biocLite.R");biocLite("DDRTree")};library(DDRTree)
	if(!require(HSMMSingleCell)){source("http://bioconductor.org/biocLite.R");biocLite("HSMMSingleCell")};library(HSMMSingleCell)
	## import data from Seurat obj
	gbm <- newimport(seuratObj, import_all = TRUE)
	### creat CellDataSet object
	my_cds <- newCellDataSet(exprs(gbm), 
    phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
    featureData = new("AnnotatedDataFrame", data = fData(gbm)),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size())
    #############################################
	### standardlization
	## Next, weâ€™ll perform some normalisation and variance estimation steps, 
	# which will be used in the differential expression analyses later on.
	my_cds <- estimateSizeFactors(my_cds)
	my_cds <- estimateDispersions(my_cds)
	# First we will need to run the detectGenes() function, 
	# which tallies the number of cells expressing a gene and the number 
	# of genes expressed among all cells. This will be clearer when you see the tables.
	# The num_cells_expressed column is a tally of the number of cells expressing a particular
	# gene (a gene is â€œexpressedâ€? if there is at least one count since we set min_expr = 0.1). 
	my_cds <- detectGenes(my_cds, min_expr = 0.1)
	my_cds_subset <- my_cds
    #############################################
	### plot data information
	#standardise to Z-distribution
	#x <- pData(my_cds)$nFeature_RNA
	#x_1 <- (x - mean(x)) / sd(x)
	#df <- data.frame(x = x_1)
	#p <- ggplot(df, aes(x)) + geom_histogram(bins = 50) + theme_gray()+ geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
	#CairoPNG(file="01.Norm_his.png",width=800,height=500)
	#plot(p)
	#dev.off()
	#p <- ggplot(pData(my_cds),aes(nCount_RNA,num_genes_expressed))+geom_point()
	#CairoPNG(file='02.GeneUMI_corPlot.png',width=800,height=800)
	#plot(p)
	#dev.off()
	#############################################
	## Constructing Single Cell Trajectories
	##-------------------1-Choose genes that define progress
	##----------PCA reduce dimensions
	#expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= 10))
	#my_cds_subset <- my_cds[expressed_genes, ]
	#my_cds_subset <- detectGenes(my_cds_subset, min_expr = 0.1)
	#fData(my_cds_subset)$use_for_ordering <- fData(my_cds_subset)$num_cells_expressed > 0.05 * ncol(my_cds_subset)
	# how many genes are used?
	#table(fData(my_cds_subset)$use_for_ordering)
	#FALSE  TRUE
	# 9141  6305
	#CairoPNG(file="03.subset_PseuElbow.png",width=800,height=500)
	#print(plot_pc_variance_explained(my_cds_subset, return_all = FALSE))
	#dev.off()
	##----------tsne reduce dimensions
	## tsne reduce dimension is used for clustering
	# my_cds_subset <- reduceDimension(my_cds_subset,
									 # max_components = 2,
									 # norm_method = 'log',
									 # scaling = TRUE,
									 # num_dim = 3,
									 # reduction_method = 'tSNE',
									 # verbose = TRUE)

	##-------------clustering cells 
	# HSMM_myo <- clusterCells(HSMM_myo, verbose = F)
	## After the clustering, we can check the clustering results.
	# plot_cell_clusters(HSMM_myo, color_by = 'as.factor(seurat_clusters)')
	# plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)')
	## We also provide the decision plot for users to check the Ρ, Δ for each cell and decide the threshold 
	## for defining the cell clusters.
	# plot_rho_delta(HSMM_myo, rho_threshold = 2, delta_threshold = 4 )
	# table(pData(my_cds_subset)$seurat_clusters)
	# #   1    2    3    4 
	# # 229  405  175 1251
	##-------------compare for diff exp gene of each cluster
	##After we confirm the clustering makes sense, 
	##we can then perform differential gene expression test as a way to extract the genes that distinguish them. 
	if(anchor==FALSE){
		clustering_DEG_genes <- differentialGeneTest(my_cds_subset,
													 fullModelFormulaStr = '~Cluster',
													 cores = 8)
		## differential gene expression table
		diff_exp_table <- clustering_DEG_genes %>% arrange(qval)
		write.csv(diff_exp_table,file="Pseu_diffGene.csv",row.names=F)
		## We'll use the top 1,000 most significantly differentially expressed genes as the set of ordering 
		## genes and perform the dimension reduction and the trajectory analysis (using the orderCells() function).
		## this should be change your marker genes
		## this is auto detected diff.Exp.Gene
		my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
	}else{
		##-------------- using the gene that seurat used in clustering
		##-------------- these gene is selected by top2000 from the FindVariableFeatures() by set nfeatures = 2000 or more
		##--------------this is not the same as single sample 
		my_ordering_genes <- anchor
	}

	my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = my_ordering_genes)
	##----------------plot for order genes
	CairoPNG(file="01.subset_dispGene.png",width=800,height=500)
	print(plot_ordering_genes(my_cds_subset))
	dev.off()

	#############################################
	# 2-Reduce the dimensionality of the data
	my_cds_subset <- reduceDimension(my_cds_subset,
									 max_components = 2,
									 norm_method = 'log',
									 scaling = TRUE,
									 num_dim = 3,
									 reduction_method = 'DDRTree',
									 verbose = TRUE)
	# my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree')

	#############################################
	##---------------------------3-Order cells in pseudotime
	####
	# the warnings were for use of deprecated code
	my_cds_subset <- orderCells(my_cds_subset)
	## plot
	CairoPNG(file="02.subset_ClusTrajectory.png",width=800,height=500)
	print(plot_cell_trajectory(my_cds_subset, color_by = "seurat_clusters"))
	dev.off()
	CairoPNG(file="02.subset_SplitClusTrajectory.png",width=1000,height=1000)
	print(plot_cell_trajectory(my_cds_subset, color_by = "orig.ident",cell_size=0.5) + facet_wrap(~seurat_clusters))
	dev.off()
	CairoPNG(file="03.subset_StateTrajectory.png",width=800,height=500)
	print(plot_cell_trajectory(my_cds_subset, color_by = "State"))
	dev.off()
	CairoPNG(file="03.subset_SplitStateTrajectory.png",width=1000,height=500)
	print(plot_cell_trajectory(my_cds_subset, color_by = "State") + facet_wrap(~State))
    dev.off()
    ### only for combined Data
        sampleNumer <- length(levels(as.factor(my_cds_subset$orig.ident)))
        if(sampleNumer > 1 ){
		CairoPNG(file="subset_SampleTrajectory.png",width=800,height=500)
		print(plot_cell_trajectory(my_cds_subset, color_by = "orig.ident"))
		dev.off()
		CairoPNG(file="subset_SplitSampleTrajectory.png",width=1000,height=500*round(sampleNumer/3))
		print(plot_cell_trajectory(my_cds_subset, cell_size = 0.5, color_by = "seurat_clusters")+ facet_wrap(~orig.ident))
		dev.off()
	}

	
	##################################################################
    #plot the cell dense distribution
    # my_cds_state <- function(cds){
    # if (length(unique(pData(cds)$State)) > 1){
    #     T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    #     return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    # }else {
    #     return (1)
    # }
    #}
    #my_cds_subset <- orderCells(my_cds_subset, root_state = my_cds_state(my_cds_subset))
    ##
    #CairoPNG(file="07.subset_DensTrajectory.png",width=1000,height=500)
    #print(plot_cell_trajectory(my_cds_subset, color_by = "Pseudotime"))
    #dev.off()	
	##########################################################
	## Finding Genes that Change as a Function of Pseudotime, and plot the top gene tracks
	## Once we have a trajectory, we can use differentialGeneTest() to find genes that have an expression 
	## pattern that varies according to pseudotime.
	my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 8)
	# get the top genes
	my_pseudotime_de %>% arrange(qval) %>% head(6) -> my_pseudotime_gene
	#write.csv(my_pseudotime_gene,file="pseudotime_gene.csv")
	my_pseudotime_gene <- my_pseudotime_gene$gene_short_name
	top_gene_num <- vector()
	for(gene in my_pseudotime_gene){
	    num <- grep(pattern=paste("^",gene,"$",sep=""),x=rownames(fData(my_cds_subset)))
	    top_gene_num <- c(top_gene_num,num)
	}
	# plot
	CairoPNG(file="04.TopGenePseu.png",width=1000,height=1000)
	print(plot_genes_in_pseudotime(my_cds_subset[top_gene_num,],color_by = "seurat_clusters",label_by_short_name = TRUE))
	dev.off()
	##########################################################
	### Clustering Genes by Pseudotemporal Expression Pattern
	# cluster the top 50 genes that vary as a function of pseudotime
	my_pseudotime_de %>% arrange(qval) %>% head(50) -> gene_to_cluster
	gene_to_cluster <- gene_to_cluster$gene_short_name
	top_gene_num <- vector()
	for(gene in gene_to_cluster){
	    num <- grep(pattern=paste("^",gene,"$",sep=""),x=rownames(fData(my_cds_subset)))
	    top_gene_num <- c(top_gene_num,num)
	} 
	#pdf('tsst.pdf')
	CairoPNG(file="05.Top50Pseu_Heatmap.png",width=800,height=1000)
	my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[top_gene_num,],
	                                                 num_clusters = 3,
	                                                 cores = 8,
	                                                 show_rownames = TRUE,
	                                                 return_heatmap = TRUE)
	dev.off()
	#CairoPDF(file='05.Top50Pseu_Heatmap.pdf',width=8.264,height=11.688)
	#print(my_pseudotime_cluster)
	#dev.off()
	write.csv(pData(my_cds_subset),file='Pseu_Info.Summary.csv',row.names=T)
	##########################################################
	return(list(my_cds_subset,my_cds))
	cat("Monocle work has done!!!","\n")
	cat("=============================","\n")
}
#############################################################################################################
#GO and pathway: analysis biological function for genes. 
#2019
#############################################################################################################
## shorten the names of GO_KEGG terms
shorten_names <- function(x, n_word=10, n_char=100){	
	if(!is.na(x)){
		if( length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > n_char) ){
			if(nchar(x) > n_char ) {
				x <- substr(x,1,n_char)
				x <- paste(paste(strsplit(x," ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)], collapse=" "), "...", sep="")
			} 
		}
	}
	return(x)
}
######################################################################
getGOTest <- function(geneFC,sampleName,organism,type){
	#############################################################
	# clusterProfile don GO and KEGG analysis
	if(!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager") 
	if(organism == 'human') {
		if(!require(org.Hs.eg.db)){ BiocManager::install('org.Hs.eg.db',ask = F,update = F) }; library(org.Hs.eg.db)
		organismAnn <- 'org.Hs.eg.db'
		EntrezAnnEnv <- org.Hs.eg.db
	}else if(organism == 'mouse') {
		if(!require(org.Mm.eg.db)){ BiocManager::install('org.Mm.eg.db',ask = F,update = F) }; library(org.Mm.eg.db)
		organismAnn <- 'org.Mm.eg.db'
		EntrezAnnEnv <- org.Mm.eg.db
	}else{
		cat('\n\n','!!! this organism is not found in annotation library or not to be setted !!!','\t',organism,'\n\n')
	}
	###########################################################################
	if(!require(clusterProfiler)){ BiocManager::install('clusterProfiler',ask = F,update = F) };library(clusterProfiler)
	if(!require(dplyr)){ install.packages("dplyr") };library(dplyr)
	###########################################################################
	# input gene from markerFind result
	entrezIDs.df <- bitr(geneFC$gene, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = EntrezAnnEnv)
	geneFC$SYMBOL <- geneFC$gene
	geneFC <- merge(geneFC,entrezIDs.df,by='SYMBOL',all.x=T)
	# entrezIDs <- mget(rownames(geneFC), EntrezAnnoEnv2, ifnotfound=NA)
	# geneFC$ENTREZID <- as.character(entrezIDs)
	geneFC.df <- na.omit(geneFC)
	display_number =vector()
	#GO enrichment with clusterProfiler analyse
	## keyType:keytype of input gene
	## ont： One of "MF", "BP", and "CC" subontologies.
	## pvalueCutoff:Cutoff value of pvalue.
	## pAdjustMethod: one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
	## universe: background genes
	## qvalueCutoff: qvalue cutoff
	## minGSSize: minimal size of genes annotated by Ontology term for testing.
	## maxGSSize: maximal size of genes annotated for testing
	## readable: whether mapping gene ID to gene Name(Symbol)
	## pool: If ont=’ALL’, whether pool 3 GO sub-ontologies
	## go ont: BP, CC, MF, ALL
	go_ont=c('BP','CC','MF')
	for(ont in go_ont) {
		ego <- enrichGO(gene = geneFC.df$ENTREZID,
		                 OrgDb=organismAnn,
		                 ont = ont,
		                 pAdjustMethod = "BH",
		                 minGSSize = 1,
		                 pvalueCutoff = 0.05,
		                 qvalueCutoff = 0.05,
		                 readable = TRUE)
		if(!is.null(ego)) {
			result <- na.omit(as.data.frame(ego@result))
			if (nrow(result) >=3){
				all_Annoterm <- vector()
				for(count in 1:nrow(result)){
				  all_Annoterm[count] <- as.numeric(strsplit(result$BgRatio[count],split = '/')[[1]][1])
				}
				result$RichFactor <- result$Count / all_Annoterm
				write.table(as.data.frame(result), file=paste(sampleName,type,ont,"GO.xls",sep="_"),sep="\t",row.names=FALSE)
				##############################################
				if(nrow(result)>20){ 
				  display_number=20 
				}else{ 
				  display_number=nrow(result) 
				}
				## rank
				result %>% arrange((p.adjust),desc(Count)) %>% head(display_number) -> ego_result
			    if(nrow(ego_result) >= 3){
				    ## GO bar figure
					for (count in 1:nrow(ego_result)){
					  ego_result$shortenTerm[count] <- shorten_names(as.character(ego_result$Description[count]))
					}
					## reorder the data to plot
				    require(ggplot2)
				    p <- ggplot(data=ego_result,mapping=aes(x=reorder(shortenTerm,-p.adjust),y=Count,fill=p.adjust)) + 
				    scale_fill_continuous(low='red',high ='blue') + xlab("GO Term") + ylab("Gene Number")+
				    geom_bar(stat="identity", width=0.8) + coord_flip() + 
				    labs(title=paste('Top',display_number,'of',ont,'Enrichment',sep=' ')) + theme_bw() +
				    theme(axis.text=element_text(face="bold",color="gray50"),plot.title=element_text(hjust=0.5))
				    ## plot and save
				    CairoPNG(file=paste(sampleName,type,ont,"GO.png",sep='_'),width=9,height=8.264,units='in',dpi=300)
				    print(p)
				    dev.off()

	                            if(!file.exists("pdf")) {dir.create('pdf')}
        	                    setwd("pdf")

				    CairoPDF(file=paste(sampleName,type,ont,"GO.pdf",sep='_'),width=9,height=9)
				    print(p)
				    dev.off()
				    
				    setwd('../')
				    ## GO dotplot enrichment analysis
				    pbubble = ggplot(ego_result,aes(x=reorder(shortenTerm,-p.adjust),y=RichFactor))+ 
				    geom_point() + coord_flip()+ geom_point(aes(size=Count,color=p.adjust)) + 
				    scale_colour_gradient(low="red",high="blue")
				    # plot GO enrichment dot graph
				    pgo = pbubble + labs(color=expression(p.adjust),size="Gene number",y="Rich Factor",title=paste('Top',display_number,'of',ont,'Enrichment',sep=' ')) +
				    theme_bw()+theme(plot.title = element_text(hjust = 0.5)) + xlab('GO Term')
				    ## plot and save
				    CairoPNG(file=paste(sampleName,type,ont,"GO_Enrichment.png",sep='_'),width=9,height=8.264,units='in',dpi=300)
				    print(pgo)
				    dev.off()

                                    if(!file.exists("pdf")) {dir.create('pdf')}
                                    setwd("pdf")

				    CairoPDF(file=paste(sampleName,type,ont,"GO_Enrichment.pdf",sep='_'),width=9,height=9)
				    print(pgo)
				    dev.off()  

				    setwd('../')
			    }else{
			    	cat('GO term is too short!!!','\n',sampleName,'\t',type,'\t',ont,'\n')
			    }
			}else{
				cat('there is no valuable GO term output!','\n',sampleName,'\t',type,'\t',ont,'\n')
			}
		}else{
			cat('there is no GO term output!','\n',sampleName,'\t',type,'\t',ont,'\n')
		}
	}
}
#############################################################################################################
##species: human(Hs), mouse(Mm), rat(Rn), chicken(Gg), pig(Ss), chimp(Pt), fly(Dm), worm(Ce), zebrafish(Dr)
## kegg pathview analysis
## just for all markers
getKEGGTest <- function(geneFC,sampleName,organism,type){
	#############################################################
	# clusterProfile don GO and KEGG analysis
	if(!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager") 
	if(organism == 'human') {
		if(!require(org.Hs.eg.db)){ BiocManager::install('org.Hs.eg.db',ask = F,update = F) }; library(org.Hs.eg.db)
		organismAnn <- 'org.Hs.eg.db'
		EntrezAnnEnv <- org.Hs.eg.db
	}else if(organism == 'mouse') {
		if(!require(org.Mm.eg.db)){ BiocManager::install('org.Mm.eg.db',ask = F,update = F) }; library(org.Mm.eg.db)
		organismAnn <- 'org.Mm.eg.db'
		EntrezAnnEnv <- org.Mm.eg.db
	}else{
		cat('\n\n','!!! this organism is not found in annotation library or not to be setted !!!','\t',organism,'\n\n')
	}

	###########################################################################
	if(!require(clusterProfiler)){ BiocManager::install('clusterProfiler',ask = F,update = F) };library(clusterProfiler)
	if(!require(dplyr)){ install.packages("dplyr") };library(dplyr)
	###########################################################################
	# input gene from markerFind result
	# entrezIDs.df <- bitr(rownames(geneFC), fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
	# geneFC$ENTREZID <- entrezIDs.df$ENTREZID
	# entrezIDs <- mget(rownames(geneFC), EntrezAnnoEnv2, ifnotfound=NA)
	# geneFC$ENTREZID <- as.character(entrezIDs)
	entrezIDs.df <- bitr(geneFC$gene, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = EntrezAnnEnv)
	geneFC$SYMBOL <- geneFC$gene
	geneFC <- merge(geneFC,entrezIDs.df,by='SYMBOL',all.x=T)
	# entrezIDs <- mget(rownames(geneFC), EntrezAnnoEnv2, ifnotfound=NA)
	# geneFC$ENTREZID <- as.character(entrezIDs)
	geneFC.df <- na.omit(geneFC)
	#############################################################
	## KEGG analysis
	kegg <- enrichKEGG(gene = geneFC.df$ENTREZID,
		organism =organism,
		pvalueCutoff = 0.05,
		qvalueCutoff = 0.05,
		minGSSize = 1,
		#readable = TRUE,
		use_internal_data =FALSE)
	kegg <- setReadable(kegg, EntrezAnnEnv, keyType = "ENTREZID")
	if(!is.null(kegg)){
		result <- na.omit(as.data.frame(kegg@result))
		if(nrow(result) >= 3){
			all_Annoterm <- vector()
			for(count in 1:nrow(result)){
				all_Annoterm[count] <- as.numeric(strsplit(result$BgRatio[count],split = '/')[[1]][1])
			}
			result$RichFactor <- result$Count / all_Annoterm
			write.table(as.data.frame(result), file=paste(sampleName,type,"Pathway.xls",sep='_'),sep="\t",row.names=FALSE)
			##########################################
			## pathway富集分析bar
			if(nrow(result) >= 20 ){
				display_number=20
			}else{
				display_number <- nrow(result)
			}
			## rank
			na.omit(result) %>% arrange((p.adjust),desc(Count)) %>% head(display_number)-> kegg_result
			if(nrow(kegg_result) >= 3){
			    for (count in 1:nrow(kegg_result)){
			      kegg_result$shortenTerm[count] <- shorten_names(as.character(kegg_result$Description[count]))
			    }
			    require(ggplot2)
			    p <- ggplot(data=kegg_result, aes(x=reorder(shortenTerm,-p.adjust), y=Count, fill=p.adjust)) + 
			      scale_fill_continuous(low = 'red', high = 'blue') + geom_bar(stat="identity", width=0.8) + 
			      theme_bw() + coord_flip() + xlab("Pathway Term") + ylab("Gene Number") +
			      labs(title=paste('Top',display_number,'of Pathway Enrichment',sep=' ')) + 
			      theme(axis.text=element_text(face = "bold", color="gray50"),plot.title = element_text(hjust = 0.5))
			    #p
			    CairoPNG(file=paste(sampleName,type,"Pathway.png",sep='_'),width=9,height=8.264,units='in',dpi=300)
			    plot(p)
			    dev.off()

                            if(!file.exists("pdf")) {dir.create('pdf')}
                            setwd("pdf")

			    CairoPDF(file=paste(sampleName,type,"Pathway.pdf",sep='_'),width=8,height=8.5)
			    plot(p)
			    dev.off()

                            setwd('../')

			    #############################################
			    ## Pathway富集分析气泡?
			    pbubble = ggplot(data=kegg_result,aes(x=reorder(shortenTerm,-p.adjust),y=RichFactor)) + geom_point() +
			    theme_bw() + coord_flip() + geom_point(aes(size=Count,color=p.adjust)) + 
			    scale_colour_gradient(low="red",high="blue") 
			    pkegg = pbubble + labs(color=expression(p.adjust),
			    size="Gene number", x="Pathway Term",y="Rich Factor",title=paste('Top',display_number,'of Pathway Enrichment',sep=' ')) + 
			    theme(plot.title = element_text(hjust = 0.5)) #
			    
			    #CairoPNG(file=paste(sampleName,type,"Pathway_Enrichment.png",sep='_'),width =par('din')[1], height = par('din')[2], units='in', dpi=600)
			    CairoPNG(file=paste(sampleName,type,"Pathway_Enrichment.png",sep='_'),width=9,height=8.264,units='in',dpi=300)
			    plot(pkegg)
			    dev.off()

			    if(!file.exists("pdf")) {dir.create('pdf')}	
		            setwd("pdf")

			    #CairoPDF(file=paste(sampleName,type,"Pathway_Enrichment.png",sep='_'),width =par('din')[1], height = par('din')[2], units='in', dpi=600)
			    CairoPDF(file=paste(sampleName,type,"Pathway_Enrichment.pdf",sep='_'),width=9,height=8.264)
			    plot(pkegg)
			    dev.off()

			    setwd('../')

			}else{
				cat('kegg term is too short!!!','\n',sampleName,'\t',type,'\n')
			}
		}else{
			cat('there is no valuable kegg term output!','\n',sampleName,'\t',type,'\n')
		}
	}else{
		cat('there is no kegg term output!','\n',sampleName,'\t',type,'\n')
	}
}
