#############################################################################################################
#seurat functions including: getInitialPNG, getFiltObj, getOCA_TSNE, getGraphBased, markerPlot, markerFind
#by jiangxc
#2019.3.12
#############################################################################################################
getInitialPNG <- function(seuratObj){
	## view and evaluation the initial data
	mito.genes1 <- grep(pattern = "^mt-",x = rownames(seuratObj@data),ignore.case = TRUE ,value = TRUE)
	mito.genes2 <- grep(pattern = "^mrp",x = rownames(seuratObj@data),ignore.case = TRUE ,value = TRUE)
	mito.genes <- c(mito.genes1,mito.genes2)
	percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ])/Matrix::colSums(seuratObj@raw.data)
	#add the percent.mito to metagene matrix
	seuratObj <- AddMetaData(object=seuratObj,metadata=percent.mito,col.name="percent.mito")
	#VlnPlot the nGene,nUMI and percent.mito
	CairoPNG(file="01.vinPlot.png",width=800,height=500)
	print(VlnPlot(object = seuratObj, features.plot = c("nUMI","nGene","percent.mito"), nCol = 3))
	dev.off()
	CairoPNG(file="02.corPlot.png",width=900,height=500)
	par(mfrow = c(1, 2))
	GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "percent.mito")
	GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "nGene")
	dev.off()
	## 
	return(seuratObj)
	cat("getInitialPNG work has done!!!","\n")
	cat("=============================","\n")
}
#############################################################################################################
getFiltObj <- function(seuratObj, LowThres.nUMI=0, HighThres.nUMI=Inf, LowThres.nGene=200, HighThres.nGene=Inf, LowThres.MitoPercent=0,HighThres.MitoPercent=0.1){
	path_step1 <- getwd()
	############################################################
	## remove low qulality cells
	############################################################
	## filter data by UMI,Gene and percent.mito 
	# seuratObj <- FilterCells(object = seuratObj, subset.names = c("nUMI"), low.thresholds =0 , high.thresholds = umiThres)
	# # filter by Gene number
	# seuratObj <- FilterCells(object = seuratObj, subset.names = c("nGene"), low.thresholds =200 , high.thresholds = geneThres)
	# # filter by mito.percent
	# seuratObj <- FilterCells(object = seuratObj, subset.names = c("percent.mito"), low.thresholds =-Inf , high.thresholds = 0.05)
	seuratObj <- FilterCells(object = seuratObj, subset.names = c('nUMI','nGene','percent.mito'), 
		low.thresholds = c(LowThres.nUMI, LowThres.nGene, LowThres.MitoPercent) , 
		high.thresholds = c(HighThres.nUMI, HighThres.nGene, HighThres.MitoPercent))

	CairoPNG(file="03.filt_vinPlot.png",width=800,height=500)
	print(VlnPlot(object = seuratObj, features.plot = c( "nUMI","nGene", "percent.mito"), nCol = 3))
	dev.off()
	CairoPNG(file="04.filt_corPlot.png",width=900,height=500)
	par(mfrow = c(1, 2))
	GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "percent.mito")
	GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "nGene")
	dev.off()
	############################################################
	## remove ribosome and mito genes
	############################################################
	setwd('../')
	filter.mito <- grep(pattern = "^mt-" , x = rownames(seuratObj@data), ignore.case = TRUE, value = TRUE, invert = TRUE)
	filter.mito <- grep(pattern = "^mrp" , x = filter.mito, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filter.mito.rp <- grep(pattern = "^rpl" , x = filter.mito, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filter.mito.rp <- grep(pattern = "^rps" , x = filter.mito.rp, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filter.mito.rp <- grep(pattern = "^rpf" , x = filter.mito.rp, ignore.case = TRUE, value = TRUE, invert = TRUE)
	filter.mito.rp <- grep(pattern = "^rpn" , x = filter.mito.rp, ignore.case = TRUE, value = TRUE, invert = TRUE)
	write.csv(filter.mito.rp, file = "All.genes.use.csv", row.names = FALSE, col.names = FALSE)
	## save the filtered matrix
	genes.used <- filter.mito.rp
	data.use <- t(data.frame(FetchData(seuratObj, vars.all = genes.used)))
	write.table(data.use, file = "Filtered_bc_matrix.xls", sep = "\t")
	#data.use <- read.table(file = "Filtered_bc_matrix.xls", header=T, row.names=1, sep = "\t")
	# using data that filtered low-quality barcodes and mito-ribosome genes 
	###########################################################################################
	setwd(path_step1)
	seuratObj <- CreateSeuratObject(raw.data = data.use, project =  seuratObj@project.name)
	CairoPNG(file="05.final_vinPlot.png",width=550,height=500)
	print(VlnPlot(object = seuratObj, features.plot = c( "nUMI","nGene"), nCol = 2))
	dev.off()
	CairoPNG(file="06.final_corPlot.png",width=600,height=500)
	GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "nGene")
	dev.off()
	############################################################################################
	return(seuratObj)
	cat("getFiltObj1 work has done!!!","\n")
	cat("=============================","\n")
}
#############################################################################################################
getStandardization <- function(seuratObj){
	# data standardization
	seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
	#FindVariableGenes
	CairoPNG(file="07.VarGenes_disp.png",width=800,height=600)
	seuratObj <- FindVariableGenes(object = seuratObj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
	dev.off()
	# center, scale and regress
	seuratObj <- ScaleData(object = seuratObj, do.scale = TRUE, do.center = TRUE)
	return(seuratObj)
	cat("getStandardization work have done!!!","\n")
	cat("====================================","\n")
}
#############################################################################################################
getPCA_TSNE <- function(seuratObj,PCA_dims=20){
	# PCA reduce the dimension
	# do PCA
	seuratObj <- RunPCA(object = seuratObj, do.print = TRUE, pcs.print = 1:3, genes.print = 3)
	# Visualization: distribution
	CairoPNG(file="08.PCA.png",width=700,height=600)
	print(PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 2))
	dev.off()
	#PCs 
	CairoPNG(file="09.PC1-2.png",width=700,height=600)
	print(VizPCA(object = seuratObj, pcs.use = 1:2))
	dev.off()
	# PCs elbow line
	CairoPNG(file="10.Elbow.png",width=600,height=400)
	print(PCElbowPlot(object = seuratObj))
	dev.off()
	#
	# t-SNE plot
	seuratObj <- RunTSNE(object = seuratObj, dims.use = 1:PCA_dims, do.fast = TRUE, reduction.use = "pca")
	CairoPNG(file="11.tSNE.png",width=700,height=600)
	print(TSNEPlot(object = seuratObj))
	dev.off()
	return(seuratObj)
	cat("getPCA_TSNE work have done!!!","\n")
	cat("==============================","\n")
}
#############################################################################################################
getGraphBased <- function(seuratObj,PCA_dims=20,resolution=0.6){
	# 3k cells: res=(0.6,1.2) is best
	# >3k cells: res-should to test
	# PCs= 10-20
	# set the parameters 
	# w3 pcs=20, w7 pcs=8
	seuratObj <- FindClusters(object = seuratObj, reduction.type = "pca", dims.use = 1:PCA_dims, resolution = resolution, print.output = 0, save.SNN = TRUE)
	CairoPNG(file="GraphBased_TSNE.png",width=700,height=600)
	print(TSNEPlot(object = seuratObj))
	dev.off()
	CairoPNG(file="GraphBased_TSNELabel.png", width = 700, height = 600)
	print(TSNEPlot(object = seuratObj, do.label = TRUE,label.size = 5))
	dev.off()
	##########################################################
	# export metaData
	metaData<- data.frame(barcode=rownames(seuratObj@meta.data),nGene=seuratObj@meta.data$nGene,nUMI=seuratObj@meta.data$nUMI,sample=seuratObj@meta.data$orig.ident,cluster=seuratObj@meta.data$res)
	metaData$barcode <- str_sub(metaData$barcode,start=-16)
	write.csv(metaData,file='meta_data.csv',row.names=F)
	# count the cell number in subCluster and sample
	count_cell <- table(metaData[,c(4,5)])
	write.csv(t(count_cell),file='count_cell_num.csv',row.names=T)
	count_cell <- read.csv(file='count_cell_num.csv',header=T)
	colnames(count_cell)[1]='clusterID'
	# count_cell$sum <- rowSums(count_cell[,2:ncol(count_cell)])
	write.csv(count_cell,file='count_cell_num.csv',row.names=F)
	####################################################################
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
	for (counterCol in colN){
		dir.create(marker_type[counterCol])
		path.loop <- paste(path_save,marker_type[counterCol],sep="/")
		setwd(path.loop)
		rowN <- 1:length(marker_list[,counterCol])
		for(counterRow in rowN){
			gene_name <- as.character(marker_list[counterRow,counterCol])
			if(nchar(gene_name)>0){
				if( length( grep(pattern=paste("^",gene_name,"$",sep=""),x=rownames(seuratObj@data),ignore.case = TRUE) ) ){
					gene_name <- grep(pattern=paste("^",gene_name,"$",sep=""),x=rownames(seuratObj@data),ignore.case = TRUE, value=TRUE)
					name.vinPlot <- paste(gene_name,"vin.png",sep="_")
					name.tsnePlot <- paste(gene_name,"tsne.png",sep="_")
					CairoPNG(file=name.vinPlot,width=700,height=600)
					print(VlnPlot(object = seuratObj, features.plot = gene_name, use.raw = TRUE, y.log = TRUE))
					dev.off()
					CairoPNG(file=name.tsnePlot,width=700,height=600)
					print(FeaturePlot(object = seuratObj, features.plot = gene_name, cols.use = c("lightgrey", "blue"), reduction.use = "tsne"))
					dev.off()
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
	gene_list <- unique(as.character(gene_list))
	CairoPNG(file = "DotPlot_allGene.png",width=200+50*length(gene_list),height=300+100*length(levels(seuratObj@ident)))
	print(DotPlot(seuratObj, genes.plot = rev(gene_list), cols.use = c("blue","red"), dot.scale = 8, plot.legend = T, x.lab.rot = T, do.return = T))
	dev.off()
	if(length(levels(as.factor(seuratObj@meta.data$stim)))>1) {
		CairoPNG(file = "splitDotPlot_allGene.png",width=200+50*length(gene_list),height=500+150*(length(levels(seuratObj@ident))+length(levels(as.factor(seuratObj@meta.data[,4])))))
		print(SplitDotPlotGG(seuratObj, genes.plot = rev(gene_list), cols.use = c("blue","red","green","purple","coral","gold2","deepskyblue","darkorange"), x.lab.rot = T, plot.legend = T, dot.scale = 8, do.return = F, grouping.var = "stim"))
		dev.off()
	}

	cat("markerPlot work has done!!!","\n")
	cat("=============================","\n")
}
#############################################################################################################
markerFind <- function(seuratObj,test_method="MAST",min.pct=0.25,logFC=0.25){
	if (min.pct>=0 && min.pct<=1 ){
		ClusterNum <- levels(seuratObj@ident)
		for (numC in ClusterNum ){
			name.cluster = paste("cluster",numC,test_method,"markers.csv",sep="_")
		  	cluster.markers <- FindMarkers(object = seuratObj, ident.1 = numC, test.use = test_method, min.pct = min.pct, logfc.threshold=logFC)
		  	write.csv(cluster.markers, file = name.cluster)
		}
		# find markers for every cluster compared to all remaining cells
		# only the positive ones
		seurat.markers <- FindAllMarkers(object = seuratObj, only.pos = TRUE, test.use = test_method, min.pct = min.pct, logfc.threshold=logFC)
		# write.csv(seurat.markers, file = "AllMarkers.csv")
		# plot top gene heatmap
		top20 <- seurat.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
		CairoPNG(file="Top20Heatmap.png",width=4000,height=3800)
		print(DoHeatmap(object = seuratObj, genes.use = top20$gene, slim.col.label = TRUE, remove.key = TRUE, cex.row= 20, group.cex= 50))
		dev.off()
	}
	else{
		print("Error!")
		print("Please check input data. The input data may out of range or with incorrect format!")
	}
	cat("markerFind work has done!!!","\n")
	cat("=============================","\n")
}
#############################################################################################################
#monocle:analysis pesudotime using seurat result. 
#by jiangxc
#2019.3.12
#############################################################################################################
monoclePlot <- function(seuratObj,clusRes){
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
	gbm <- importCDS(seuratObj, import_all = TRUE)
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
    #############################################
	### plot data information
	#standardise to Z-distribution
	x <- pData(my_cds)$nGene
	x_1 <- (x - mean(x)) / sd(x)
	df <- data.frame(x = x_1)
	p <- ggplot(df, aes(x)) + geom_histogram(bins = 50) + theme_gray()+ geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
	CairoPNG(file="01.Norm_his.png",width=800,height=500)
	plot(p)
	dev.off()
	p <- ggplot(pData(my_cds),aes(nUMI,num_genes_expressed))+geom_point()
	CairoPNG(file='02.GeneUMI_corPlot.png',width=800,height=800)
	plot(p)
	dev.off()
	#############################################
	## Constructing Single Cell Trajectories
	## 1-Choose genes that define progress
	expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= 10))
	my_cds_subset <- my_cds[expressed_genes, ]
	my_cds_subset <- detectGenes(my_cds_subset, min_expr = 0.1)
	fData(my_cds_subset)$use_for_ordering <- fData(my_cds_subset)$num_cells_expressed > 0.05 * ncol(my_cds_subset)
	# how many genes are used?
	table(fData(my_cds_subset)$use_for_ordering)
	#FALSE  TRUE
	# 9141  6305
	CairoPNG(file="03.subset_PseuElbow.png",width=800,height=500)
	print(plot_pc_variance_explained(my_cds_subset, return_all = FALSE))
	dev.off()
	#############################################
	# 2-Reduce the dimensionality of the data
	my_cds_subset <- reduceDimension(my_cds_subset,
	                                 max_components = 2,
	                                 norm_method = 'log',
	                                 scaling = TRUE,
	                                 num_dim = 30,
	                                 reduction_method = 'DDRTree',
	                                 verbose = TRUE)
	getClusSite <- grep(pattern=paste('res',clusRes,sep='.'),x=colnames(pData(my_cds_subset)))
	colnames(pData(my_cds_subset))[getClusSite]<- 'Cluster'
	table(pData(my_cds_subset)$Cluster)
	#   1    2    3    4 
	# 229  405  175 1251
	## Now weâ€™ll perform the differential gene expression analysis as before but across all cell clusters.
	##After we confirm the clustering makes sense, 
	##we can then perform differential gene expression test as a way to extract the genes that distinguish them. 
	clustering_DEG_genes <- differentialGeneTest(my_cds_subset,
	                                             fullModelFormulaStr = '~Cluster',
	                                             cores = 8)
	## differential gene expression table
	diff_exp_table <- clustering_DEG_genes %>% arrange(qval)
	write.csv(diff_exp_table,file="Pseu_diffGene.csv",row.names=F)
	#############################################
	# 3-Order cells in pseudotime
	## Weâ€™ll use the top 1,000 most significantly differentially expressed genes as the set of ordering 
	## genes and perform the dimension reduction and the trajectory analysis (using the orderCells() function).
	my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
	my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = my_ordering_genes)
	####
	## plot 
	CairoPNG(file="04.subset_dispGene.png",width=800,height=500)
	print(plot_ordering_genes(my_cds_subset))
	dev.off()
	####
	my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree')
	# the warnings were for use of deprecated code
	my_cds_subset <- orderCells(my_cds_subset)
	## plot
	CairoPNG(file="05.subset_ClusTrajectory.png",width=800,height=500)
	print(plot_cell_trajectory(my_cds_subset, color_by = "Cluster"))
	dev.off()
	CairoPNG(file="06.subset_SplitClusTrajectory.png",width=1000,height=500)
	print(plot_cell_trajectory(my_cds_subset, color_by = "Cluster") + facet_wrap(~State))
    dev.off()
	CairoPNG(file="07.subset_StateTrajectory.png",width=800,height=500)
	print(plot_cell_trajectory(my_cds_subset, color_by = "State"))
	dev.off()
	CairoPNG(file="08.subset_SplitStateTrajectory.png",width=1000,height=500)
	print(plot_cell_trajectory(my_cds_subset, color_by = "State") + facet_wrap(~State))
    dev.off()
    ### only for combined Data
    if(length(levels(as.factor(my_cds_subset$orig.ident))) > 1 ){
		CairoPNG(file="subset_SampleTrajectory.png",width=800,height=500)
		print(plot_cell_trajectory(my_cds_subset, color_by = "orig.ident"))
		dev.off()
		CairoPNG(file="subset_SplitSampleTrajectory.png",width=800,height=500)
		print(plot_cell_trajectory(my_cds_subset, color_by = "orig.ident")+ facet_wrap(~State))
		dev.off()
	}
    ##################################################################
    #plot the cell dense distribution
 #    my_cds_state <- function(cds){
	#     if (length(unique(pData(cds)$State)) > 1){
	#         T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
	#         return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
	#     }else {
	#         return (1)
	#     }
 #    }
 #    my_cds_subset <- orderCells(my_cds_subset, root_state = my_cds_state(my_cds_subset))
 #    ##
 # CairoPNG(file="09.subset_DensTrajectory.png",width=1000,height=500)
 #    print(plot_cell_trajectory(my_cds_subset, color_by = "Pseudotime"))
 #    dev.off()	
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
	CairoPNG(file="09.TopGenePseu.png",width=1000,height=1000)
	print(plot_genes_in_pseudotime(my_cds_subset[top_gene_num,],color_by = "Cluster",label_by_short_name = TRUE))
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
	CairoPNG(file="10.Top50Pseu_Heatmap.png",width=800,height=1000)
	my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[top_gene_num,],
	                                                 num_clusters = 3,
	                                                 cores = 8,
	                                                 show_rownames = TRUE,
	                                                 return_heatmap = TRUE)
	dev.off()
	CairoPDF(file='10.Top50Pseu_Heatmap.pdf',width=8.264,height=11.688)
	print(my_pseudotime_cluster)
	dev.off()
	write.csv(pData(my_cds_subset),file='Pseu_Info.Summary.csv',row.names=T,append=F)
	##########################################################
	return(list(my_cds_subset,my_cds))
	cat("Monocle work has done!!!","\n")
	cat("=============================","\n")
}
#############################################################################################################
#GO and pathway: analysis biological function for genes. 
#by jiangxc
#2019.3.18
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
	entrezIDs.df <- bitr(rownames(geneFC), fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = EntrezAnnEnv)
	geneFC$SYMBOL <- rownames(geneFC)
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
			if (nrow(result) > 5){
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
			    if(nrow(ego_result) >= 5){
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
				    CairoPDF(file=paste(sampleName,type,ont,"GO.pdf",sep='_'),width=9,height=9)
				    print(p)
				    dev.off()

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
				    CairoPDF(file=paste(sampleName,type,ont,"GO_Enrichment.pdf",sep='_'),width=9,height=9)
				    print(pgo)
				    dev.off() 				    
				 #    ## plot GO enrichment DAG graph
				 #    if((type == 'all' || type == 'All' || type == 'ALL') && nrow(ego_result) >= 10){
				 #    	## top 10 node		    
					#     CairoPDF(file=paste(sampleName,type,ont,"GO_DAG.png",sep='_'),width=30,height=25)
					#     goplot(ego,showCategory = 10)
					#     dev.off()
					# }

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
	entrezIDs.df <- bitr(rownames(geneFC), fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = EntrezAnnEnv)
	geneFC$SYMBOL <- rownames(geneFC)
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

	if(!is.null(kegg)){
		result <- na.omit(as.data.frame(kegg@result))
		if(nrow(result) > 5){
			all_Annoterm <- vector()
			for(count in 1:nrow(result)){
				all_Annoterm[count] <- as.numeric(strsplit(result$BgRatio[count],split = '/')[[1]][1])
			}
			result$RichFactor <- result$Count / all_Annoterm
			write.table(as.data.frame(result), file=paste(sampleName,type,"KEGG.xls",sep='_'),sep="\t",row.names=FALSE)
			##########################################
			## pathway富集分析bar
			if(nrow(result) >= 20 ){
				display_number=20
			}else{
				display_number <- nrow(result)
			}
			## rank
			na.omit(result) %>% arrange((p.adjust),desc(Count)) %>% head(display_number)-> kegg_result
			if(nrow(kegg_result) >= 5){
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
			    CairoPDF(file=paste(sampleName,type,"Pathway.pdf",sep='_'),width=8,height=8.5)
			    plot(p)
			    dev.off()
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
			    #CairoPDF(file=paste(sampleName,type,"Pathway_Enrichment.png",sep='_'),width =par('din')[1], height = par('din')[2], units='in', dpi=600)
			    CairoPDF(file=paste(sampleName,type,"Pathway_Enrichment.pdf",sep='_'),width=9,height=8.264)
			    plot(pkegg)
			    dev.off()
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
######################################################################
#pathview plot
#input data: path_id and gene_id and gene_exp
#
getPathview <- function(geneFC,sampleName,organism,type){
	if(type == 'all' | type == 'All' | type == 'ALL'){
		gEns <- unlist(geneFC.df$ENTREZID)
		gene.data <- geneFC.df$avg_logFC
		names(gene.data) <- gEns
		## remove no annotation kegg ID
		kegg.df <- kegg@result[kegg@result$pvalue < 0.05,]
		kegg.df <- kegg.df[-grep(pattern='04215$',x=kegg.df$ID),]
		kegg.df <- kegg.df[-grep(pattern='04723$',x=kegg.df$ID),]
		if(nrow(kegg.df) >= 5){
		    ## creat folder
		    path_view =paste(sampleName,'Pathview',type,sep='_')
		    if(!file.exists(path_view)) {dir.create(path_view)}
		    setwd(path_view)
		    require(pathview)
		    ## pathview colour bar set threshold(deafult): limit = list(gene=c(-1,1), cpd=c(-1,1))
		    ## set colour bar align style: keys.align = "y"
		    ## set colour bar position: key.pos = "topright"
		    ## set whether the expression data and gene show in one graph: match.data = F
		    ## multiple sample show in one graph: multi.state = T
		    ## pathway color(deafult): low = list(gene = "green"), mid =list(gene = "gray"), high = list(gene = "red"),
		    for(i in 1:nrow(kegg.df)){ 
		      pv.out <- pathview::pathview(gene.data, pathway.id=as.character(kegg.df$ID)[i], species=organism, limit = list(gene=c(floor(min(gene.data)),ceiling(max(gene.data)))),
		                        kegg.native=TRUE,min.nnodes=3,Map.nul=FALSE,expand.node=FALSE, map.symbol=TRUE,gene.annotpkg=NULL,
		                        split.group=FALSE,new.signature=TRUE,plot.col.key=TRUE,key.pos='topright',keys.align='y')
		    }
		    setwd('../')
		}else{
			cat('pathview cannot be plot, because of too short kegg ID!!!','\n',sampleName,'\t',type,'\n')
		}
	}
}
#############################################################################################################
#topGO for GO enrichment DAG graph
#input data: geneList 
getTopGO <- function(geneFC,sampleName,organism,type){
	if(type == 'all' | type == 'All' | type == 'ALL'){
	    require(topGO)
	    ## creat folder
	    GOenrich_view =paste(sampleName,'DAG',type,sep='_')
	    if(!file.exists(GOenrich_view)) {dir.create(GOenrich_view)}
	    setwd(GOenrich_view)
	    ## for no annotation package organism using this code
	    # EG2GO <- toTable(org.Hs.egGO)
	    # geneID2Go <- by(EG2GO$go_id,EG2GO$ENTREZID,function(x) as.character(x))	    
	    interesting_genes<-factor(geneFC.df[abs(geneFC.df$avg_logFC) >1,]$ENTREZID)
	    all_genes <- geneFC.df$ENTREZID
	    ## adjust for gene number, 
	    if(length(all_genes) > 3 && length(interesting_genes)>0 && length(all_genes)>length(interesting_genes)){
		    geneList <- factor(as.integer(all_genes %in% interesting_genes))
		    names(geneList)=all_genes    
		    go_ontology <- c('BP','CC','MF')
		    for (ont in go_ontology){
				sampleGOdata <- new("topGOdata", 
				                  ontology = ont,
				                  allGenes = geneList, 
				                  nodeSize = 10,
				                  annot = annFUN.org,
				                  mapping = organismAnn,
				                  ID = 'entrez')
				if(!is.null(sampleGOdata)){
					resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher") 
					#resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks") 
					#resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks") 
					if(length(sampleGOdata@graph@nodes) >= 5){
						nodes = length(sampleGOdata@graph@nodes)
						if(length(sampleGOdata@graph@nodes) >= 100){ nodes = 100 }
						allRes <- GenTable(sampleGOdata, 
						#classicKS = resultKS, 
						#elimKS = resultKS.elim, 
						#orderBy = "classicFisher", 
						#ranksOf = "classicFisher", #,"elim","elimKS")
						classicFisher = resultFisher,
						topNodes = nodes)
						write.table(allRes, file=paste(sampleName,type,ont,'DAG.xls',sep='_'),row.names=F,sep="\t")
						CairoPNG(file=paste(sampleName,type,ont,'DAG.png',sep='_'),width=800,height=1000)
						showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes=5, useInfo='all')
						dev.off()
						CairoPDF( file=paste(sampleName,type,ont,'DAG.pdf',sep='_'),width=8.264,height=11.688) 
						showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes=5, useInfo='all') 
						dev.off()
					}else{
						cat('the topGO-DAG node is too short, no result output!!!','\n',sampleName,'\t',type,'\t',ont,'\n')
					}
				}else{
					cat('there is no sampleGOdata to be found!!!','\n',sampleName,'\t',type,'\t',ont,'\n')
				}
		    }
	    }else{
	    	cat('the gene or interesting_genes number is too short, or they are equal, there is no DAG result output!!!',
	    		'\n',sampleName,'\t',type,'\n')
	    }
	    setwd('../')
	}
}


############################################################
# for sample combination
CCAqc <- function(seuratObj,PCs){
	setwd(path.step1)
	p1 <- DimPlot(object = seuratObj, reduction.use = "cca", group.by = "stim", pt.size = 0.5, do.return = TRUE)
	CairoPNG(file="01.CCA_dim.png",width=800,height=700)
	print(plot_grid(p1))
	dev.off()
	p2 <- VlnPlot(object = seuratObj, features.plot = "CC1", group.by = "stim", do.return = TRUE)
	CairoPNG(file="02.CCA_vin.png",width=800,height=700)
	print(plot_grid(p2))
	dev.off()

	CairoPNG(file="03.CCA_Meta.png",width=800,height=500)
	p3 <- MetageneBicorPlot(seuratObj, grouping.var = "stim", dims.eval = 1:PCs, display.progress = FALSE)
	dev.off()
	CairoPNG(file="04.CCA_DimHeat.png",width=1500,height=1500) 
	print(DimHeatmap(object = seuratObj, reduction.type = "cca", cells.use = 500, dim.use = 1:PCs, do.balanced = TRUE))
	dev.off()

	seuratObj <- AlignSubspace(seuratObj, reduction.type = "cca", grouping.var = "stim", dims.align = 1:PCs)

	p1 <- VlnPlot(object = seuratObj, features.plot = "ACC1", group.by = "stim", 
	    do.return = TRUE)
	p2 <- VlnPlot(object = seuratObj, features.plot = "ACC2", group.by = "stim", 
	    do.return = TRUE)

	CairoPNG(file="05.CCA_DimVin.png",width=1000,height=500)
	print(plot_grid(p1, p2))
	dev.off()
	return(seuratObj)

}

CCAcluster <- function(seuratObj,PCs,clusRes){
	seuratObj <- RunTSNE(seuratObj, reduction.use = "cca.aligned", dims.use = 1:PCs, do.fast = T)
	seuratObj <- FindClusters(seuratObj, reduction.type = "cca.aligned", resolution = clusRes, dims.use = 1:PCs)
	# Visualization
	CairoPDF(file="Cluster_TSNE.pdf")
	p1 <- TSNEPlot(seuratObj, do.return = F, pt.size = 0.5, group.by = "stim")
	p2 <- TSNEPlot(seuratObj, do.label = T, do.return = F, pt.size = 0.5)
	dev.off()
	CairoPNG(file="Cluster_TSNE.png",width=1000,height=500)
	print(plot_grid(p1, p2))
	dev.off()
	# export the metaData
	resCol <- grep(pattern=paste('res',clusRes,sep='.'),x=colnames(seuratObj@meta.data))
	metaData<- data.frame(barcode=rownames(seuratObj@meta.data),nGene=seuratObj@meta.data$nGene,
		nUMI=seuratObj@meta.data$nUMI,sample=seuratObj@meta.data$orig.ident,cluster=seuratObj@meta.data[,resCol])	
	metaData$barcode <- str_sub(metaData$barcode,start=-16)
	write.csv(metaData,file='meta_data.csv',row.names=F)
	# count the cell number in subCluster and sample
	count_cell <- table(metaData[,c(4,5)])
	write.csv(t(count_cell),file='count_cell_num.csv',row.names=T)
	count_cell <- read.csv(file='count_cell_num.csv',header=T)
	colnames(count_cell)[1]='clusterID'
	count_cell$sum <- rowSums(count_cell[,2:ncol(count_cell)])
	write.csv(count_cell,file='count_cell_num.csv',row.names=F)

	return(seuratObj)

}




