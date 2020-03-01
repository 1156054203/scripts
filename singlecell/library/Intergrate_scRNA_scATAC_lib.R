# help: 
https://www.yuque.com/docs/share/5c170e2d-36aa-4819-825e-a0004b3346f1#

{
  library(pheatmap)
	library(dplyr)
	if(!require('Cairo')){  BiocManager::install('Cairo',ask = F,update = F)  }; library(Cairo)
	if(!require('ggplot2')){  BiocManager::install('ggplot2',ask = F,update = F)  }; library(ggplot2)
	if(!require('Seurat')){  BiocManager::install('Seurat',ask = F,update = F)  }; library(Seurat)
	if(!require('Signac')){  BiocManager::install('Signac',ask = F,update = F)  }; library(Signac)
	if(!require('EnsDb.Mmusculus.v75')){  BiocManager::install('EnsDb.Mmusculus.v75',ask = F,update = F)  }; library(EnsDb.Mmusculus.v75)
	library(EnsDb.Hsapiens.v86)
	library(GenomeInfoDb)
}
Outdir<- function(step2){
        if(!file.exists(step2)) {dir.create(step2,recursive=T)}
        setwd(step2)
}


COVERAGE_PLOT <- function(pbmc.atac.filtered,G_VER=EnsDb.Mmusculus.v75,gene_name,gene_region,ymax=NULL){
	cat("n---------------------------------------------------\n [START] COVERAGE_PLOT at ",date(),"\n")
  DefaultAssay(pbmc.atac.filtered)='ATAC'
	pdf(paste(gene_name,".pdf",sep=""),height=8*length(gene_region))
	print(CoveragePlot(object = pbmc.atac.filtered,ymax=ymax,region = gene_region,
	  group.by="predicted.id", sep = c("_", "_"), annotation = G_VER, extend.upstream = 20000,  extend.downstream = 20000, ncol = 1 ))
	dev.off()
	cat(" [END] COVERAGE_PLOT end at ",date(),"\n\n")
}
# Usage:
#   pbmc.atac.filtered=readRDS("scRNA_scATAC.combind.highscore.rds")
#   COVERAGE_PLOT(pbmc.atac.filtered,G_VER=EnsDb.Mmusculus.v75,'Cd8a','chr6_71373427_71379173'，ymax=1)


HEATMAP_PLOT <- function(seuratObj,gene_list){
	cat("n---------------------------------------------------\n [START] HEATMAP_PLOT at ",date(),"\n")
	data = DotPlot(seuratObj, features = rev(gene_list), cols = c("blue","red"), dot.scale = 8)$data
	new_data=data.frame()
	for(id in levels(as.factor(data$id))){
		#id = levels(as.factor(data$id))[2]
		for(gene in  levels(as.factor(data$features.plot)) ){
			# gene=levels(as.factor(data$features.plot))[1]
			new_data[gene,id]=data[(data$id==id) & (data$features.plot == gene),"avg.exp.scaled"]
		}
	}

	height=500
	if(length(gene_list)>30){height=1800}
	CairoPNG(file="heatmap_cellType.png",width=600,height=height)
	pheatmap(new_data,border_color=NA,cluster_cols=FALSE,color=colorRampPalette(c("blue", "white", "red"))(100)[c(20:80)])
	dev.off()
	
	gene_list=gene_list[c(1:20)]
	data = DotPlot(seuratObj, features = rev(gene_list), cols = c("blue","red"), dot.scale = 8)$data
	new_data=data.frame()
	for(id in levels(as.factor(data$id))){
		#id = levels(as.factor(data$id))[2]
		for(gene in  levels(as.factor(data$features.plot)) ){
			# gene=levels(as.factor(data$features.plot))[1]
			new_data[gene,id]=data[(data$id==id) & (data$features.plot == gene),"avg.exp.scaled"]
		}
	}

	CairoPNG(file="heatmap_cellType_top20.png",width=600,height=500)
	pheatmap(new_data,border_color=NA,cluster_cols=FALSE,color=colorRampPalette(c("blue", "white", "red"))(100)[c(20:80)])
	dev.off()	
	cat(" [END] HEATMAP_PLOT end at ",date(),"\n\n")
}

markerPlot <- function(seuratObj,marker_list,group,QUICK=TRUE){
	cat("n---------------------------------------------------\n [START] markerPlot at ",date(),"\n")
	# group="CellType"
	text=theme(text=element_text(size=20, color="black"))

	##introduction for parameters
	path_save <- getwd()
	colN <- 1:ncol(marker_list)
	marker_type <- colnames(marker_list)
	#initial the value
	gene_list <- data.frame()
	## start ploting
	#DefaultAssay(object = seuratObj) <- "RNA"
	for (counterCol in colN){
		Outdir(marker_type[counterCol])
		path.loop <- paste(path_save,marker_type[counterCol],sep="/")
		setwd(path.loop)
		rowN <- 1:length(marker_list[,counterCol])
		for(counterRow in rowN){
			gene_name <- as.character(marker_list[counterRow,counterCol])
			if(nchar(gene_name)>0 && !is.na(gene_name)){
				if( length( grep(pattern=paste("^",gene_name,"$",sep=""),x=rownames(x=seuratObj),ignore.case = TRUE) ) ){
					gene_name <- grep(pattern=paste("^",gene_name,"$",sep=""),x=rownames(x=seuratObj),ignore.case = TRUE, value=TRUE)
					#QUICK=TRUE
					if(QUICK){
						vinPlot_name <- paste(gene_name,"vin.png",sep="_")
						tsnePlot_name <- paste(gene_name,"tsne.png",sep="_")
						umapPlot_name <- paste(gene_name,"umap.png",sep="_")
						CairoPNG(file=umapPlot_name,width=700,height=600)
						print(FeaturePlot(object = seuratObj, features = gene_name, reduction = "umap") + text	)	
						dev.off()
						CairoPNG(file=vinPlot_name,width=700,height=600)
						print(VlnPlot(object = seuratObj, features = gene_name, slot = "counts", log = TRUE,group.by=group) + text)
						dev.off()
						CairoPNG(file=tsnePlot_name,width=700,height=600)
						print(FeaturePlot(object = seuratObj, features = gene_name, reduction = "tsne") + text)
						dev.off()
					}
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

	#text=theme(text=element_text(size=40, color="black"))
	## summary gene dot plot
	gene_list <- unique(as.character(gene_list))
	write.csv(gene_list,file='gene_plot.csv',row.names=FALSE)
	
	HEATMAP_PLOT(seuratObj,gene_list)
	CairoPNG(file = "DotPlot_allMarker.png",width=200+40*length(gene_list),height=100+30*length(levels(seuratObj$seurat_clusters)))
	print(DotPlot(seuratObj, features = rev(gene_list), cols = c("blue","red"), dot.scale = 8)+ RotatedAxis() + text )
	dev.off()
	
	gene_list_20=gene_list[c(1:20)]
	CairoPNG(file = "DotPlot_allMarker_top20.png",width=200+40*length(gene_list_20),height=100+30*length(levels(seuratObj$seurat_clusters)))
	print(DotPlot(seuratObj, features = rev(gene_list_20), cols = c("blue","red"), dot.scale = 8)+ RotatedAxis() + text )
	dev.off()
	
	if(length(levels(as.factor(seuratObj$orig.ident)))>1) {
		CairoPNG(file = "DotPlot_allSample_allMarker.png",width=200+30*length(gene_list),height=100+40*( length(levels(seuratObj$seurat_clusters))+length(levels(as.factor(seuratObj$orig.ident))) ) )
		print(DotPlot(seuratObj, features = rev(gene_list), cols = c("blue","red","green","purple","coral","gold2","deepskyblue","darkorange"), dot.scale = 8, split.by = "orig.ident") + RotatedAxis() + text)
		dev.off()
		
		CairoPNG(file = "DotPlot_allSample_allMarker_top20.png",width=200+30*length(gene_list_20),height=100+40*( length(levels(seuratObj$seurat_clusters))+length(levels(as.factor(seuratObj$orig.ident))) ) )
		print(DotPlot(seuratObj, features = rev(gene_list_20), cols = c("blue","red","green","purple","coral","gold2","deepskyblue","darkorange"), dot.scale = 8, split.by = "orig.ident") + RotatedAxis() + text)
		dev.off()
	}
	cat("markerPlot work has done!!!","\n")
	cat(" [END] markerPlot end at ",date(),"\n\n")
}

markerPlot_gene_name <- function(gene_name,seuratObj,group){
	cat("n---------------------------------------------------\n [START] markerPlot_gene_name at ",date(),"\n")
	# group="CellType"
	text=theme(text=element_text(size=20, color="black"))

	gene_name <- grep(pattern=paste("^",gene_name,"$",sep=""),x=rownames(x=seuratObj),ignore.case = TRUE, value=TRUE)
	vinPlot_name <- paste(gene_name,"vin.png",sep="_")
	tsnePlot_name <- paste(gene_name,"tsne.png",sep="_")
	umapPlot_name <- paste(gene_name,"umap.png",sep="_")
	CairoPNG(file=umapPlot_name,width=700,height=600)
	print(FeaturePlot(object = seuratObj, features = gene_name, reduction = "umap") + text	)	
	dev.off()
	CairoPNG(file=vinPlot_name,width=700,height=600)
	print(VlnPlot(object = seuratObj, features = gene_name, slot = "counts", log = TRUE,group.by=group) + text)
	dev.off()
	CairoPNG(file=tsnePlot_name,width=700,height=600)
	print(FeaturePlot(object = seuratObj, features = gene_name, reduction = "tsne") + text)
	dev.off()

	cat("markerPlot work has done!!!","\n")
	cat(" [END] markerPlot_gene_name end at ",date(),"\n\n")
}

TSNE_UMAP_Plot_CellType_tmp <- function(Step2,pbmc_atac_filtered,pbmc_rna){
	cat("n---------------------------------------------------\n [START] TSNE_UMAP_Plot_CellType_tmp at ",date(),"\n")
	#  umap
	p1 <- DimPlot(pbmc_atac_filtered, reduction = "umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + scale_colour_hue(drop = FALSE)
	p2 <- DimPlot(pbmc_rna, reduction = "umap", group.by = "CellType", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") 
		
	file1_png=paste(Step2,"/umap_scRNA_scATAC.CellType.png",sep="")
	CairoPNG(file=file1_png,width=480*2.5,height=480)
	print(CombinePlots(plots = list(p1, p2)))
	dev.off()

	#  tsne
	p1 <- DimPlot(pbmc_atac_filtered,reduction = "tsne",  group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + scale_colour_hue(drop = FALSE)
	p2 <- DimPlot(pbmc_rna,reduction = "tsne",  group.by = "CellType", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") 
		
	file1_png=paste(Step2,"/tsne_scRNA_scATAC.CellType.png",sep="")
	CairoPNG(file=file1_png,width=480*2.5,height=480)
	print(CombinePlots(plots = list(p1, p2)))
	dev.off()
	cat(" [END] TSNE_UMAP_Plot_CellType_tmp end at ",date(),"\n\n")
 }

markerFind <- function(seuratObj,test_method="MAST",min.pct=0.25,logFC=0.25) {

#修改坐标轴刻度标签
	axis.text=theme(axis.text = element_text(size = 20, color = "black",angle=0,hjust=1))
	#修改坐标轴标签
	axis.title = theme(axis.title = element_text(size = 20, color = "black", face = "bold"))
	#对legend的内容做修改
	legend.text=theme(legend.text= element_text(size=20, color="black",  vjust=0.5, hjust=0.5))
	#对legend的title做修改
	legend.title=theme(legend.title= element_text(size=20, color="black", face = "bold", vjust=0.5, hjust=0.5))

	modify = axis.text + axis.title + legend.text + legend.title
	
	if( min.pct>=0 && min.pct<=1 ) {
		# ClusterNum <- levels(seuratObj$seurat_clusters)
		# for (numC in ClusterNum ){
		# 	name.cluster = paste("cluster",numC,test_method,"markers.csv",sep="_")
		#   	cluster.markers <- FindMarkers(object = seuratObj, ident.1 = numC, test.use = test_method, min.pct = min.pct, logfc.threshold=logFC)
		#   	write.csv(cluster.markers, file = name.cluster)
		# }
		# find markers for every cluster compared to all remaining cells
		# only the positive ones
		seurat.markers <- FindAllMarkers(object = seuratObj, only.pos = FALSE, test.use = test_method, min.pct = min.pct, logfc.threshold=logFC)
		write.csv(data.frame(Name=rownames(seurat.markers),seurat.markers), file = "AllMarkers.csv",row.names=F)
		# seurat.markers <- read.csv(file='AllMarkers.csv',header=T,row.names=1)
		# plot top gene heatmap
		seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
		#top10 <- seurat.markers %>% group_by(cluster) %>% arrange((p_val_adj),desc(avg_logFC)) %>% top_n(n = 10, wt = avg_logFC)
		#CairoPNG(file="Top10Heatmap.png",width=50*length(levels(seuratObj$seurat_clusters)),height=50*length(levels(seuratObj$seurat_clusters)))
		pdf(file="Top10Heatmap.pdf",width=1.5*length(levels(seuratObj$seurat_clusters)),height=1.5*length(levels(seuratObj$seurat_clusters)))
		print(DoHeatmap(object = seuratObj, features = top10$gene, size=20) )
		dev.off()
		
		#------------------------------------------
		# wlp add
		#------------------------------------------
		level=0
		tmp=data.frame(subset(top10,cluster==level)$gene)
		colnames(tmp)=level
		top_list=tmp
		for(level in levels(top10$cluster)[c(2:length(levels(top10$cluster)))]){
			tmp=data.frame(subset(top10,cluster==level)$gene)
			colnames(tmp)=level
			top_list=cbind(top_list,tmp)
		}
		write.csv(top_list, file = "top10_list.csv",row.names=F)
		
		seurat.markers2=subset(seurat.markers,seurat.markers[,4]<0.2)
		
		seurat.markers2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10_2

		top_list=data.frame()
		for(level in levels(top10_2$cluster)){
			tmp=data.frame(subset(top10_2,cluster==level)$gene)
			colnames(tmp)=level
			top_list[rownames(tmp),colnames(tmp)]=as.character(tmp[,1])
		}
		
		write.csv(top_list, file = "top10_specific_list.csv",row.names=F)
		
		#------------------------------------------
		# wlp add
		#------------------------------------------
	}
	else{
		print("Error!")
		print("Please check input data. The input data may out of range or with incorrect format!")
	}
	time_log=""
	#return(seuratObj)
	TIME_END(paste("AllMarkers.csv","\n\t\t","Top10Heatmap.png","\n\t\t","top10_list.csv","\n\t\t","top10_specific_list.csv"),time_log)
	TIME_END(getwd(),time_log)
	TIME_END("markerFind",time_log)
}


















