source("library/scRNASeq_lib.R")
#------------------------------------------------------------
# 0.input
config <- data.frame(testMethod='MAST',PCs=50,clusRes=0.2,species='mouse')#,sampleCounter=sampleNum)
time_log=""
#------------------------------------------------------------
pathHome="combind"
sampleInfo="sampleSheet2.csv"
markerFile <- "marker_list.csv"
clusterfile="cluster.txt"
#------------------------------------------------------------
# outdir
Outdir(pathHome)

#------------------------------------------------------------
# 1.read data
#------------------------------------------------------------
sampleSheet <- read.csv(file=sampleInfo, header=T,colClasses="character")
# sampleName <- as.character(sampleSheet$sampleName)
# sampleGroup <- as.character(sampleSheet$sampleGroup)
# sampleData <- as.character(sampleSheet$sampleData)

pathData=as.character(sampleSheet$SeuratData)

time_log=TIME_START("Read",time_log)
# impordata
sample_info=list()
for(i in c(1:length(pathData))){
  # i=2
  sample1Data <- read.table(file = pathData[i], header=T, row.names=1, sep = "\t")
  sample1 <- CreateSeuratObject(counts = sample1Data, project = sampleSheet$sampleName[i], min.cells = 3,min.features = 200)
  sample1$stim <- sampleSheet$sampleGroup[i]
  sample1 <- subset(x = sample1, subset = nFeature_RNA > 200)
  sample1 <- NormalizeData(object = sample1, verbose = FALSE)
  sample1 <- FindVariableFeatures(object = sample1, selection.method = "vst", nfeatures = 2000)
  sample_info[sampleSheet$sampleName[i]]=sample1
}
time_log=TIME_END("Read",time_log)


#----------------------------------------------
# step1: 001.Combined_Cluster
#----------------------------------------------
time_log=TIME_START("step1",time_log)

step1 <- paste(pathHome,'/001.Combined_Cluster/',sep="")
Outdir(step1)
#### Perform integration
sample_anchors <- FindIntegrationAnchors(object.list = sample_info, dims = 1:config$PCs)
seuratObj <- IntegrateData(anchorset = sample_anchors, dims = 1:config$PCs)

#-----------save----------------
nameObj = 'sampleCombined_IntegrateData.rds'
saveRDS(seuratObj,file = paste(pathHome,nameObj,sep='/'))
# seuratObj=readRDS(paste(pathHome,nameObj,sep='/'))

#----------------------------------------------
#                  integrated                 
#----------------------------------------------
# Perform an integrated analysis
DefaultAssay(object = seuratObj) <- "integrated"
# Run the standard workflow for visualization and clustering
seuratObj <- ScaleData(object = seuratObj, verbose = FALSE, do.scale = TRUE, do.center = TRUE)
seuratObj <- RunPCA(object = seuratObj, npcs = config$PCs, verbose = FALSE)
# t-SNE and Clustering
seuratObj <- RunUMAP(object = seuratObj, reduction = "pca", dims = 1:config$PCs)
seuratObj <-RunTSNE(object = seuratObj,reduction='pca',dims=1:config$PCs,tsne.method='Rtsne',features=rownames(seuratObj))
seuratObj <- FindNeighbors(object = seuratObj, reduction = "pca", dims = 1:config$PCs)
#seuratObj <- FindClusters(seuratObj, resolution = config$clusRes,resolution = 0.1)
seuratObj <- FindClusters(seuratObj,resolution = config$clusRes)

#-----------save----------------
nameObj = "sampleCombined_Cluster.rds"
saveRDS(seuratObj,file = paste(pathHome,nameObj,sep='/'))
# seuratObj=readRDS(paste(pathHome,nameObj,sep='/'))

#----------------------------------------------
# get cell count, and visualization
#----------------------------------------------
if(length(levels(as.factor(seuratObj@meta.data$stim))) > 1) {
  metaData<- data.frame(barcode=rownames(seuratObj@meta.data),sample=seuratObj@meta.data$orig.ident,group=seuratObj@meta.data$stim, nGene=seuratObj@meta.data$nFeature_RNA,nUMI=seuratObj@meta.data$nCount_RNA,cluster=seuratObj@meta.data$seurat_clusters)
  metaData$barcode <- str_sub(metaData$barcode,end=16)
  write.csv(metaData,file='filtered_meta_data.csv',row.names=F)

  # count the cell number in subCluster and sample
  count_cell_all=TMP_SAVE(metaData,"cluster")
  write.csv(data.frame(clusterID=rownames(count_cell_all),count_cell_all),file='cellNum_InCluster.csv',row.names=F)
  PIE_PLOT(paste(step1,'/cellNum_InCluster.csv',sep=""))
}else{
  metaData<- data.frame(barcode=rownames(seuratObj@meta.data),sample=seuratObj@meta.data$orig.ident, nGene=seuratObj$nFeature_RNA,nUMI=seuratObj$nCount_RNA,cluster=seuratObj$seurat_clusters)
  metaData$barcode <- str_sub(metaData$barcode,end=16)
  write.csv(metaData,file='filtered_meta_data.csv',row.names=F)

  # count the cell number in subCluster and sample
  count_cell <- table(metaData[,c(2,5)])
  write.csv(t(count_cell),file='cellNum_InCluster.csv',row.names=T)
  count_cell <- read.csv(file='cellNum_InCluster.csv',header=T)
  colnames(count_cell)[1]='clusterID'
  
  # count_cell$sum <- rowSums(count_cell[,2:ncol(count_cell)])
  write.csv(count_cell,file='cellNum_InCluster.csv',row.names=F)
}


# plot tsne and umap
TSNE_UMAP_PLOT(seuratObj,outdir=step1,group="seurat_clusters")


#----------------------------------------------
# Differential analysis
#----------------------------------------------
time_log=TIME_START("markerFind",time_log)
markerFind(seuratObj,as.character(config$testMethod))
time_log=TIME_END("step1",time_log)

#------------------------------------------------------------
# DEGSeq Gene visualization
#------------------------------------------------------------
## marker plot
step2 <- paste(pathHome,'/002.Marker_Plot_v4/',sep="")

DefaultAssay(object = seuratObj) <- "RNA"
Marker_Plot(markerFile,step2,step1,seuratObj,"seurat_clusters")

#------------------------------------------------------------
# step4: add CellType
#------------------------------------------------------------
step4=paste(pathHome,"/004.CellType/",sep="")
Outdir(step4)
seuratObj=RENAME_CELL_TYPE(seuratObj,pathHome,clusterfile)

Marker_Plot(markerFile,step4,step1,seuratObj,"CellType")

##-------------- add CellType----------------------------
Outdir(step4)
metaData=data.frame(metaData,CellType=seuratObj@meta.data$CellType)
count_cell_all=TMP_SAVE(metaData,"CellType")
write.csv(data.frame(clusterID=rownames(count_cell_all),count_cell_all),file=paste(step4,'/cellNum_InCluster.csv',sep=""),row.names=F)
PIE_PLOT(paste(step4,'/cellNum_InCluster.csv',sep=""))

#------------------------------------------------------------
# step3: 003.Monocle_Analysis
#------------------------------------------------------------
step3 <- paste(pathHome,'/003.Monocle_Analysis/',sep="")
Outdir(step3)
monocleObj <- monoclePlot_wlp(seuratObj,config$clusRes)
# monocleObj=list(my_cds_subset,my_cds)
#################################################################################
nameObj = paste("sampleCombined","monocle.rds", sep = "_")
Outdir(pathHome)
saveRDS(monocleObj,file = paste(pathHome,nameObj,sep='/'))## step3 monocle end

#------------------------------------------------------------
# step5 & step6, GO and KEGG analysis
#------------------------------------------------------------
step5 <- paste(pathHome,'/005.GO_Analysis/',sep="")
step6 <- paste(pathHome,'/006.Pathway_Analysis/',sep="")
Outdir(step5)
Outdir(step6)

#################################################
seurat_all_marker<-read.csv(file=paste(step1,'AllMarkers.csv',sep='/'),header=T)
clusterNum <- levels(as.factor(seurat_all_marker$cluster))
markerSubset <- seurat_all_marker[seurat_all_marker$p_val_adj < 0.05,]
for (count in clusterNum[c(2:length(clusterNum))]) {
  ## get cluster name
  #clusterName <- str_extract(file,'cluster_[0-9]*')
  #markerData <- read.csv(file=paste(path.step2,file,sep='/'),header=T,row.names=1)
  ## filt for all gene: p_val_adj < 0.05
  clusterName <- paste('cluster',count,sep='_')
  markerData <- markerSubset[markerSubset$cluster == count,]
  if (nrow(markerData) > 10){
    ## GO analysis
    Outdir(paste(step5,"/",clusterName,sep=""))
    getGOTest(geneFC=markerData,sampleName=clusterName,organism=as.character(config$species),type='all')
    ## KEGG analysis
    Outdir(paste(step6,"/",clusterName,sep=""))
    getKEGGTest(geneFC=markerData,sampleName=clusterName,organism=as.character(config$species),type='all')
    #getPathview()
  }else{
    cat('the all gene number is too low !','\t',clusterName,'\n')
  }
  ## filter for up gene:  avg_logFC > 1
  markerUP <- markerData[markerData$avg_logFC > 1,]
  if (nrow(markerUP) > 5){
    ## GO analysis
    Outdir(paste(step5,"/",clusterName,sep=""))
    getGOTest(geneFC=markerUP,sampleName=clusterName,organism=as.character(config$species),type='up')
    ## KEGG analysis
    Outdir(paste(step6,clusterName,sep=""))
    getKEGGTest(geneFC=markerUP,sampleName=clusterName,organism=as.character(config$species),type='up')
  }else{
    cat('the up gene number is too low !','\t',clusterName,'\n')
  }
  ## filt for down gene:  avg_logFC > 1
  markerDOWN <- markerData[markerData$avg_logFC < -1,]
  if (nrow(markerDOWN) > 5){
    ## GO analysis
    Outdir(paste(step5,"/",clusterName,sep=""))
    getGOTest(geneFC=markerDOWN,sampleName=clusterName,organism=as.character(config$species),type='down')
    
    ## KEGG analysis
    Outdir(paste(step6,"/",clusterName,sep=""))
    getKEGGTest(geneFC=markerDOWN,sampleName=clusterName,organism=as.character(config$species),type='down')
  }else{
    cat('the down gene number is too low !','\t',clusterName,'\n')
  }
}

