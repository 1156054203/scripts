
{
  library(dplyr)
  if(!require('Cairo')){  BiocManager::install('Cairo',ask = F,update = F)  }; library(Cairo)
  if(!require('ggplot2')){  BiocManager::install('ggplot2',ask = F,update = F)  }; library(ggplot2)
  if(!require('Seurat')){  BiocManager::install('Seurat',ask = F,update = F)  }; library(Seurat)
  if(!require('Signac')){  BiocManager::install('Signac',ask = F,update = F)  }; library(Signac)
}

vlnPlot_corPlot_tmp <- function(Step1,atac){
  cat("\n---------------------------------------------------\n[START]vlnPlot_corPlot_tmp at ",date(),"\n")
  Outdir(Step1)
  features = c("nCount_ATAC","nFeature_ATAC","nCount_ACTIVITY","nFeature_ACTIVITY","passed_filters","peak_region_fragments","peak_region_cutsites")
  
  file1_png=paste(Step1,"/vlnPlot.png",sep="")
  file1_pdf=paste(Step1,"/vlnPlot.pdf",sep="")
  plot1=VlnPlot(object = atac, features =features ,pt.size=0.01,ncol=4,group.by="orig.ident")

  CairoPNG(file=file1_png,width=480*4,height=480*2)
  print(plot1)  
  ggsave(file=file1_pdf,width=7*6,height=7*1)
  dev.off()
  cat(" [END] vlnPlot_corPlot_tmp end at ",date(),"\n\n")
}

Extend_self <- function (x, upstream = 0, downstream = 0) 
{
    if (any(strand(x = x) == "*")) {
        warning("'*' ranges were treated as '+'")
    }
    on_plus <- strand(x = x) == "+" | strand(x = x) == "*"

    new_start <-  ifelse(test = on_plus, yes = start(x = x) -upstream, 
        no = end(x = x))
    new_end <-  ifelse(test = on_plus, yes = start(x = x), 
        no = end(x = x) +upstream)
    
  ranges(x = x) <- IRanges(start = new_start, end = new_end)
    x <- trim(x = x)
    return(x)
}

EXT_COORDS <- function(G_VER=EnsDb.Hsapiens.v86,fragment_path,pbmc_red,region="activity"){
  cat("\n---------------------------------------------------\n [START] EXT_COORDS ","end at ",date(),"\n\n")
  cat("------G_VER:\n")
  print(G_VER)
  cat("------fragment_path:\n")
  print(fragment_path)
  cat("------pbmc_red:\n")
  print(pbmc_red)
  cat("------region:\n")
  print(region)
  # extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object 
  gene.coords <- genes(G_VER, filter = ~ gene_biotype == "protein_coding")
  seqlevelsStyle(gene.coords) <- 'UCSC'
  
  # create a gene by cell matrix

  if(region=="activity"){
    # convert rownames from chromsomal coordinates into gene names
    genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)
    gene_activities <- FeatureMatrix( fragments = fragment_path, features = genebodyandpromoter.coords, cells = colnames(pbmc_red), chunk = 10 )
    gene_key <- genebodyandpromoter.coords$gene_name
    names(gene_key) <- GRangesToString(grange = genebodyandpromoter.coords)
    rownames(gene_activities) <- gene_key[rownames(gene_activities)]
  }
  if(region=="genebody"){
    genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
    gene_activities <- FeatureMatrix( fragments = fragment_path,  features = genebody.coords, cells = colnames(pbmc_red), chunk = 10)
    genebody_key <- genebody.coords$gene_name
    names(genebody_key) <- GRangesToString(grange = genebody.coords)
    rownames(gene_activities) <- genebody_key[rownames(gene_activities)]
  }
  if(region=="promoter"){
    promoter.coords <- Extend_self(x = gene.coords, upstream = 2000, downstream = 0)
    gene_activities <- FeatureMatrix( fragments = fragment_path, features = promoter.coords, cells = colnames(pbmc_red), chunk = 10)
    promoter_key <- promoter.coords$gene_name
    names(promoter_key) <- GRangesToString(grange = promoter.coords)
    rownames(gene_activities) <- promoter_key[rownames(gene_activities)]
  }
  cat(" [END] EXT_COORDS end at ",date(),"\n\n")
  return(gene_activities)
}

Outdir<- function(step2){
        if(!file.exists(step2)) {dir.create(step2,recursive=T)}
        setwd(step2)
}

markerFind_scRNA_scATACSeq <- function(seuratObj,Step3,calculate="yes") {
  cat("\n---------------------------------------------------\n[START]markerFind_scRNA_scATACSeq at ",date(),"\n")
#修改坐标轴刻度标签
  axis.text=theme(axis.text = element_text(size = 20, color = "black",angle=0,hjust=1))
  #修改坐标轴标签
  axis.title = theme(axis.title = element_text(size = 20, color = "black", face = "bold"))
  #对legend的内容做修改
  legend.text=theme(legend.text= element_text(size=20, color="black",  vjust=0.5, hjust=0.5))
  #对legend的title做修改
  legend.title=theme(legend.title= element_text(size=20, color="black", face = "bold", vjust=0.5, hjust=0.5))

  modify = axis.text + axis.title + legend.text + legend.title
  
  if( 1 ) {
    #calculate="yes"
    file_markers=paste(Step3,"/seurat_markers.csv",sep="")
    if(calculate == "yes"){
      seurat_markers <- FindAllMarkers(object = seuratObj, only.pos = F, test.use = 'LR')
      write.csv(data.frame(Peak=rownames(seurat_markers),seurat_markers), file = file_markers,row.names=F)
    }else{
      seurat_markers=read.csv(file_markers,row.names=1)
    }
    
    seurat_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
    #seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
    #CairoPNG(file="Top10Heatmap.png",width=50*length(levels(seuratObj$seurat_clusters)),height=50*length(levels(seuratObj$seurat_clusters)))
    pdf(file=paste(Step3,"/Top10Heatmap.pdf",sep=""),width=1.5*length(levels(seuratObj$seurat_clusters)),height=1.5*length(levels(seuratObj$seurat_clusters)))
    print(DoHeatmap(object = seuratObj, features = as.character(top10$gene), size=10) )
    dev.off()
    
    #------------------------------------------
    # wlp add: save top10_list as marker for plot scatter
    #------------------------------------------
    if(length(levels(top10$cluster))==0){
      level=0
      tmp=data.frame(subset(top10,cluster==level)$gene)
      colnames(tmp)=level
      top_list=tmp
      for(level in unique(top10$cluster)[c(2:length(unique(top10$cluster)))]){
        tmp=data.frame(subset(top10,cluster==level)$gene)
        colnames(tmp)=level
        top_list=cbind(top_list,tmp)
      }
    }else{
      level=0
      tmp=data.frame(subset(top10,cluster==level)$gene)
      colnames(tmp)=level
      top_list=tmp
      for(level in levels(top10$cluster)[c(2:length(levels(top10$cluster)))]){
        tmp=data.frame(subset(top10,cluster==level)$gene)
        colnames(tmp)=level
        top_list=cbind(top_list,tmp)
      }
    }
    write.csv(top_list, file = paste(Step3,"/top10_list.csv",sep=""),row.names=F)

    # wlp add: save top10_list specif expressed in each cluster as marker for plot scatter

    seurat.markers2=subset(seurat_markers,seurat_markers[,4]<0.2)
    
    #seurat.markers2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10_2

    #top_list=data.frame()
    #for(level in levels(top10_2$cluster)){
    # tmp=data.frame(subset(top10_2,cluster==level)$gene)
    # colnames(tmp)=level
    # top_list[rownames(tmp),colnames(tmp)]=as.character(tmp[,1])
    #}
    #write.csv(top_list, file = paste(Step3,"/top10_specific_list.csv",sep=""),row.names=F)
    
    #------------------------------------------
    # wlp add
    #------------------------------------------
  }
  else{
    print("Error!")
    print("Please check input data. The input data may out of range or with incorrect format!")
  }
  cat(" [END] markerFind_scRNA_scATACSeq end at ",date(),"\n\n")
}

PLOT_red <- function(pbmc_red,outdir,resolution=0.8){
  cat("n---------------------------------------------------\n [START] PLOT_red at ",date(),"\n\n")
  p1_1=DimPlot(object = pbmc_red, label = TRUE,group.by="seurat_clusters",reduction="lsi")
  p2_1=DimPlot(object = pbmc_red, label = TRUE,group.by="seurat_clusters",reduction="umap")
  p3_1=DimPlot(object = pbmc_red, label = TRUE,group.by="seurat_clusters",reduction="tsne")

  p1_2=DimPlot(object = pbmc_red, label = TRUE,group.by="orig.ident",reduction="lsi")
  p2_2=DimPlot(object = pbmc_red, label = TRUE,group.by="orig.ident",reduction="umap")
  p3_2=DimPlot(object = pbmc_red, label = TRUE,group.by="orig.ident",reduction="tsne")
  p1 <- CombinePlots(plots = list(p1_1,p2_1,p3_1), ncol = 3)
  p2 <- CombinePlots(plots = list(p1_2,p2_2,p3_2), ncol = 3)
  
  file_png=paste(outdir,"/03_reduction_",resolution,".png",sep="")
  file_pdf=paste(outdir,"/03_reduction_",resolution,".pdf",sep="")
  CairoPNG(file=file_png,width=480*3,height=480*2*0.8)
  print(CombinePlots(list(p1,p2),ncol = 1))
  ggsave(file=file_pdf,width=7*3,height=7*2*0.8)
  dev.off()
  cat(" [END] PLOT_red ",file_png,"end at ",date(),"\n\n")
}

ADD_SAMPLE <- function(sample_name_order,pbmc_red){
  cat("n---------------------------------------------------\n [START] ADD_SAMPLE at ",date(),"\n\n")
  #ident
  if(length(sample_name_order)>0){
    merge_name=unlist(strsplit(sample_name_order,","))
    ident=sub("^[A-Z]+-","",colnames(pbmc_red))
    ident_tmp=factor(ident,levels=c(1:length(merge_name)),labels=merge_name)
    pbmc_red$orig.ident <- ident_tmp
  }
  return(pbmc_red)
  cat(" [END] ADD_SAMPLE end at ",date(),"\n\n")
}

PLOT_FRAG <-function(pbmc,outdir){
  cat("\n---------------------------------------------------\n[START]PLOT_FRAG at ",date(),"\n")
  file_png=paste(outdir,"/02_fragment_size.png",sep="")
  file_pdf=paste(outdir,"/02_fragment_size.pdf",sep="")
  CairoPNG(file=file_png,width=480*2*0.8,height=480*0.8*0.8)
  print(PeriodPlot(object = pbmc, group.by = 'nucleosome_group',region='chr1-1-200000000'))
  ggsave(file=file_pdf,width=7*0.8*2,height=7*0.8*0.8)
  dev.off()
  dev.off()
  cat(" [END] PLOT_FRAG ",file_png," end at ",date(),"\n\n")
}

PLOT_QC <- function(pbmc,type,outdir){
  cat("n---------------------------------------------------\n [START] PLOT_QC at ",date(),"\n")
  plot1 <- VlnPlot(object = pbmc, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal'), pt.size = 0.01) + NoLegend()

  plot2_a <- VlnPlot( object = pbmc, features = 'peak_region_fragments', pt.size = 0.01, log = TRUE) + NoLegend()
  plot2_b <- FeatureScatter(pbmc,"peak_region_fragments",'nucleosome_signal', pt.size = 0.01) + NoLegend()
  plot2_c <- FeatureScatter(pbmc,"peak_region_fragments",'blacklist_ratio', pt.size = 0.01) + NoLegend()
  plot2 <- CombinePlots(plots = list(plot2_a,plot2_b,plot2_c), ncol = 3)

  file_png=paste(outdir,"/01_","QC_",type,".png",sep="")
  file_pdf=paste(outdir,"/01_","QC_",type,".pdf",sep="")
  CairoPNG(file=file_png,width=480*3*0.8,height=480*2*0.8)
  print(CombinePlots(list(plot1,plot2),ncol = 1))
  ggsave(file=file_pdf,width=7*3*0.8,height=7*2*0.8)
  dev.off()
  cat(" [END] PLOT_QC ",file_png,"end at ",date(),"\n\n")
}




