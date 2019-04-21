library(pheatmap)
#library(pheatmap)
#pdf("pheatmap.pdf",width=6,height=15)
#png("pheatmap.png",width=480,height=480)                  
#pdf("pheatmap.pdf",width=6,height=15)
p <- read.table("test.xls",header=T,sep="\t")
rownames(p)=p$EnsemblGene_GeneSymbol
colnames(p)=head(p,1)
p <- p[,-1]
p_matrix <- as.matrix(p)
p_matrix <- log10(p_matrix + 1)

annotation_col=data.frame(Sample=factor(c("S1","S2","S3","S4","S5","S6")),time=1:6) 
rownames(annotation_col)=colnames(p)
  
annotation_row=data.frame(GeneClass=factor(rep(c("Path1","Path2","Path3"),c(9,2,2)))) 
rownames(annotation_row)=rownames(p)

ann_colors=list(Time=c("white","firebrick"),Sample=c(S1="#1B9E77",S2="#D95F02",S3="red",S4="blue",
        S5="green",S6="yellow"),GeneClass=c(Path1="#7570B3",Path2 ="#E7298A",Path3 ="#66A61E"))

pheatmap(p_matrix,color=colorRampPalette(c("blue","white","red")) (100) ,
         scale="column",border_color="grey60",show_rownames=T,show_colnames=T,
         clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",clustering_method="complete",
         fontsize_col=6,fontsize_row=6,fontsize_number=8,fontsize=9,main="Test pheatmap",legend=T,
         display_numbers=T,number_color="black",annotation_colors=ann_colors,cellwidth=40,
         cellheight=15,annotation_col=annotation_col,annotation_row=annotation_row,
         annotation_legend=T,cutree_rows=3,cutree_cols=2,filename="C:/Users/chen_yulong/Desktop/R/Example_pheatmap.png")
#dev.off()
