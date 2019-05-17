argv=commandArgs(TRUE)
#argv=c("FC.txt","pdf")
if(length(argv)==0){
  cat("USAGE:\n")
  cat(" $0 <common_pathway> <metabolites.FC> <protein.FC>\n\n")
  q()
  n
}
library(ggplot2)

file1=argv[1]
file2=argv[2]
file3=argv[3]

library(pathview)

data=read.table(file1,head=T,sep="\t")
exp_cpd=read.table(file2,head=T,sep="\t",quote="\"")
exp_gene=read.table(file3,head=T,sep="\t")

exp_gene=subset(exp_gene,grepl("K",KO))
exp_gene=exp_gene[,c(4,2)]
exp_cpd=exp_cpd[,c(1,3)]

exp_gene1=exp_gene[,2]
names(exp_gene1)=exp_gene[,1]

exp_cpd1=exp_cpd[,2]
names(exp_cpd1)=exp_cpd[,1]

pathway=as.character(data[,1])
pathway=sub("ko","",pathway)

pv=pathview(gene.data=exp_gene1,cpd.data=exp_cpd1,pathway.id=pathway,
	cpd.idtype="kegg",gene.idtype ="kegg",kegg.native = T,
	species = "ko",kegg.dir ="/home/R05/bin/Metabolics_Proteomics/KEGG_Orthologs",
	low = list(gene = "green", cpd = "blue"), 
	mid = list(gene = "gray", cpd = "gray"), 
	high = list(gene = "red", cpd =	"orange"))
url_all=c()
for (p in names(pv)){
	url="https://www.kegg.jp/kegg-bin/show_pathway?";

	col1=subset(pv[[p]][[1]],mol.col!="#FFFFFF")
	col2=subset(pv[[p]][[2]],mol.col!="#FFFFFF")
	col=rbind(col1,col2)
	col[,"mol.col"]=sub("#","",col[,"mol.col"])
	col=paste(col[,"kegg.names"],col[,"mol.col"],sep="%09%23")
	col=paste(col,collapse="/")
	p=sub("ko","map",p)
	url=paste(url,p,"/",sep="")
	url=paste(url,col,sep="")
	url_all[p]=url

}
write.table(url_all,"url.list",sep="",quote=F,col.names=F,row.names=F)
pathview(gene.data=exp_gene1,cpd.data=exp_cpd1,pathway.id=pathway,
	cpd.idtype="kegg",gene.idtype ="kegg",kegg.native = F,sign.pos="topright", 
	species = "ko",kegg.dir ="/home/R05/bin/Metabolics_Proteomics/KEGG_Orthologs",
	low = list(gene = "green", cpd = "blue"), 
	mid = list(gene = "gray", cpd = "gray"), 
	high = list(gene = "red", cpd =	"yellow"))
