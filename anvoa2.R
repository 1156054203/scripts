#install.packages("multcomp")
library(multcomp)
library(impute)
library(pheatmap)
library(car)

setwd("C:\\Users\\admin\\Desktop\\GSE21422\\healthy_DCIS_DIC\\Diff\\Diff_mRNA")
data_int <- read.table("GSE21422_series_matrix.txt",header=T,sep="\t",row.names = 1)
data_int <- data_int[rowMeans(data_int)>0.001,]

#数据准备
standard <- as.data.frame(t(data_int)+1)

#standalization
standard2 <- as.data.frame(scale(standard))
sta_mean <- as.data.frame(sapply(standard2,mean)) #均值都非常小
sta_sd <- as.data.frame(sapply(standard2,sd)) #标准差都等于1
standard =as.matrix(standard)
standard=t(standard)
write.table(standard,file = "probeid_stand_expr.txt",quote = F,sep = "\t",row.names = T)


#标准正态分布检验
par(mfrow=c(2,1))
qqnorm(standard2[,1],main="qq图")
qqline(standard2[,1])
hist(standard2[,1])

#主成分分析
pdf(file="screeplot.pdf")
standard2.pca <- prcomp(standard2)
(pca_result <- summary(standard2.pca))
#standard2 deviation代表每个主成分的标准差
#Proportion of Variance代表每个主成分的贡献率
#Cumulative Proportion代表各个主成分的累积贡献率
screeplot(standard2.pca, type="lines") #碎石图
dev.off()

#单因素多元方差分析,两个前提假设，一个是多元正态性，一个是方差—协方差矩阵同质性
condition <- factor(c(rep("healthy",5),rep("dcis",9),rep("dic",5)))
standard2 <- as.matrix(standard2)
fit <- manova(standard2~condition) #检验组间差异
fit1 <- summary.aov(fit)

#edit(fit1)查看fit1数据结构
#提取结果并输出
probeid <- colnames(standard2)
ncols <- ncol(standard2)
Df <- c()
ReDf <- c()
SumSq <- c()
ReSum <- c()
MeanSq <- c()
ReMean <- c()
fvalue <- c()
pvalue <- c()
for(i in 1:ncols)
{ Df[i] <- fit1[[i]][1,1]
  ReDf <- fit1[[i]][2,1]
  SumSq[i] <- fit1[[i]][1,2]
  ReSum[i] <- fit1[[i]][2,2]
  MeanSq[i] <- fit1[[i]][1,3]
  ReMean[i] <- fit1[[i]][2,3]
  fvalue[i] <- fit1[[i]][1,4]
  pvalue[i] <- fit1[[i]][1,5]
  }
 diff <- cbind(probeid,Df,ReDf,SumSq,ReSum,MeanSq,ReMean,fvalue,pvalue)
 diff <- as.data.frame(diff)
 pvalue = as.vector(diff$pvalue)
 adjusted.p.value = p.adjust(pvalue,method = "BH",n=length(pvalue)) #BH或者FDR
 diff <- cbind(diff,adjusted.p.value)
 rownames(diff) <- probeid
 write.table(diff,"Probeid_pvalue.txt",quote = F,sep = "\t",row.names = F)
 
 #probeid to genesymbol
 diff2 <- read.table("Probeid_pvalue.txt",header = T,sep = "\t",row.names = 1)
 diff3 <- diff2[diff2[,"adjusted.p.value"]<0.05,]#选择合适的cutoff值使得差异基因数目合适
 
 diff3=diff3[order(diff3$adjusted.p.value),]
 diff3 <- standard[rownames(diff3),]
 
 match.name <- read.table("probeid_genesymbol.txt", header = T, sep="\t", quote="",row.names=1)
 cross <- intersect(rownames(diff3),rownames(match.name))
 genesymbol <- as.character(match.name[cross,1])
 data_temp <- data.frame(genesymbol,diff2[cross,])
 data <- data.frame(genesymbol,diff3[cross,])
 
 data <- data[!duplicated(data$genesymbol), ] 
 data_temp <- data_temp[!duplicated(data_temp$genesymbol), ] 
 
 
 write.table(data_temp,"mRNA_diff_pvalue.txt",sep = "\t", quote = F, append = F, row.names = F)
 
 write.table(data, "mRNA_diff_expr.txt", sep = "\t", quote = F, append = F, row.names = F)
 #cross <- intersect(rownames(diff3),rownames(annot))
 #genesymbol <- as.character(annot[cross,1])
 #data <- data.frame(genesymbol,diff3[cross,])
 #data <- data[!duplicated(data$genesymbol), ] 
 #rownames(data) <- data[,1]
 #write.table(data,file = "Genesymbol_diff.txt",quote = F,sep = "\t",row.names = F)
 #data = data[,-1]
 
 #normlization
 #matData = betaqn(data) #matdata的热图结果不理想,所以注释掉了
 #data <- as.matrix(data)
 #library("wateRmelon")
 #rlogData <- rlogTransformation(data)
 #rownames(rlogData) <- row.names(data)
 
 pdf("databox.pdf")
 par(mfrow=c(2,1))
 boxplot(data,col = "blue",xaxt = "n",outline = F)
 #boxplot(matData,col = "green",xaxt = "n",outline = F)
 #boxplot(rlogData,col = "red",xaxt = "n",outline = F)
 dev.off()
 
 ###out of R,use log.pl when you need
 
 #分组信息
 diff4 <- read.table("log_mRNA_diff_expr.txt",header = T,sep = "\t",row.names = 1)
 
 data <- diff4
 group_info <- as.data.frame(condition)
 colnames(group_info) <- "group_info"
 rownames(group_info) <- colnames(data)
 ann_colors = list( group_info = c( healthy = "tan2",dcis = "pink2",dic = "lightblue"))
 
 #pdf(file="pheatmap_mat.pdf",width=20,height=10)
 #pheatmap(matData,color=colorRampPalette(c("green","black","red"))(100),fontsize_row=10
  #        ,fontsize_col=3,scale="row",border_color=NA,cluster_col = T,annotation_col = group_info
  #        ,annotation_colors = ann_colors)
 #dev.off()
 
 pdf(file="pheatmap_top_10.pdf",width=20,height=10)
 top10_down_up <- data[c(1:10,(nrow(data)-9):nrow(data)),]
 write.table(top10_down_up,file="top10_down_up.txt",quote=F,sep = "\t",row.names=T,col.names = NA)
 pheatmap(top10_down_up,color=colorRampPalette(c("green","black","red"))(100),fontsize_row=6
          ,fontsize_col=6,scale="row",border_color=NA,cluster_col = F,annotation_col = group_info
          ,annotation_colors = ann_colors)
 dev.off()
 ##matdata的热图结果不理想，选用rlog函数进行标准化绘制热图
 
 