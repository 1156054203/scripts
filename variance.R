setwd('C:\\Users\\cloudhealth\\Desktop\\1123')
a=read.table("log_norm.txt",header=T,sep="\t",row.names=1) 
a
#预生成2个长度与输入文件行数相同的全为0的向量，将用于存储p value和差异倍数（log2FC）
Pvalue<-c(rep(0,nrow(a))) 
log2_FC<-c(rep(0,nrow(a))) 
# 2~4列是处理组1,5~7列是处理组2；
#将使用循环对每一行进行t检验
#如果某一行两组的标准差都等于0，将无法进行t检验，所以这些行使用NA替代
#每一行计算得到p value和log2FC将被加入原文件的后两列；
#计算log2FC时，每组均值加0.001，是为了防止分母为0导致bug；
for(i in 1:nrow(a)){
  if(sd(a[i,1:9])==0&&sd(a[i,10:18])==0){
    Pvalue[i] <-"NA"
    log2_FC[i]<-"NA"
  }else{
    y=var.test(as.numeric(a[i,1:9]),as.numeric(a[i,10:18]))
    Pvalue[i]<-y$p.value
    log2_FC[i]<-log2((mean(as.numeric(a[i,1:9]))+0.00001)/(mean(as.numeric(a[i,10:18]))+0.00001)) 
  }
}

##############################################
############分布检验##########################
##############################################
for(i in 1:nrow(a)){
  if(sd(a[i,1:9])==0&&sd(a[i,10:18])==0){
    Pvalue[i] <-"NA"
    log2_FC[i]<-"NA"
  }else{
    b=matrix(c(as.numeric(a[i,1:9]),as.numeric(a[i,10:18])),ncol=2)
    x=mvnorm.etest(b,R=1000)
    Pvalue[i]<-x$p.value
    log2_FC[i]<-log2((mean(as.numeric(a[i,1:9]))+0.00001)/(mean(as.numeric(a[i,10:18]))+0.00001)) 
  }
}

# 对p value进行FDR校正
fdr=p.adjust(Pvalue, "BH") 
# 在原文件后面加入log2FC，p value和FDR,共3列；
out<-cbind(a[,0],log2_FC,Pvalue,fdr) 
genesymbol<-rownames(out)
out<-cbind(genesymbol,out)
write.table(out,file="energy-test.xls",quote=FALSE,sep="\t",row.names=FALSE)


setwd('C:\\Users\\cloudhealth\\Desktop\\1122')
a=read.table("log_norm.txt",header=T,sep="\t",row.names=1) 

#单因素多元方差分析,两个前提假设，一个是多元正态性，一个是方差—协方差矩阵同质性
type<-factor(c(rep('c',6),rep('m',6),rep('n',6)))
#预生成2个长度与输入文件行数相同的全为0的向量，将用于存储p value和差异倍数（log2FC）
Pvalue<-c(rep(0,nrow(a))) 
log2_FC<-c(rep(0,nrow(a))) 
# 2~4列是处理组1,5~7列是处理组2；
#将使用循环对每一行进行方差检验
#两组表达量都是0的基因，不检验；
#每一行计算得到p value和log2FC将被加入原文件的后两列；
#计算log2FC时，每组均值加0.001，是为了防止分母为0导致bug；
for(i in 1:nrow(a)){
  if(sum(a[i,1:6])==0&&sum(a[i,7:12])==0&&sum(a[i,13:18])==0){
    Pvalue[i] <-"NA"
    log2_FC[i]<-"NA"
  }else{
    y=aov(as.numeric(a[i,1:18])~type)
    Pvalue[i]<-summary(y)[[1]][,5][1]
    log2_FC[i]<-log2((mean(as.numeric(a[i,1:9]))+0.00001)/(mean(as.numeric(a[i,10:18]))+0.00001)) 
  }
}
# 对p value进行FDR校正
fdr=p.adjust(Pvalue, "BH") 
# 在原文件后面加入log2FC，p value和FDR,共3列；
out<-cbind(a[,0],log2_FC,Pvalue,fdr) 
genesymbol<-rownames(out)
out<-cbind(genesymbol,out)
write.table(out,file="aov.xls",quote=FALSE,sep="\t",row.names=FALSE)
