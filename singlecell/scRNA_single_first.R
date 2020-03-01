# https://www.yuque.com/docs/share/4af2f15f-418c-44c7-bc24-526a9afa1591#

source("/Data2/wanglp/bin/Seurat/Self/library/scRNASeq_lib.R")
#----------------------------------------------
# 0.input
#----------------------------------------------
config <- data.frame(testMethod='MAST',PCs=50,clusRes=0.8,species='mouse')
##species: human(Hs), mouse(Mm), rat(Rn), chicken(Gg), pig(Ss), chimp(Pt), fly(Dm), worm(Ce), zebrafish(Dr)

#source("/Data2/wanglp/bin/scRNASeq/Seurat/scRNA_library_wlp.R")
pathHome="/Data/wanglp/project/JYSH1903LHT27_Zhengyx_2scATAC_20190328/scRNASeq/seurat/single"
sampleInfo="/Data/wanglp/project/JYSH1903LHT27_Zhengyx_2scATAC_20190328/scRNASeq/seurat/sampleSheet2.csv"
markerFile <- "/Data/wanglp/project/JYSH1903LHT27_Zhengyx_2scATAC_20190328/scRNASeq/seurat/marker_list.csv"

#----------------------------------------------
# 1.read data
#----------------------------------------------
sampleSheet <- read.csv(file=sampleInfo, header=T,colClasses="character")
sampleSheet$sampleName <- as.character(sampleSheet$sampleName)
sampleSheet$sampleData <- as.character(sampleSheet$sampleData)

# get marker list
markerList <- read.csv(file=markerFile, header=T)
sampleNum <- length(sampleSheet$sampleName)
config$sampleCounter=sampleNum

#----------------------------------------------
## 2.start basis analysis
#----------------------------------------------

for (Scount in 1:config$sampleCounter) {
  sample=sampleSheet$sampleName[Scount]
  cat("===================================================","\n")
  cat("start do analysis for --- ",sample,"at ",date(),"\n")
  if(!file.exists(sample)) {dir.create(paste(pathHome,"/",sample,sep="")}
  pathSample = paste(pathHome,sample,sep = "/")
  seuratData <- Read10X(data.dir = sampleSheet$sampleData[Scount])
  seuratObj <- CreateSeuratObject(counts = seuratData, project = sample, min.cells = 3, min.features = 200)

  ## step 001
  step1<-paste(pathSample,'/001.Basic_Analysis',sep="")
  if(!file.exists(step1)) {dir.create(step1)} 
  setwd(step1)
  seuratObj <- getInitialPNG(seuratObj)
  nameObj = paste(pathSample,"/",sample,"_orignal.rds", sep = "")
  saveRDS(seuratObj, file = nameObj)
  # seuratObj <- readRDS(file=nameObj)
  
  ## getFiltObj <- function(seuratObj, LowThres.nUMI=0, HighThres.nUMI=Inf, LowThres.nGene=200, HighThres.nGene=Inf, LowThres.MitoPercent=0,HighThres.MitoPercent=0.1)
  seuratObj <- getFiltObj(seuratObj)
  cat("start do analysis for --- ",sample,"at ",date(),"\n")
}