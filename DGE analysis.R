setwd("C:/PEducation/PRLD/2nd Year Work/Research Data/Brassica Spring 2016/RNA Seq Analysis/3-8-17 fastqc results/kallisto output files/kallisto_output/abundance files")

# abundance1.tsv <- read.delim("abundance_sh_1_3_9c_S23.tsv")
# abundance2.tsv <- read.delim("abundance_sh_1_7_5c_S24.tsv")
# 
# merged_data2 <- rbind(abundance1.tsv, abundance2.tsv)
# 
# #merge files to create 1 dataframe
# merged_data <- merge(abundance1.tsv, abundance2.tsv, by="target_id")
# head(merged_data)
# data.frame(merged_data)
# 
# ?merge
# ?cbind
# class(abundance1.tsv)

files <- list.files()
files

# 1) import all the tsv files
myfiles = lapply(files, read.delim)
length(myfiles) 

# 2) we want a master read count file, use cbind OR rbind 
dim(myfiles[[1]])

tmp <- do.call(cbind, myfiles)
head(tmp)
dim(tmp)

# 3) extract the est.count column (dplyr)
# install.packages("tidyverse")
library("tidyverse")

read.count.colnames <- grep("est_counts", colnames(tmp), value = TRUE)
read.count <- tmp[,c(read.count.colnames)]
head(read.count)

# get cleaner sample ID 
sample.ID <- gsub("(abundance_)([[:print:]]+)(.tsv)", "\\2", files)
sample.ID

# assign sample ID to colnames 
colnames(read.count) <- sample.ID
head(read.count)

# get geneID
gene.ID <- tmp[,"target_id"]

# combine geneID & read count to make read count 
read.count.final <- cbind(gene.ID, read.count)
dim(read.count.final)
View(read.count.final)

#make histogram to view count distribution
library(ggplot2)
hist1 <- qplot(read.count.final$sh_1_3_9c_S23) + scale_x_log10()
hist1

#keep only genes with >10 reads in more than 3 samples
read.count.final2 <- data.matrix(read.count.final)
read.count.final2 <- data.frame(read.count.final2)
read.count.final2 <- read.count.final[rowSums(read.count.final2 > 10) >= 3, ]

#make sample description and set reference to sun treatment
sample.description <- data.frame(
  leaf=regmatches(colnames(read.count.final2), regexpr("leaf3|leaf7|leaf5", colnames(read.count.final2))),
  trt=regmatches(colnames(read.count.final2), regexpr("sh|shsu|sun|sush", colnames(read.count.final2)))
)
sample.description$group <- paste(sample.description$leaf, sample.description$trt,sep="_")
sample.description$trt <- relevel(sample.description$trt,ref="sun")

#calculate normalization factors
library(edgeR)
read.count.final3 <- subset(read.count.final2, select= -c(gene.ID))
dge.data <- DGEList(counts=read.count.final3, group=sample.description$group)
dim(dge.data) 
dge.data <- calcNormFactors(dge.data, method = "TMM")
dge.data$samples

#BCV plot
plotMDS(dge.data, method = "bcv")

#make design matrix
design <- model.matrix(~leaf+trt,data = sample.description)
rownames(design) <- sample.description$sample
design

#calculate dispersion factors
dge.data <- estimateGLMCommonDisp(dge.data,design,verbose = TRUE)
