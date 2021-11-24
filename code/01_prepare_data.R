# ----------------------------------------------------------------
# Script name:
# 01_prepare_data.R
#
# Script Description:
# Prepares MiXCR alignment files for clone calling
# Adapted from Silvia Pineda
#
# Notes:
# This script takes in the MiXCR alignment files and generates 
# the following two output files for each tumor type:
# 1) data_for_cloneInfered_Ig_[cancer].txt
# 2) data_BCR_[cancer].rds"
# ----------------------------------------------------------------
library(stringr)
library(data.table)

cancers = c("ACC" , "BLCA", "BRCA", "CESC" ,"CHOL", "COAD" ,"GBM","HNSC" ,"KICH", "KIRC", "KIRP" ,"LGG" ,
            "LIHC" ,"LUAD", "LUSC", "MESO" ,"PAAD" ,"PCPG", "PRAD" ,"READ" ,"SARC", "SKCM" ,"TGCT", "THCA", "THYM",
            "UCEC" ,"UCS" , "UVM")

for (cancer in cancers){
  # Read in the MiXCR alignment files
  files <- list.files(paste("../data/", cancer, "/tmp/", sep = ""))
  files <- files[grep(pattern = "*alignments.txt", files)]
  
  data<-c()
  for(file_name in files) {
    t <- read.delim(paste0("../data/", cancer, "/tmp/",file_name))
    t$sample<-substr(file_name, 1, nchar(file_name)-15)
    t <- t[t$nSeqCDR3 != "",]
    t <- t[!is.na(t$readId), ]
    data <- rbind(data, t)
  }
  
  data$chainType <- ifelse(data$bestVGene!="", substr(data$bestVGene,1,3),ifelse(data$bestVGene=="",substr(data$bestJGene,1,3),NA))
  data$CDR3_length <- nchar(as.character(data$nSeqCDR3))
  data <- data[which(data$CDR3_length>3),] 
  data$seqID <- seq(1,nrow(data))
  data$SEQUENCE_ID <- paste(data$sample,data$seqID,sep="_")
  data$V_J_lengthCDR3 <- paste(data$bestVGene,data$bestJGene,data$CDR3_length,sep="_")
  
  # Prepare data for clone calling
  data_BCR <- data[which(data$chainType=="IGH" |
                                   data$chainType=="IGK" |
                                   data$chainType=="IGL"),]
  
  data_clonesInference_BCR <- data_BCR[,c("SEQUENCE_ID","sample","nSeqCDR3","CDR3_length","bestVGene","bestJGene","V_J_lengthCDR3")]
  write.table(data_clonesInference_BCR,file=paste("../data/", cancer, "/data_for_cloneInfered_Ig_", cancer, ".txt", sep=""),row.names = F,sep="\t")
  saveRDS(data_BCR, paste("../data/", cancer, "/data_BCR_", cancer, ".rds", sep=""))
}
