# ----------------------------------------------------------------
# Script name:
# 03_generate_summary_metrics.R
#
# Script Description:
# Calculates diversity and network metrics for each sample
# Adapted from Silvia Pineda
#
# Notes:
# This script takes in the clone calling and BCR files and 
# generates the following output file for each tumor type:
# 1) [cancer]_FullData.Rdata
# 2) [cancer]_repertoire_diversity.rds
# ----------------------------------------------------------------
library(data.table)
library(stringr)
library(ineq)
library(GenomicDataCommons)

calculate_diversity <- function(data_sample, chain){
  clones_sample <- data_sample[data_sample$chainType == chain, ]
  clones_sample <- clones_sample$V_J_lengthCDR3_CloneId
  clones_proportions <-as.numeric(table(clones_sample))/length(clones_sample)
  top_clone_proportion <- max(clones_proportions)
  shannon_entropy <- -sum(clones_proportions*log2(clones_proportions))
  unique_clones <- length(unique(clones_sample))
  evenness <- shannon_entropy/log2(unique_clones)
  return(list(shannon_entropy, evenness, top_clone_proportion))
}

TCGAtranslateID = function(file_ids, legacy = TRUE) {
  info = files(legacy = legacy) %>%
    filter( ~ file_id %in% file_ids) %>%
    select('cases.samples.submitter_id') %>%
    results_all()
  id_list = lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  barcodes_per_file = sapply(id_list,length)
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    submitter_id = unlist(id_list)))
}

cancers = c("ACC" , "BLCA", "BRCA", "CESC" ,"CHOL", "COAD" ,"GBM","HNSC" ,"KICH", "KIRC", "KIRP" ,"LGG" ,
            "LIHC" ,"LUAD", "LUSC", "MESO" ,"PAAD" ,"PCPG", "PRAD" ,"READ" ,"SARC", "SKCM" ,"TGCT", "THCA", "THYM",
            "UCEC" ,"UCS" , "UVM")

for (cancer in cancers){
  # Read in the files generated in 01_prepare_data.R and 02_clone_calling.py
  data_BCR <- readRDS(paste("../data/", cancer, "/data_BCR_", cancer, ".rds", sep=""))
  nucleotides <-read.csv(paste("../data/", cancer, "/ClonesInfered_Ig_", cancer, ".csv", sep=""))

  data_merge <- merge(data_BCR,nucleotides[,c("SEQUENCE_ID","CloneId")],by=c("SEQUENCE_ID"))
  data_merge$V_J_lengthCDR3_CloneId = paste(data_merge$V_J_lengthCDR3,data_merge$CloneId,sep="_")
  
  # Count reads per chain
  read_count <- table(data_merge$sample)
  read_count_chain <- table(data_merge$sample, data_merge$chainType)
  reads <- data.frame(cbind(read_count,read_count_chain))
  reads$sample <- rownames(reads)
  
  # Normalize reads by total reads per sample
  totalReads = readRDS(paste("../data/total_counts/", cancer, "_total_counts.rds", sep=""))
  totalReads$total_counts = as.numeric(as.character(totalReads$total_counts))
  rownames(totalReads) = totalReads$sample
  id<-match(rownames(reads),totalReads$sample)
  reads$totalReads<-totalReads[id,2]
  
  ##Total reads
  reads$Ig_Reads<-reads$IGH+reads$IGK+reads$IGL

  ####Normalize the number of reads
  reads$IG_expression<-(reads$IGH+reads$IGK+reads$IGL)/reads$totalReads
  reads$IGH_expression<-reads$IGH/reads$totalReads
  reads$IGK_expression<-reads$IGK/reads$totalReads
  reads$IGL_expression<-reads$IGL/reads$totalReads
  
  ###Ratio
  reads$KappaLambda_ratio_expression <- (reads$IGK_expression / reads$IGL_expression)
  
  clones_count<- unique(data_merge[,c("sample","V_J_lengthCDR3_CloneId","chainType")])

  clones<-data.frame(cbind(table(clones_count$sample,clones_count$chainType)))
  colnames(clones)<-c("clones_IGH","clones_IGK","clones_IGL")
  clones <- clones[match(reads$sample, rownames(clones)), ]
  
  ##Diversity measures
  entropy_IGH<-NULL
  entropy_IGK<-NULL
  entropy_IGL<-NULL

  evenness_IGH<-NULL
  evenness_IGK<-NULL
  evenness_IGL<-NULL

  top_clone_prop_IGH<-NULL
  top_clone_prop_IGK<-NULL
  top_clone_prop_IGL<-NULL
  
  sample_ids = reads$sample
  
  for (i in 1:length(sample_ids)){
    data_sample <- data_merge[data_merge$sample == sample_ids[i],]

    diversity_measurements <- calculate_diversity(data_sample, "IGH")
    entropy_IGH[i] = diversity_measurements[[1]]
    evenness_IGH[i] = diversity_measurements[[2]]
    top_clone_prop_IGH[i] = diversity_measurements[[3]]
    
    diversity_measurements <- calculate_diversity(data_sample, "IGK")
    entropy_IGK[i] = diversity_measurements[[1]]
    evenness_IGK[i] = diversity_measurements[[2]]
    top_clone_prop_IGK[i] = diversity_measurements[[3]]
    
    diversity_measurements <- calculate_diversity(data_sample, "IGL")
    entropy_IGL[i] = diversity_measurements[[1]]
    evenness_IGL[i] = diversity_measurements[[2]]
    top_clone_prop_IGL[i] = diversity_measurements[[3]]
    
    print(paste0(cancer, " diversity ", i, "/", length(sample_ids)))
  }
  
  diversity<-cbind(clones,entropy_IGH,entropy_IGK,entropy_IGL,
                   evenness_IGH, evenness_IGK, evenness_IGL,
                   top_clone_prop_IGH, top_clone_prop_IGK, top_clone_prop_IGL)
  
  gini_index <- data.frame()
  for(i in 1:length(sample_ids)){
    gini_index_chains <- NULL
    sample_name <- sample_ids[i]
    sample_data <- data_merge[data_merge$sample == sample_name, ]
    sample_data$V_J_CDR3seq <- paste0(sample_data$V_J_lengthCDR3, "_", sample_data$nSeqCDR3)
    for (chain in c("IGH", "IGK", "IGL")){
      sample_data_chain <- sample_data[sample_data$chainType == chain, ]

      sample_vertex <- data.frame(table(sample_data_chain$V_J_CDR3seq))
      sample_cluster <- table(unique(sample_data_chain[,c("V_J_CDR3seq", "V_J_lengthCDR3_CloneId")])$V_J_lengthCDR3_CloneId)
      
      gini_vertex <- ineq::Gini(sample_vertex$Freq)
      gini_cluster <- ineq::Gini(sample_cluster)
      tmp <- c(gini_vertex, gini_cluster)
      names(tmp) <- c(paste0("gini_vertex_", chain), paste0("gini_cluster_", chain))
      gini_index_chains <- c(gini_index_chains, tmp)
    }
    
    gini_index <- rbind(gini_index, data.frame(t(gini_index_chains)))
    print(paste0(cancer, " gini ", i, "/", length(sample_ids)))
  }
  diversity = cbind(diversity,gini_index)
  
  repertoire_diversity = cbind(reads, diversity)
    
  # Annotate samples with TCGA barcodes
  file_ids = repertoire_diversity$sample
  annotation = TCGAtranslateID(file_ids)
  annotation$patient_barcode = substr(annotation$submitter_id,1,12)
  
  annotation$tumor_type = substr(annotation$submitter_id, 14, 15)
  tcga_code = data.table::fread("data/tcga_sample_code.txt", colClasses = "character")
  tcga_code$Definition = gsub("-", "", tcga_code$Definition)
  tcga_code$Definition = gsub(" ", "_", tcga_code$Definition)
  annotation = merge(annotation, tcga_code[,1:2], by.x = "tumor_type", by.y = "Code")
  annotation$tumor_type = NULL
  
  repertoire_diversity = merge(repertoire_diversity, annotation, by.x = "sample", by.y = "file_id")

  save(data_merge,repertoire_diversity,file=paste("../data/", cancer, "_FullData.Rdata", sep = ""))
  saveRDS(repertoire_diversity, file=paste("../data/", cancer, "_repertoire_diversity.rds", sep = ""))
  print(paste0(cancer, " is done!"))s
  
}


