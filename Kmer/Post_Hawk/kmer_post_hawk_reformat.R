# RScript for Kmer-pval-alignment format

# Command on Command line
# Rscript kmer_post_hawk_reformat.R 'auto' numCores inputdir 
  # Input Directory must only have the pval, abyss, and bwamem output of one comparison
  
# Input Directory (Manual)
  indir<-"~/Desktop/Fish/2020_01_09_sex/"
  case.pval.path<-"pvals_case_top_merged_sorted.txt"
  case.abyss.path<-"case_abyss.25_49.fa"
  case.bwamem.path<-"case_abyss.25_49.fBetSpl5.2.sam"

  control.pval.path<-"pvals_control_top_merged_sorted.txt"
  control.abyss.path<-"control_abyss.25_49.fa"
  control.bwamem.path<-"control_abyss.25_49.fBetSpl5.2.sam"
  
  numCores<-4

# Arguments (Terminal); if you want to run without manually describing path
  args = commandArgs(trailingOnly=TRUE)
  if(args[1]=="auto"){
    numCores<-args[2]
    indir<-args[3]
    setwd(indir)
    files<-list.files()
    
    case.pval.path<-files[grep("sorted", files[grep("pvals_case", files)])]
    case.abyss.path<-files[grep("case_abyss", files)]
    case.bwamem.path<-files[grep("sam", files[grep("case_abyss", files)])]
    
    control.pval.path<-files[grep("sorted", files[grep("pvals_control", files)])]
    control.abyss.path<-files[grep("control_abyss", files)]
    control.bwamem.path<-files[grep("sam", files[grep("control_abyss", files)])]
  }
    
  

# Libraries
  library(foreach)
  library(doParallel)
  library(dplyr)
  library(plyr)

registerDoParallel(numCores)


# Functions:
  # Internal
    rev.comp<-function(x,rev=TRUE){
      x<-toupper(x)
      y<-rep("N",nchar(x))
      xx<-unlist(strsplit(x,NULL))
      for (bbb in 1:nchar(x))
      {
        if(xx[bbb]=="A") y[bbb]<-"T"    
        if(xx[bbb]=="C") y[bbb]<-"G"    
        if(xx[bbb]=="G") y[bbb]<-"C"    
        if(xx[bbb]=="T") y[bbb]<-"A"
      }
      if(rev==FALSE) 
      {
        for(ccc in (1:nchar(x)))
        {
          if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
        }
      }
      if(rev==T)
      {
        zz<-rep(NA,nchar(x))
        for(ccc in (1:nchar(x)))
        {
          zz[ccc]<-y[nchar(x)+1-ccc]
          if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
        }
      }
      return(yy)  
    }
    sam_kmer<-function(data, bwa, rev.idx=NULL){
      
      if(is.na(rev.idx)){
        bwa<-cbind(data, bwa[match(data$abyss_kmer, bwa$V10),])}else{
        bwa<-cbind(data, bwa[match(data$abyss_kmer, rev.idx),])
      }
      bwa<-bwa[which(!is.na(bwa$V1)),]
      bwa$pos<-apply(bwa[,c(1,15)], 1, function(x){
        regexpr(x[1], x[2])[1]})
      bwa$genome_start<-bwa$V4+case_complete_for$pos
      bwa$genome_end<-bwa$V4+bwa$pos+nchar(as.character(bwa$kmerseq))-1
      
      return(bwa)
    }

  # Wrapper
    hawk_data<-function(pval.path, abyss.path, bwamem.path, indir, newdir){
    
    # Import Kmer and Pval
    setwd(indir)
    pvals_case<- read.delim(pval.path, 
                            header=FALSE, stringsAsFactors=FALSE)
    
    case_abyss <- read.csv(abyss.path, 
                           header=FALSE, stringsAsFactors=FALSE)
    
    pvals_case$present<-"no"
    case_abyss$pvals<-""
    case_abyss$idx<-""
    
    # Find matching Cases
    data<-foreach (x=1:nrow(pvals_case), .combine=rbind) %dopar% {
      idx<-grep(pvals_case$V2[x], case_abyss$V1)
      case_abyss$pvals[idx]<-paste0(case_abyss$pvals[idx],";", pvals_case$V1[x])
      case_abyss$idx[idx]<-paste0(case_abyss$idx[idx],";", x)
      if(length(idx)>0){
        cbind(pvals_case$V2[x],  case_abyss$V1[idx], pvals_case$V1[x], x)}
    }
    colnames(data)<-c("kmerseq", "abyss_kmer", "pval", "idx")
    data$length<-nchar(as.character(data$abyss_kmer))
    data$pval<-as.numeric(as.character(data$pval))
    
    # Save Progress Data
    setwd(newdir)
    save(data, file = paste0(ifelse(grepl("case", pval.path), "case", "control"), "_data", ID, ".RData"))
    
    
    # Import Alignment
    bwamem<- read.delim(bwamem.path, 
                             header=FALSE, stringsAsFactors=FALSE)
   
    # Combine Alignment info with Kmer-pval data
    # Forward
    complete_for<-sam_kmer(data, bwamem[which(bwamem$V2 %in% c(0)),])
    
    # Do the same for the complement
    rev.idx<-sapply(bwamem_rev$V10, rev.comp)
    complete_rev<-sam_kmer(data, bwamem[which(bwamem$V2 %in% c(0)),])
    
    # Combine Forward and Reverse
    hawk<-rbind(complete_for, complete_rev)
    hawk$type<-ifelse(grepl("case", pval.path), "case", "control")
    
    return(hawk)
  }


    
# Generate Unique ID
  ID<-sample(1000:9999,1,replace=T)

  newdir<-paste0("Hawk_", Sys.Date(), "_", ID)
  dir.create(newdir)

  
  
# Match Alignment and Pval
  hawk<-hawk_data(case.pval.path, case.abyss.path, case.bwamem.path, indir, newdir)
  hawk2<-hawk_data(control.pval.path, control.abyss.path, control.bwamem.path, indir, newdir)
  

  
# Save Data
  save(hawk, hawk2, file = paste0(newdir, ".RData"))