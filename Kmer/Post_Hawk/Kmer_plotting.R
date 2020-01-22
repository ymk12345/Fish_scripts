####- KMER Analysis and Plot -####

indir<-"/Users/labby/Desktop/Fish"
chr_lengthfile.path<-"~/Desktop/Fish/Hudson_2020-01-07_color/chr_length_fBetSpl5.2.txt"
Rdatafromkmer_post_hawk_reformat.R.path<-"sex_1_20_20.RData"


# Load Data:
  setwd(indir)
  load(Rdatafromkmer_post_hawk_reformat.R.path)

# Import Chromosome Length File  (Chrom, pos)
  chr_lengths <- read.delim(chr_lengthfile.path, header=FALSE, stringsAsFactors=FALSE)
  chr_lengths$wg_pos<-c(0,cumsum(as.numeric(chr_lengths$V2))[-nrow(chr_lengths)])

  chr_lengths$wg_pos_label<-NA
  for(x in 1:nrow(chr_lengths)){
    chr_lengths$wg_pos_label[x]<-chr_lengths$wg_pos[x]+(as.numeric(as.character(chr_lengths$V2[x]))/2)} 

  
#a. Reformat RData for plotting
  hawk_refmt<-function(hawk, category, chr_lengths, color=c("darkblue", "lightblue")){
    
    #a. Chromosome Name and Add Category info
      hawk$chr<-hawk$V3
      hawk$type<-category
  
    #b.  Add a column that creates a genome coordinate for the SNP
      hawk$chr<-as.character(hawk$chr)
      hawk$marker<-hawk$genome_start+chr_lengths$wg_pos[match(hawk$chr, chr_lengths$V1)]
      hawk$marker2<-hawk$genome_end+chr_lengths$wg_pos[match(hawk$chr, chr_lengths$V1)]
  
  
    #c. Relevel the Chromosome and contigs to the appropriate order (critical to be in the right order!)
      ##NOTE: this is for the WILD Betta Genome
        hawk$CHROM<-factor(hawk$chr, 
                         levels = c(unique(hawk$chr[!grepl("CAAAFW0200", hawk$chr)][order(as.numeric(hawk$chr[!grepl("CAAAFW0200", hawk$chr)]))]),
                                    unique(hawk$chr[grepl("CAAAFW0200", hawk$chr)][order(as.numeric(gsub("CAAAFW0200", "", hawk$chr[grepl("CAAAFW0200", hawk$chr)])))])))
  
      #i) Convert the names of the chromosome/contigs to numeric
          hawk$chr<-as.numeric(hawk$CHROM)
          hawk<-hawk[order(hawk$chr),]
  
  
    #d. Transform the p-values to the -log10 scale.
        hawk <- transform(hawk,pval = -log10(pval))
  
  
    #e. Add column "odd.chr" to the table, and find the positions of the chromosomes along the x-axis.
        hawk <- transform(hawk,odd.chr = (chr %% 2) == 1)
        hawk$odd.chr[which(hawk$odd.chr==TRUE)]<-1
        hawk$odd.chr[which(hawk$odd.chr==FALSE)]<-2
  
  
    # f. Assign color depending on the chromosome number (interchanges between dark red and red)
        hawk$color<-ifelse(hawk$chr%%2, color[1], color[2])
        
    return(hawk)
  }
  
  case<-hawk_refmt(hawk, "case", chr_lengths, c("darkred", "lightpink"))
  control<-hawk_refmt(hawk2, "control", chr_lengths)


# b. Generate Plots using output from gemma_plot_df
  #NOTE: REQUIRES concatenated whole genome start and end coordinates!!
  # Whole Genome: wgstart = 0; wgend = end of Chromosomes, not contigs
  wgstart<-0
  wgend<-435314794 


pdf("sex_1_20_20.pdf", height = 10, width = 20)
par(mfrow=c(2,1))
plot(case$marker,case$pval, pch=20, 
     xlim = c(wgstart,wgend), 
     ylim = c(2, 3), ylab= "-log10(pval)",
     col = case$color,
     cex.axis = 1.5, xlab = "", 
     xaxt='n', 
     #yaxs="i",
     xaxs = "i",
     cex.lab = 1.9,
     cex.axis=1.8)
title(main=unique(case$type), cex.main=1)
abline(v=chr_lengths$wg_pos[-grep("CAAAFW0200", chr_lengths$V1)]-1, col = "orange")
axis(1, at =chr_lengths$wg_pos_label[-grep("CAAAFW0200", chr_lengths$V1)]-1, labels = chr_lengths$V1[-grep("CAAAFW0200", chr_lengths$V1)], cex.axis = 1.8)


plot(control$marker,control$pval, pch=20, 
     xlim = c(wgstart,wgend), 
     ylim = c(2,3), ylab= "-log10(pval)",
     col = control$color,
     cex.axis = 1.5, xlab = "", 
     xaxt='n', 
     #yaxs="i",
     xaxs = "i",
     cex.lab = 1.9,
     cex.axis=1.8)
title(main=unique(control$type), cex.main=1)
abline(v=chr_lengths$wg_pos[-grep("CAAAFW0200", chr_lengths$V1)]-1, col = "orange")
axis(1, at =chr_lengths$wg_pos_label[-grep("CAAAFW0200", chr_lengths$V1)]-1, labels = chr_lengths$V1[-grep("CAAAFW0200", chr_lengths$V1)], cex.axis = 1.8)

dev.off()






