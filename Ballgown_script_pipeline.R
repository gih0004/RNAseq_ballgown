#You need to have directory named ballgown in which all subfolders with stringtie output will be arranged. Along with .gtf, each subfolder representing each sample will have 5 .ctab files. Please do not touch them i.e. delete or rename them, or else following scripts will not work.



devtools::install_github('alyssafrazee/RSkittleBrewer')


library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(ggplot2)
library(gplots)
library(devtools)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# change directories to the correct one for your own computer 
setwd("C:/Users/guill/Rstudio/OmicsInAG/ballgown/")
list.files()
list.files("ballgown/")


#you want to note this order of files in the ballgown folder that you have donwloaded from HPC. Now in the .csv file, it should be the same order.
#we want our pheno_df matrix to have two columns -- the first column being the samples as they appear in the ballgown directory, and the second being the condition of each sample in the column before it (say resistant vs susceptible).
# where <> exist, those strings must be changed to there descriptors and the <> eliminated for proper running of this script 
sample<-c("<file.gtf>", "<file.gtf>", "<file.gtf>","<file.gtf>" )
<time point distinguisher><-c(rep("<timepoint>", <number of reps>), rep("<timepoint>", <number of reps>), rep("<timepoint>", <number or reps>)) 

#sample1=c("HA_10D_2_ACCTTGGC-GGCCTCAT_S39_L002.sunflower.sam.sort" , "HA_10D_3_ATATCTCG-ATCTTAGT_S40_L002.sunflower.sam.sort", "HA_20D_2_GCGCTCTA-GCTCCGAC_S41_L002.sunflower.sam.sort", "HA_20D_3_GCGCTCTA-GCTCCGAC_S41_L002.sunflower.sam.sort")
#dpf=c(rep("10D", 2), rep("20D", 2))


pheno_df<-data.frame("sample"=sample,"<dpf>"=<dpf>) #timepoint destinguisher for samples. <dpf>=days post flowering, should be changed accordingly to the user's sample timepoints

pheno_df


#check if the table looks correct with correct columns.
# change the directory to your computer
bg <- ballgown(dataDir = "C:/Users/guill/Rstudio/OmicsInAG/ballgown/", pData=pheno_df, samplePattern = "HA")

bg

save(bg, file='bg.rda1')

#this ballgown object contains expr (tables containing expression data for genomic features), structure and indexes. 
#components of bg object: expr(), structure(), indexes()
#meas- types of measurements contained in the data slots ('rcount', 'ucount', 'cov', 'mcov','FPKM')


#A ballgown object has six slots: structure, expr, indexes, dirs, mergedDate, and meas.

#extract sampleNames
sampleNames(bg)

#Exon, intron, and transcript structures are easily extracted from the main ballgown object:
structure(bg)$exon

structure(bg)$intron


structure(bg)$trans

#expr
#The expr slot is a list that contains tables of expression data for the genomic features. These tables are very similar to the *_data.ctab Tablemaker output files. Ballgown implements the following syntax to access components of the expr slot:
  
 # *expr(ballgown_object_name, <EXPRESSION_MEASUREMENT>)

#where * is either e for exon, i for intron, t for transcript, or g for gene, and is an expression-measurement column name from the appropriate .ctab file. Gene-level measurements are calculated by aggregating the transcript-level measurements for that gene. 
#All of the following are valid ways to extract expression data from the bg ballgown object:

# FPKM - Fragments per kilobase of transcript per million mapped reads
# normalized reads to account for the difference in sequencing depth - i.e., more sequencing will have more reads mapping to each gene
# normalizing for the length of the gene - i.e., longer genes will have more reads mapping to them

transcript_fpkm = texpr(bg, 'FPKM') # fragments per kilobase of transcript per million mapped reads
transcript_cov = texpr(bg, 'cov') # average per base read coverage
whole_tx_table = texpr(bg, 'all') # all transcript level expression data with extra metadata
exon_mcov = eexpr(bg, 'mcov') # multi map corrected average per base read coverage
junction_rcount = iexpr(bg) # intron expression
whole_intron_table = iexpr(bg, 'all') #intron level expression data with extra metadata
gene_expression_bg = gexpr(bg) # gene expression in FPKM

exon_transcript_table_bg = indexes(bg)$e2t # exons to transcripts table
transcript_gene_table_bg = indexes(bg)$t2g # transcripts to genes table


write.csv(transcript_gene_table_bg, "transcript_gene_table_bg.csv")

head(transcript_gene_table_bg)

head(gene_expression_bg)
#setwd("C:/Users/guill/Rstudio/OmicsInAG/ballgown/HA_10D_2_ACCTTGGC-GGCCTCAT_S39_L002.sunflower.sam.sort")
#genes=read.table(file="t_data.ctab", sep="\t", header=TRUE)


#Filter to remove low-abundance genes. One common issue with RNA-seq data is that genes often have very few or zero counts. A common step is to filter out some of these. Another approach that has been used for gene expression analysis is to apply a variance filter. Here we remove all transcripts with a variance across samples less than one:
bg_filt = subset(bg, "rowVars(texpr(bg))>1" ,genomesubset=TRUE)
bg_filt

save(bg_filt, file='bg_filt.rda')



#To find transcripts that are deferentially expressed under variation. We use stattest function from Ballgown. We set getFC-TRUE so we that we can look at the confounder adjusted fold change between two groups.
results_transcripts_bg = stattest(bg_filt, feature="transcript", meas = "FPKM" , covariate = "<dpf>", getFC="TRUE" ) #covariate is what your comparing between samples (could be genotype or time points )


#Identify genes that show statistically significant differences between groups. For this we can run the same function that we used to identify differentially expressed transcripts, but here we set feature=”gene” in the stattest command:

results_genes_bg = stattest(bg_filt, feature="gene" , covariate = "<dpf>" , meas = "FPKM" , getFC="TRUE")
results_genes_order_bg<- results_genes_bg[ order (-results_genes_bg$fc),]
write.csv(results_genes_order_bg, "results_gene_order_bg.csv")

#Add gene names and gene IDs to the results_transcripts data frame:
results_transcripts_bg = data.frame(geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt), results_transcripts_bg)


#Sort the results from the smallest P value to the largest:
results_transcripts_bg = arrange(results_transcripts_bg,pval)
results_genes_bg = arrange(results_genes_bg,pval)

#The log fold change (logFC) can be computed.
results_transcripts_bg$logFC <- log2(results_transcripts_bg$fc)
results_genes_bg$logFC <- log2(results_genes_bg$fc)


#Write the results to a csv file that can be shared and distributed: Identify transcripts and genes with a q value <0.05:
write.csv(results_transcripts_bg, "transcript_results_bg.csv", row.names=FALSE)
write.csv(results_genes_bg, "gene_results_bg.csv", row.names=FALSE)


#Finally, in order to keep only the statistically significant genes, it is possible to filter the data using a Q-Value threshold of 0.05.
results_transcripts_bg_q0.05 <- subset(results_transcripts_bg,results_transcripts_bg$qval<0.05)
results_genes_bg_q0.05 <- subset(results_genes_bg,results_genes_bg$qval<0.05)

#Write the results to a csv file that can be shared and distributed: Identify transcripts and genes with a q value <0.05:
write.csv(results_transcripts_bg_q0.05, "transcript_results__bg_q0.05.csv", row.names=FALSE)
write.csv(results_genes_bg_q0.05, "gene_results_bg_q0.05.csv", row.names=FALSE)

#Visualisation-Volcano plot

volcano_plot <- function(data){
  logfc.threshold <- 1
  with(data, plot(logFC, -log10(qval), pch=20, main="Volcano plot"))
  with(subset(data, qval<.05 ), points(logFC, -log10(qval), pch=20, col="red"))
  with(subset(data, abs(logFC)>logfc.threshold), points(logFC, -log10(qval), pch=20, col="orange"))
  with(subset(data, qval<.05 & abs(logFC)>logfc.threshold), points(logFC, -log10(qval), pch=20, col="green"))
}

volcano_plot(results_genes_bg)





#now with ggplot because that looks good. 
library(ggrepel)

logfc.threshold <- 5
ggplot() + 
  geom_point(data = subset(results_genes_bg, qval>.05 & abs(logFC)<logfc.threshold), aes(x = logFC, y = -log10(qval)), fill = "grey", shape = 21, alpha = 0.2) + 
  geom_point(data = subset(results_genes_bg, qval<0.05), aes(x = logFC, y = -log10(qval)), fill = cbbPalette[[7]], shape = 21) + 
  geom_point(data = subset(results_genes_bg, abs(logFC)>logfc.threshold), aes(x = logFC, y = -log10(qval)), fill = cbbPalette[[5]], shape = 21) + 
  geom_point(data = subset(results_genes_bg, qval<.05 & abs(logFC)>logfc.threshold), aes(x = logFC, y = -log10(qval)), fill = cbbPalette[[3]], shape = 21) +
  theme_classic() +
  geom_text_repel(data = subset(results_genes_bg, -log10(qval)>10 & abs(logFC)>2), aes(x = logFC, y = -log10(qval), label = id), size = 2)
  
  
#Distribution

tropical= c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
fpkm = texpr(bg_filt,meas="FPKM")
fpkm = log2(fpkm+1)
par(mar = c(10,5,4,2) + 0.1)
boxplot(fpkm,col=as.numeric(pheno_df$<dpf>),las=2,ylab='log2(FPKM+1)')
#The boxplot gives an overview of expression of fpkm values of different genes across different samples. We have log10 transformed the values to visualise it better and added 1 gene_expression+1 to avoid errors if the fpkm values are 0.

# do it with ggplot
library(tidyr)
fpkm2 <- data.frame(fpkm)
fpkm2$id <- rownames(fpkm2)
pivot_fpkm2 <- pivot_longer(fpkm2, cols = 1:length(colnames(fpkm2))-1, names_to = "sample", values_to = "fpkm")

pivot_fpkm3 <- separate(pivot_fpkm2, col = 2, sep = "_", into = c("sample", "genotype", "time", "rep"))
pivot_fpkm3$time <- factor(pivot_fpkm3$time, levels = c("10D", "20D", "35D"))
ggplot(pivot_fpkm3, aes(x = time, y = fpkm)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_manual(values = cbbPalette) + 
  ylab("log2 FPKM") + 
  xlab("")

bg_filt@structure


#Data variation for one gene
gene.data.frame <- data.frame(fpkm[11,]); colnames(gene.data.frame) <- "fpkm"
gene.data.frame$sample <- rownames(gene.data.frame)
gene.data.frame <- separate(gene.data.frame, col = 2, sep = "_", into = c("sample", "<dpf>", "time", "rep"))
gene.data.frame$time <- factor(gene.data.frame$time, levels = c("10D", "20D"))

ggplot(gene.data.frame, aes(x = <dpf>, y = fpkm, fill = <dpf>)) + 
  geom_boxplot() + 
  geom_jitter() +
  #facet_wrap(~genotype) + 
  theme_classic() +
  scale_fill_manual(values = cbbPalette) + 
  ylab("log2 FPKM") + 
  xlab("")
  



plot(fpkm[11,] ~ pheno_df$<dpf>, border=c(1,2), main=paste(ballgown::geneNames(bg_filt)[11],' : ', ballgown::transcriptNames(bg_filt)[11]),pch=19, xlab="<dpf>", ylab='log2(FPKM+1)')
points(fpkm[11,] ~ jitter(as.numeric(pheno_df$<dpf>)), col=as.numeric(pheno_df$<dpf>))

#Plot transcripts
sampleNames(bg_filt)

plotTranscripts(ballgown::geneIDs(bg_filt)[ballgown::geneIDs(bg_filt) == "MSTRG.13012"], bg_filt, main=c('Gene CLV3'), sample=c('HA_10D_2_ACCTTGGC-GGCCTCAT_S39_L002.sunflower.sam.sort', 'HA_10D_3_ATATCTCG-ATCTTAGT_S40_L002.sunflower.sam.sort', 'HA_35D_2_TGGTGGCA-TCCTGTAA_S45_L002.sunflower.sam.sort', 'HA_35D_3_AGGCAGAG-AGAATGCC_S46_L002.sunflower.sam.sort'), labelTranscripts = TRUE)


plotTranscripts(ballgown::geneIDs(bg_filt)[ballgown::geneIDs(bg_filt) == "MSTRG.13012"], bg_filt, main=c('CLV3 10D_2'), sample=c('HA_10D_2_ACCTTGGC-GGCCTCAT_S39_L002.sunflower.sam.sort'), legend = TRUE, labelTranscripts = TRUE)
plotTranscripts(ballgown::geneIDs(bg_filt)[ballgown::geneIDs(bg_filt) == "MSTRG.13012"], bg_filt, main=c('CLV3 35D_2'), sample=c('HA_35D_2_TGGTGGCA-TCCTGTAA_S45_L002.sunflower.sam.sort'), legend = TRUE, labelTranscripts = TRUE)

plotTranscripts(ballgown::geneIDs(bg_filt)[ballgown::geneIDs(bg_filt) == "MSTRG.13012"], bg_filt, main=c('CLV3 10D_3'), sample=c('HA_10D_3_ATATCTCG-ATCTTAGT_S40_L002.sunflower.sam.sort'), legend = TRUE, labelTranscripts = TRUE)
plotTranscripts(ballgown::geneIDs(bg_filt)[ballgown::geneIDs(bg_filt) == "MSTRG.13012"], bg_filt, main=c('CLV3 35D_3'), sample=c('HA_35D_3_AGGCAGAG-AGAATGCC_S46_L002.sunflower.sam.sort'), legend = TRUE, labelTranscripts = TRUE)


# different isoforms of differntially expressed genes
plotMeans('MSTRG.13012', bg_filt, groupvar ="<dpf>", legend=TRUE) 


#sampleNames(bg_filt)

#levels <- c(1,0)
#pData(bg_filt) <- data.frame(sampID=sampleNames(bg_filt), group=rep(levels, each=9))
#pData(bg_filt)

##plotMeans(MSTRG.23241, bg_filt, groupvar='<genotype', meas='FPKM', colorby='transcript')

#Heatmap
library(RColorBrewer)
library(gplots)
my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)
mat_data <- data.matrix(fpkm)
heatmap.2(mat_data,
          main = "FPKM",
          notecol="black",
          # density.info="none",
          trace="none",
          margins =c(12,9),
          col=my_palette,
          # breaks=col_breaks,
          dendrogram="both")


#plotly volcano plots
library("ggplot2")
library("gridExtra")
library("plotly")

input_file <- ("C:/Users/guill/Rstudio/OmicsInAG/ballgown/transcript_results_bg.csv")

#read the file into a dataframe
diff_df_bg <- read.delim(file=input_file,header=TRUE,sep=',')

#check some attributes of the data
colnames(diff_df_bg)

dim(diff_df_bg)

# keep only the fields needed for the plot
# FDR = false discovery rate = adjusted p value = significance 
diff_df_bg <- diff_df_bg[c("geneIDs", "logFC", "qval")]

# preview the dataset; data required for the plot
head(diff_df_bg)

# add a grouping column; default value is "not significant"
diff_df_bg["group"] <- "NotSignificant"

# for our plot, we want to highlight 
# FDR < 0.05 (significance level)
# Fold Change > 1.5

# change the grouping for the entries with significance but not a large enough Fold change
diff_df_bg[which(diff_df_bg['qval'] < 0.05 & abs(diff_df_bg['logFC']) < 1.5 ),"group"] <- "Significant"

# change the grouping for the entries a large enough Fold change but not a low enough p value
diff_df_bg[which(diff_df_bg['qval'] > 0.05 & abs(diff_df_bg['logFC']) > 1.5 ),"group"] <- "FoldChange"

# change the grouping for the entries with both significance and large enough fold change
diff_df_bg[which(diff_df_bg['qval'] < 0.05 & abs(diff_df_bg['logFC']) > 1.5 ),"group"] <- "Significant&FoldChange"


# Find and label the top peaks..
top_peaks_bg <- diff_df_bg[with(diff_df_bg, order(logFC, qval)),][1:8,]
top_peaks_bg <- rbind(top_peaks_bg, diff_df_bg[with(diff_df_bg, order(-logFC, qval)),][1:8,])


# Add gene labels to the plot
# Single Gene Annotation example
# m <- diff_df[with(diff_df, order(fc, pval)),][1,]
# a <- list(
#   x = m[["fc"]],
#   y = -log10(m[["pval"]]),
#   text = m[["id"]],
#   xref = "x",
#   yref = "y",
#   showarrow = TRUE,
#   arrowhead = 7,
#   ax = 20,
#   ay = -40
# )

# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly
a <- list()
for (i in seq_len(nrow(top_peaks_bg))) {
  m <- top_peaks_bg[i, ]
  a[[i]] <- list(
    x = m[["logFC"]],
    y = -log10(m[["qval"]]),
    text = m[["geneIDs"]],
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 0.5,
    ax = 20,
    ay = -40
  )
}

diff_df_bg
save(diff_df_bg, file='diff_df10to35')
# make the Plot.ly plot
p <- plot_ly(data = diff_df_bg, x = ~logFC, y = ~-log10(qval), text = geneIDs, color = ~group) %>% add_trace(type = "scatter", mode = "markers")%>%
  layout(title ="Volcano Plot 10 to 35 D") %>%
  layout(annotations = a)
p

save(p, file='diffexpplot 20 to 35 D ')
# significant and fold change=== positive logFC MSTRG.6223 (Pepper1.55ch03	StringTie	transcript	25701551	25707456	.	+	.	gene_id "MSTRG.6223", CA03g08350),
#19659, 20079, 25377,10483; negative == 6323, 5189, 20079, 4017

