#####
# Chupp_RNAseqData_ComparetoMicro.R - R Script to compare RNAseq Data with Microarray Data
#####

### New_Ideas ###
# +Isolate Cool Genes

### Variables ###
#Sample Name for Comparison.
samp_name="Sample_000_554"
#Column number for RNAseq data and Microarray data.
RnaCol=6
MicroCol=72
#Correlation Value.
Correl="r=0.618"

### Data Preparation ###

##(1) Read in data files and generate data frames
#Read in the gencode table (for gene name and type), list of RNAseq loci names, and list of sample names
gencode=read.table("/home/brian/Bioinformatics/Chupp_CellComp/gencode.v19.dat")
gene_list=read.table("/home/brian/Bioinformatics/Chupp_CellComp/GeneIDList.txt", header=TRUE)
sample_names=read.table("/home/brian/Bioinformatics/Chupp_CellComp/SE_ChuppData/Sample_Names.txt")

#Read in the RNA-seq expression data. Add the gencode name & gene type as separate columns. Label rows and columns. Remove early reads that do not correspond to genes.
tpm_data=read.table("/home/brian/Bioinformatics/Chupp_CellComp/SE_ChuppData/Expression_Analysis.txt", header=TRUE)
colnames(tpm_data)=sample_names$V1
tpm_data$gene_name=gencode[match(gene_list$gene_id, gencode$V1), 2]
tpm_data$gene_type=gencode[match(gene_list$gene_id, gencode$V1), 8]
rownames(tpm_data)=gene_list$gene_id
tpm_data=tpm_data[-(1:624),]

#Read in the microarray data (col -> sample names and rows -> gene names with values). Match for gene name.
microarray=read.table("/home/brian/Bioinformatics/Chupp_CellComp/sputum_protein_coding_gene_expression.txt", header=TRUE)
microarray$gene_name=rownames(microarray)

##(2) Align data for comparison. Create a new data frame with values for both RNAseq and microarray data (no expression recorded as "NA") 
#Match (by gene name) the expression data for the RNAseq data with the microarray data (using match) 
RNAseq=tpm_data[, RnaCol]
Micro=microarray[match(tpm_data$gene_name, microarray$gene_name), MicroCol]

compare_data=data.frame(RNAseq, Micro, tpm_data$gene_name)
compare_data$RNAseq[compare_data$RNAseq==0.00] <- NA

### Data Interpretation

#(1) Determine which genes have no detectable expression in the RNAseq data or Microarray data.
RNAseq_noexpression_data <- compare_data[is.na(compare_data$RNAseq),]
cat(as.character(RNAseq_noexpression_data$tpm_data.gene_name), file=paste("/home/brian/Bioinformatics/Chupp_CellComp/MicroArray_Comparison/", samp_name, "_RnaseqNoExpress_GeneNames.txt", sep=""), sep="\n")

Microarray_noexpression_data <- compare_data[is.na(compare_data$Micro),]
cat(as.character(Microarray_noexpression_data$tpm_data.gene_name), file=paste("/home/brian/Bioinformatics/Chupp_CellComp/MicroArray_Comparison/", samp_name, "_MicroarrayNoExpress_GeneNames.txt", sep=""), sep="\n")

#(2) Count the number of genes with expression in either data set.
micro_length=length(na.omit(Micro))
RNAseq_length=length(na.omit(RNAseq))
cat("The total number of loci in the Microarray data is", micro_length, "\n")
cat("The total number of loci in the RNAseq data is ", RNAseq_length, "\n")

RNAseq_df=as.data.frame(RNAseq)
RNAseq_df$RNAseq[RNAseq_df$RNAseq==0.00] <- NA
RNAseq_expression=na.omit(RNAseq_df)
RNAseq_noexpression <- RNAseq_df[rowSums(is.na(RNAseq_df)) > 0,] 

diff=(length(RNAseq) - length(RNAseq_expression$RNAseq) -(length(RNAseq)-length(na.omit(RNAseq))))
cat("There are", diff, "loci with no expression in the RNAseq data \n")

#(3) Create the final dataframe with both RNAseq and microarray data for plotting.
compare_data=compare_data[,-3]
compare_data=na.omit(compare_data)

compare_data=data.frame(apply(compare_data[1],2,log), compare_data[2])
colnames(compare_data)=c("RNAseq","Microarray")

### Data Visualization ###

#(1) Generate the plot (using fields and smoothscatter) with data points from compare_data and linear regression line.
library("fields")

fudgeit <- function(){
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F)
}

par(mar = c(5,4,4,5) + .1)
smoothScatter(compare_data$Microarray, compare_data$RNAseq, nbin=64, colramp=colorRampPalette(c("white","blue","cyan","green","yellow","red","brown")), nrpoints=0, postPlotHook = fudgeit, xlab="Microarray Intensity", ylab="RNAseq Intensity", main=samp_name)

reg1=lm(compare_data$RNAseq~compare_data$Microarray)
abline(reg1)

legend(2, 11, Correl)

#(2) Calculate the correlation coefficient.
print(cor(compare_data$Microarray, compare_data$RNAseq, use="complete"))
