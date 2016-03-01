#####
#	Microarray_Correlation.Rd - Measures the correlation between RNA-seq and microarray data
#####

###Data Preparation

#Read in the gencode table (for gene name and type), list of RNAseq loci names, and list of sample names
gencode=read.table("/Users/dgrevelation/Bioinformatics/Chupp_CellComp/gencode.v19.dat")
gene_list=read.table("/Users/dgrevelation/Bioinformatics/Chupp_CellComp/GeneIDList.txt", header=TRUE)
sample_names=read.table("/Users/dgrevelation/Bioinformatics/Chupp_CellComp/SE_ChuppData/Sample_Names.txt")

#Read in the RNA-seq expression data. Add the gencode name & gene type as separate columns. Label rows and columns.
tpm_data=read.table("/Users/dgrevelation/Bioinformatics/Chupp_CellComp/SE_ChuppData/Expression_Analysis.txt", header=TRUE)
colnames(tpm_data)=sample_names$V1
tpm_data$gene_name=gencode[match(gene_list$gene_id, gencode$V1), 2]
tpm_data$gene_type=gencode[match(gene_list$gene_id, gencode$V1), 8]
rownames(tpm_data)=gene_list$gene_id

#Read in the microarray data (col -> sample names and rows -> gene names with values). Match for gene name.
microarray=read.table("/Users/dgrevelation/Bioinformatics/Chupp_CellComp/sputum_protein_coding_gene_expression.txt", header=TRUE)
microarray$gene_name=rownames(microarray)

#
###Scatterplot Visualization
#

#Extract the relevant expression data from each dataframe. 
Micro=microarray[,19]
RNAseq=tpm_data[match(microarray$gene_name, tpm_data$gene_name), 2]

#Floor the RNAseq data (assign 0 values to 0.01325)
#RNAseq=RNAseq*100
#RNAseq[RNAseq < 0.01325] <- 0.01325
RNAseq=RNAseq+1.0*10^-4
RNAseq=log(RNAseq)

reg1=lm(RNAseq~Micro)
plot(Micro, RNAseq, main="Sample_000_554")
print(cor(Micro, RNAseq, use="complete"))
abline(reg1)

###For comparison to other samples (to make sure if it's specific)
n=sample(1:113, 1)
Micro=microarray[,n]
RNAseq=tpm_data[match(microarray$gene_name, tpm_data$gene_name), 2]

#Floor the RNAseq data (assign 0 values to 0.01325)
#RNAseq=RNAseq*100
#RNAseq[RNAseq < 0.01325] <- 0.01325
RNAseq=RNAseq+1.0*10^-4
RNAseq=log(RNAseq)

print(cor(Micro, RNAseq, use="complete"))

###For comparison to other samples (to make sure if it's specific)
n=sample(1:113, 1)
Micro=microarray[,n]
RNAseq=tpm_data[match(microarray$gene_name, tpm_data$gene_name), 2]

#Floor the RNAseq data (assign 0 values to 0.01325)
#RNAseq=RNAseq*100
#RNAseq[RNAseq < 0.01325] <- 0.01325
RNAseq=RNAseq+1.0*10^-4
RNAseq=log(RNAseq)

print(cor(Micro, RNAseq, use="complete"))

###For comparison to other samples (to make sure if it's specific)
n=sample(1:113, 1)
Micro=microarray[,n]
RNAseq=tpm_data[match(microarray$gene_name, tpm_data$gene_name), 2]

#Floor the RNAseq data (assign 0 values to 0.01325)
#RNAseq=RNAseq*100
#RNAseq[RNAseq < 0.01325] <- 0.01325
RNAseq=RNAseq+1.0*10^-4
RNAseq=log(RNAseq)

print(cor(Micro, RNAseq, use="complete"))

#####
#df=data.frame(Micro, RNAseq)
#lm.out = lm(Micro ~ RNAseq + I(RNAseq^2) + I(RNAseq^3), data=df)
#print(summary(lm.out))
