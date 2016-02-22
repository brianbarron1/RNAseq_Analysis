#####
#	Chupp_CellComp_GeneExpression - R script to analyze gene expression of bulk RNAseq Data
#####

#####IDEAS
#

#####Data Preparation

#(1) Generate a data frame containing the expression values (labeled with gene name and gene type)

#Read in the gencode table (for gene name and type), list of RNAseq loci names, and list of sample names
gencode=read.table("/Users/dgrevelation/Bioinformatics/Chupp_CellComp/gencode.v19.dat")
gene_list=read.table("/Users/dgrevelation/Bioinformatics/Chupp_CellComp/GeneIDList.txt", header=TRUE)
sample_names=read.table("/Users/dgrevelation/Bioinformatics/Chupp_CellComp/SE_ChuppData/Sample_Severity.txt")

#Read in the expression data, add the gencode gene name & gene type, and label rows and columns
tpm_data=read.table("/Users/dgrevelation/Bioinformatics/Chupp_CellComp/SE_ChuppData/Expression_Analysis.txt", header=TRUE)
colnames(tpm_data)=sample_names$V1
tpm_data$gene_name=gencode[match(gene_list$gene_id, gencode$V1), 2]
tpm_data$gene_type=gencode[match(gene_list$gene_id, gencode$V1), 8]
rownames(tpm_data)=gene_list$gene_id

#Test
#print(head(tpm_data, n=2000))


###(2) Isolate the relevant gene
#By gene name
gene_expression=subset(tpm_data, gene_name %in% "IL4")
print(gene_expression)

gene_expression=subset(tpm_data, gene_name %in% "IL5")
print(gene_expression)

gene_expression=subset(tpm_data, gene_name %in% "IL13")
print(gene_expression)

gene_expression=subset(tpm_data, gene_name %in% "DEFB1")
print(gene_expression)

# By Ensembl
#tpm_data$Ensembl=gene_list$gene_id

#gene_expression=tpm_data[grep("ENSG00000133048", tpm_data$Ensembl),]
#gene_expression=gene_expression[,1:12]
#colnames(gene_expression)=sample_names$V1
#print(gene_expression)


###### MISC
#(2) Prepare a subset of the data set that corresponds to protein_coding genes
#protein_coding=subset(tpm_data, gene_type %in% "protein_coding")
#rownames(protein_coding)=make.names(protein_coding$gene_name, unique=TRUE)
#protein_coding=protein_coding[,-14]
#protein_coding=protein_coding[,-13]
#colnames(protein_coding)=sample_names$V1

#protein_coding=t(protein_coding)
#print(head(protein_coding))
#print(rownames(protein_coding))

######
# Visualization of Results 
######

#(1) Dendrogram using Euclidian Distance

#Process data by (1) adding small value to avoid -Inf and (2) taking the log of each value 
#tpm_data=tpm_data[,1:6]
#tpm_data=tpm_data[-(1:624),]

#tpm_data=tpm_data+1.0*10^-5
#tpm_data=log(tpm_data)

#tpm_data=t(tpm_data)
#print(head(tpm_data))

#test=dist(tpm_data)
#hc=hclust(test)
#plot(hc, hang=-1, cex=0.6, main="BkSQ Sample Clustering")
