#source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
library("vsn")
library(biomaRt)
library(GenomicFeatures)

# defining a mart-object containing information about the mm9 genome
mart = useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl', host="may2012.archive.ensembl.org")

# read the table with information from all RNA-Seq experiments 
setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/mapping/")
data <- read.table('fusion.txt', header = T, fill=T, stringsAsFactors = F)
data <- data[,c(-5:-9)]
colnames(data)<- sub("X","",colnames(data))

# prepare data for analysis with DESeq2 by creating a matrix from the data frame 
# and creating a coldata frame with information about the conditions and types of the experiments
countdata <- as.matrix(data)
condition <- factor(c(rep("48h",4), rep("ESC",2)))
type <- factor(c(rep("paired-end",1),c(rep("single-read",4), rep("paired-end",1))))
coldata <- data.frame(condition,type)

# calculating the DESeq-functions
dds <- DESeqDataSetFromMatrix(countdata, colData = coldata, design = ~ condition + type)
dds <- DESeq(dds)
norm <- sizeFactors(dds)
norm
sd(norm)


#Varianz mit log-transformierten Werten
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
#Varianz mit Varianz-stabilisierten Werten
vsd <- vst(dds)
meanSdPlot(assay(vsd))



# extracting results of the dds-data by only looking at the contrast between the both stages
# and with the p-value adjustment method of Benjamini-Hochberg, by testing for a significance level of 0.01
r <- results(dds, alpha = 0.01,pAdjustMethod="BH",contrast = c("condition","ESC","48h"))
summary(r)



#VOLCANO-PLOT 
cols <- densCols(r$log2FoldChange, -log10(r$pvalue))
plot(r$log2FoldChange, -log10(r$padj), panel.first=grid(), col=cols,
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(0.01), col="brown")



# differentially expressed genes
genes.selected1 <- abs(r$log2FoldChange) > 5 & r$padj < 0.01
sum(genes.selected1,na.rm=TRUE)

# differentially expressed genes:: up-regulated
genes.selected2 <- r$log2FoldChange > 5 & r$padj < 0.01
sum(genes.selected2,na.rm=TRUE)
#[1] 370: signifikant hochregulierte

# differentially expressed genes: down-regulated
genes.selected3 <- r$log2FoldChange < (-5) & r$padj < 0.01
sum(genes.selected3,na.rm=TRUE)


# annotation of ensembl-geneIDs to gene names and identify most significant genes
r.ordered <- r[order(r$padj,decreasing = F),]
un <- unlist(lapply(rownames(head(r.ordered,n=50)), function(x) getBM(attributes=c('external_gene_id')
                                                                      ,filters = 'ensembl_gene_id', values = x
                                                                      ,mart=mart)[1,1]))


# MA-plot
plotMA(r,ylim=c(-15,15))


# heatmap for gene read counts in different experiments
# sorted by adjusted p-values in increasing order just looking at the 20 most significant ones
library("pheatmap")
select <- order(r$padj,decreasing = F)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])
rows <- c()
for (k in 1:nrow(assay(ntd)[select,])){
  name <- (getBM(attributes=c('external_gene_id'),filters='ensembl_gene_id',
                 values = rownames(assay(ntd)[select,])[k],mart=mart)$external_gene_id)
  rows <- append(rows,name)
}
pheatmap(assay(ntd)[select,], cluster_rows=T, show_rownames=T,
         cluster_cols=F, annotation_col=df,labels_row=rows)

# heatmap for distances between experiments regarding the differences of their read counts
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$type, colnames(vsd), sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



