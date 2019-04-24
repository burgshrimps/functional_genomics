library(biomaRt)
library(GenomicRanges)
library(GenomicFeatures)

## read RNASeq files of 48hMEF with counts and merge into one matrix 
setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/mapping/48hMEF/RNASeq")
files <- list.files()
out_48h <- c()
cols <- c()
for (filename in files){
  a <- read.table(filename)
  a <- a[c(-4:-1),]
  out_48h <- cbind(out_48h,a[,4])
  cols <- append(cols,filename)
}
rownames(out_48h) <- a[,1]
cols <- sub("_mRNA-Seq.ReadsPerGene.out.tab","",cols)
cols <- sub("_mRNA-SeqReadsPerGene.out.tab","",cols)
cols <- sub("_mRNA-Seq_1ReadsPerGene.out.tab","",cols)
colnames(out_48h) <- cols

## read RNASeq files of ESC with counts and merge into one matrix
setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/mapping/ESC/RNASeq")
files <- list.files()
out_ESC <- c()
cols <- c()
for (filename in files){
  a <- read.table(filename)
  a <- a[c(-4:-1),]
  out_ESC <- cbind(out_ESC,a[,4])
  cols <- append(cols,filename)
}
rownames(out_ESC) <- a[,1]
cols <- sub("_mRNA-Seq.ReadsPerGene.out.tab","",cols)
cols <- sub("_mRNA-SeqReadsPerGene.out.tab","",cols)
cols <- sub("_mRNA-Seq_1ReadsPerGene.out.tab","",cols)
colnames(out_ESC) <- cols


# make mart object for the mm9 genome and get all protein coding genes
ensembl67 <- useMart(host='may2012.archive.ensembl.org',
                     biomart='ENSEMBL_MART_ENSEMBL',
                     dataset='mmusculus_gene_ensembl')
prot.gene <- getBM(attributes=c("ensembl_gene_id", "gene_biotype", "strand",  
                                "chromosome_name", "start_position", "end_position")
                   ,mart = ensembl67, filters = "biotype", values = "protein_coding")
convert.strand <- function(strand.col){
  xx <-strand.col
  xx <- as.factor(xx)
  levels(xx)[levels(xx)=="-1"] <- "-"
  levels(xx)[levels(xx)=="1"] <- "+"
  levels(xx)[!levels(xx)%in%c("-","+")] <- "*"
  return(xx)
}
###
correct.chr <- function(genes){
  genes <- genes[genes$chromosome_name%in%c(1:19,"X","Y"),]
  genes$chromosome_name <- paste0("chr",genes$chromosome_name)
  return(genes)
}

prot.gene$strand <- convert.strand(prot.gene$strand)
prot.gene <- correct.chr(prot.gene)





# RPKM-normmalization for gene-length of RNASeq-Data

# get all exons corresponding to all ensembl genes #
mm9.exons <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                                 dataset="mmusculus_gene_ensembl",
                                 transcript_ids=NULL,
                                 circ_seqs=DEFAULT_CIRC_SEQS,
                                 filter=NULL,
                                 id_prefix="ensembl_",
                                 host="may2012.archive.ensembl.org", #www.ensembl.org",
                                 port=80,
                                 taxonomyId=NA,
                                 miRBaseBuild=NA)
# now get the exons per gene (list of genomic ranges)
exonic <- exonsBy(mm9.exons, by="gene")
# reduce the exons by the union (list of genomic ranges)
red.exonic <- reduce(exonic)
# lengts of all genes as a sum of exons
exon.lengths <- sum(width(red.exonic))




## calculate and replace counts with RPKM values 
out <- cbind(out_48h,out_ESC)
for (r in 1:nrow(out)){
  g <- which(names(exon.lengths) == rownames(out)[r])
  out[r,] <- log2(((out[r,]/exon.lengths[g])*1000)+1)
}


norm.counts <- as.matrix(out)

# extract 48h bzw. ESC columns respectively
norm.counts_48h <- norm.counts[,1:4]
norm.counts_ESC <- norm.counts[,5:6]

# calculate means for every row, which means summarizing counts in an average count per gene per stage
norm.counts_48h <- rowMeans(norm.counts_48h)
norm.counts_ESC <- rowMeans(norm.counts_ESC)

m_48h <- norm.counts_48h
m_ESC <- norm.counts_ESC



## filter out rows with genes, that are also in the norm-table with histone-modifications
##   -> from prot.gene GRanges (ensembl67)
n <- c()
for (r in 1:nrow(prot.gene)){
  l <- which(names(m_48h) == prot.gene[r,1])
  n <- append(n,l)
}
n <- as.numeric(n)
m_48h <- m_48h[n]
m_ESC <- m_ESC[n]


write.table(m_48h,"med_RNA_48h_RPKM.txt", row.names=T, quote=F, col.names = F)
write.table(m_ESC,"med_RNA_ESC_RPKM.txt", row.names=T, quote=F, col.names = F)

