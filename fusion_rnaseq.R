
# read output of STAR-Mapper from files and extract their 4th column (wanted results) and combine them
# in a matrix than can be used for DESeq2 analysis 
# 1. read 48h files
setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/mapping/48hMEF/RNASeq/")
files <- list.files()
out <- matrix(ncol=12,nrow=37991)
out <- as.data.frame(out, stringsAsFactors = FALSE)
cols <- c()
rows <- c()
i <- 2

for (filename in files){
  if(!grepl("PE",filename)){
    a <- read.table(filename)
    a <- a[c(-4:-1),]
    out[,i] <- a[,4]
    cols <- append(cols,filename)
    i <- i+1
  }
  else{
    a <- read.table(filename)
    a <- a[c(-4:-1),]
    out[,i] <- a[,4]
    cols <- append(cols,filename)
    i <- i+1
  }
  
}
out[,1] <- a[,1]

# 2. read ESC files
setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/mapping/ESC/RNASeq/")
files <- list.files()
for (filename in files){
  if(!grepl("PE",filename)){
    a <- read.table(filename)
    a <- a[c(-4:-1),]
    out[,i] <- a[,4]
    cols <- append(cols,filename)
    i <- i+1
  }
  else{
    a <- read.table(filename)
    a <- a[c(-4:-1),]
    out[,i] <- a[,4]
    cols <- append(cols,filename)
    i <- i+1
  }
  
}

# trim column names 
cols <- sub("_mRNA-Seq.ReadsPerGene.out.tab","",cols)
cols <- sub("_mRNA-SeqReadsPerGene.out.tab","",cols)
cols <- sub("_mRNA-Seq_1ReadsPerGene.out.tab","",cols)
cols <- sub("48_","",cols)
rownames(out)<- out[,1]
out <- out[,-1]
colnames(out) <- cols

# write resulting fusion table containing data of all RNA-Seq experiment of both stages in one matrix to a table 
setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/mapping/")
write.table(out,"fusion.txt",quote=F)
