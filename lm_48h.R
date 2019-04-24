## read norm-Table with normalized histon-modification counts in promoter regions (already normalized by control)
## and matrix with RPKM-RNASeq counts 
setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/linear_regression/")
norm <- read.table("norm_48h.txt",header = T)
m <- read.table("med_RNA_48h_RPKM.txt")
names <- m[,1]
m <- as.numeric(m[,2])
names(m) <- names



## randomly select "training" (1), "validation" (2), "test" (3) for every index/row 
## and split corresponding data into these categories 
random <- sample(x= 1:3, size=22369, prob=c(0.5,0.25,0.25), replace=T)

training <- c()
val <- c()
test <- c()
training <- m[(random == 1)]
val <- m[(random == 2)]
test <- m[(random == 3)]




# linear model (lm) fitting for expression values dependent on histon modifications 
norm <- as.matrix(norm)
norm_glm <- norm[(random == 1),]
data  <- data.frame(cbind(training,norm_glm))

model <- lm(training ~ H3.3 + H3K27ac + H3K27me3 + H3K36me3 + H3K4me1 +  H3K4me2 +  
              H3K4me3 + H3K79me2 + H3K9ac + H3K9me3 + H3 + Hdac1 + 0,data = data)
coefs <- model$coefficients

summary(model)
q <- summary(model)
rsq <- q[["r.squared"]]


# greatest beta-coefficients and corresponding Histone modifications
max.coefs <- tail(sort(coefs),n=3)


# smallest beta-coefficients and corresponding Histone modifications
min.coefs <- head(sort(coefs),n=3)

# create scatterplots to visualize correlation between specific histone modifications and RNASeq transcript counts 
library(affy)
scatter1 <- smoothScatter((norm[,12]),(m),xlab = "Hdac1", ylab = "normalized RNA-Counts")
scatter1 <- smoothScatter((norm[,7]),(m),xlab = "H3K4me3", ylab = "normalized RNA-Counts")
scatter1 <- smoothScatter((norm[,2]),(m),xlab = "H3K27ac", ylab = "normalized RNA-Counts")
# calculate specific correlation values for every considered histone modification and print these 
for (c in 1:12){
  print(paste(colnames(norm)[c],cor(norm[,c],m),sep=": "))
}


# prediction of transcript counts by using the build model and calculate differences / residuals to actual values 
norm_test = norm[(random == 3),]
pr <- predict(model, newdata = data.frame(H3.3 = norm_test[,1],
                                          H3K27ac = norm_test[,2],
                                          H3K27me3 = norm_test[,3],
                                          H3K36me3 = norm_test[,4],
                                          H3K4me1 = norm_test[,5],
                                          H3K4me2 = norm_test[,6],
                                          H3K4me3 = norm_test[,7],
                                          H3K79me2 = norm_test[,8],
                                          H3K9ac = norm_test[,9],
                                          H3K9me3 = norm_test[,10],
                                          H3 = norm_test[,11],
                                          Hdac1 = norm_test[,12]))

# plot residuals showing the values of differences between predicted RNA-counts and actual values from the RNASeq-experiments
plot((pr-test), pch = 20 )

# average absolute difference 
mean(abs(pr-test))

# histogram to show the distribution of differences
hist((pr-test), breaks=seq(-15,16, l = 50), col = "darkblue", main ="residuals 48h model validation", 
     xlab ="residuals", cex.lab=1.4, cex.axis = 1.2, cex.main = 1.6)



