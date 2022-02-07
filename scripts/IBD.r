# IBD O. colligata
## - pairwise relatedness as measured with beta in hierfstat
## - US samples geographically measured via Chinese sample
## - dbmem by RDA

library(vcfR)
library(adegenet)
library(hierfstat)
library(geodist)
library(ecodist)
library(adespatial)
library(vegan)

# load coordinates of Oc
xyOc <- read.delim("coords.txt")
names(xyOc) <- c("ID", "latitude", "longitude")

# load SNP data
vcf <- read.vcfR("multisample_HighCov_withCN.vcf", verbose = FALSE)

# convert data to genind
genind_O <- vcfR2genind(vcf)
# add pop labels
pop(genind_O)<-substr(indNames(genind_O),1,5)

# convert to hierfstat
hfs_O <- genind2hierfstat(genind_O)

# mask more than 2-allelic sites
hfs_O[hfs_O == 3] <- NA

# calculate pair-wise relatedness
CoAncestOc <- betas(hfs_O[,-1], diploid = FALSE, betaijT = TRUE)

# measure US samples geographically via Chinese sample
xyOc[8:10,2:3] <- xyOc[1, 2:3]
Oc.dist <- as.dist(geodist(xyOc, measure = "geodesic")) # exact measure
Oc.dist.matrix <- as.matrix(Oc.dist)
Oc.dist.matrix[8:10,1:7] <- as.matrix(Oc.dist)[8:10,1:7]+1.195936e+07 # add distance from Chinese to US samples to all US distances

# linear regression for regression line in figure
summary(lm( unlist(as.dist(CoAncestOc))~unlist(as.dist(Oc.dist.matrix)))) 

#Call:
#  lm(formula = unlist(as.dist(CoAncestOc)) ~ unlist(as.dist(Oc.dist.matrix)))
#
#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.65690 -0.08406  0.08181  0.25089  0.39972 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                      5.841e-01  8.107e-02   7.205 6.51e-09 ***
#  unlist(as.dist(Oc.dist.matrix)) -6.231e-08  6.874e-09  -9.066 1.55e-11 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Residual standard error: 0.3289 on 43 degrees of freedom
#Multiple R-squared:  0.6565,    Adjusted R-squared:  0.6485 
#F-statistic: 82.19 on 1 and 43 DF,  p-value: 1.553e-11

# dbMEM by RDA
CoAncestOc.pc <- prcomp(as.dist(CoAncestOc)) # r-stats
xyOc.dbmem <- dbmem(as.dist(Oc.dist.matrix))
Oc.rda<-rda(CoAncestOc.pc$x, xyOc.dbmem)
RsquareAdj(Oc.rda)
#$r.squared
#[1] 0.5389867
#
#$adj.r.squared
#[1] 0.48136
anova(Oc.rda, perm=1000)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: rda(X = CoAncestH.pc$x, Y = xyHam.dbmem)
#Df Variance      F Pr(>F)    
#Model     2  2.34039 12.635  0.001 ***
#  Residual  7  0.64829                  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Figure with regression line and rda statistics
pdf("Oc_IBD_pacific.pdf", width = 7, height = 7)                                
plot(unlist(as.dist(Oc.dist.matrix)), unlist(as.dist(CoAncestOc)), xlab = "Pairwise distance (1,000 km)", ylab = "Pairwise relatedness", xaxt = "n", cex.lab=1.5, cex.axis=1.5, cex=1.5, pch=16)
axis(side = 1, at = c(0, 5000000, 10000000, 15000000, 20000000), labels = c("0","5","10","15","20"), tck = -0.01, cex.lab=1.5, cex.axis=1.5)
abline(5.841e-01, -6.231e-08, lwd = 1.5, col="darkgrey")
text(15000000, 0.8, expression(paste(italic("R")^2 ~ " = 0.54, ", italic("p")," = 0.001", sep="")), cex = 1.5)
dev.off()