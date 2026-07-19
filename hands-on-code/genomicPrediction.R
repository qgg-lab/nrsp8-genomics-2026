# load libraries
# BGLR is needed to read PLINK formatted file
# rrBLUP is needed to perform mixed model fitting
# ============================================================

library("BGLR")
library("rrBLUP")
setwd("geno-pheno") # also possible to change from the GUI console

# define the prefix for PLINK binary files
# documentation for plink format input can be found here
# https://www.cog-genomics.org/plink/1.9/input
# ?? work with AI to
# 1) read other formats into R
# 2) convert other formats to plink, then read into R
# ============================================================

bfile.prefix <- "final_extracted_50k"

plink.data <- read_bed(
  bed_file = paste0(bfile.prefix, ".bed"),
  bim_file = paste0(bfile.prefix, ".bim"),
  fam_file = paste0(bfile.prefix, ".fam")
)

str(plink.data)
head(plink.data$x, 20)

# raw genotypes are loaded as integers
# 0	= 00 Homozygote "1"/"1"
# 1	= 01 Heterozygote
# 2	= 10 Missing genotype
# 3	= 11 Homozygote "2"/"2"
# ============================================================

plink.data$x[plink.data$x == 2] <- NA
plink.data$x[plink.data$x == 0] <- 2
plink.data$x[plink.data$x == 1] <- 1
plink.data$x[plink.data$x == 3] <- 0

M_raw <- matrix(plink.data$x - 1, ncol = plink.data$p, nrow = plink.data$n, byrow = FALSE)
# this matrix contains n rows (number of animals)
# and p columns (number of markers)
# the plink.data$x - 1 is to turn the 0, 1, 2
# coding to -1, 0, 1 to be compatible with rrBLUP

# read phenotype
pheno <- read.table("simulated_phenotype.pheno.txt", header = FALSE, as.is = TRUE)

# genotype processing
# always do some summary statistics to check for data sanity
# ============================================================

# allele frequency
allele_freq <- (colMeans(M_raw) + 1) / 2
hist(allele_freq, breaks = seq(0, 1, 0.01))
hist(pheno[, 3])

# the population has 3000 individuals
# 1000 Durocs, 1000 Landraces, and 1000 Admixed
# we will focus on the 1000 Durocs, but the data are there
# for you to test additional scenarios
# e.g. what if all training animals are Durocs but
# the model has to be tested in Landraces

durocs.geno <- M_raw[1:1000, ]
durocs.pheno <- pheno[1:1000, 3]

# remove markers that have low MAF
durocs.freq <- (colMeans(durocs.geno) + 1) / 2
hist(durocs.freq, breaks = seq(0, 1, 0.01))
# note the difference between this distribution
# and the previous one on allele_freq for the whole population
# ?? does it matter in prediction?
durocs.geno <- durocs.geno[, durocs.freq >= 0.01]

# now perform 10 fold cross validation
# we will use GBLUP, feel free to also try SNP BLUP
# they should be equivalent
# ============================================================

# split data to 10 folds randomly
set.seed(1)
folds <- sample(rep(1:10, each = 100))

# create a vector to store the accuracy
pred.accu.10cv <- numeric(10)

# G matrix only needs to be calculated once
durocs.G <- A.mat(durocs.geno)
durocs.G[1:5, 1:5]

# before we do anything, let's just estimate h2
full.data.fit <- mixed.solve(y = durocs.pheno, K = durocs.G)
full.data.fit$Vu/(full.data.fit$Vu + full.data.fit$Ve)
# about 0.4 h2

# let's do first fold
mask.pheno <- durocs.pheno
mask.pheno[folds == 1] <- NA
this.fold.fit <- mixed.solve(y = mask.pheno, K = durocs.G)
plot(this.fold.fit$u[folds == 1], durocs.pheno[folds == 1])
pred.accu.10cv[1] <- cor(this.fold.fit$u[folds == 1], durocs.pheno[folds == 1])
  
# loop through the 10 folds
for (i in 2:10) {
  
  mask.pheno <- durocs.pheno
  mask.pheno[folds == i] <- NA
  this.fold.fit <- mixed.solve(y = mask.pheno, K = durocs.G)
  pred.accu.10cv[i] <- cor(this.fold.fit$u[folds == i], durocs.pheno[folds == i])
  cat("done fold", i, "\n")

}

# look at results
pred.accu.10cv

# if there is time, let's investigate effect
# of training set size, 100, 200, 400, 600, 900
# ============================================================

pred.accu.10cv.diff.size <- matrix(NA, nrow = 10, ncol = 5)
# each row is one fold, each column is one training size
trn.sizes <- c(100, 200, 400, 600, 900)

for (i in 1:10) {
  for (j in 1:5) {
    this.trn.size <- trn.sizes[j]
    # what animals to retain
    folds.to.retain <- c(i, setdiff(1:10, i)[1:as.integer(this.trn.size/100)])
    this.animals <- which(folds %in% folds.to.retain)
    this.folds <- folds[folds %in% folds.to.retain]
    this.pheno <- durocs.pheno[this.animals]
    this.mask.pheno <- this.pheno
    this.mask.pheno[this.folds == i] <- NA
    this.G <- durocs.G[this.animals, this.animals]
    this.fit <- mixed.solve(y = this.mask.pheno, K = this.G)
    pred.accu.10cv.diff.size[i, j] <- cor(this.fit$u[this.folds == i], this.pheno[this.folds == i])
    cat(i, j, "\n")
  }
}

# visualize results
boxplot(pred.accu.10cv.diff.size)
