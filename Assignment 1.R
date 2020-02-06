
# Q1 a) -----------------------------------------------------------

load("snps.RData")
load("mat.gtex.RData")

load("gwas.als.RData")


# Quality control checks were performed

# Conversion to dataframe
expression.df <- as.data.frame(mat.gtex)

class(snps)
class(mat.gtex)
nSNPs <- ncol(snps)
nSNPs

samples.GenoData <- rownames(snps)
samples.ExprData <- rownames(mat.gtex)
setdiff(samples.GenoData, samples.ExprData)
setdiff(samples.ExprData, samples.ExprData)

length(intersect(samples.GenoData,samples.ExprData))
length(intersect(samples.GenoData,samples.ExprData)) == length(rownames(mat.gtex))

# Use of linear regression models to identify potential eQTLs
# Identifying potential eQTLs for all genes

# Calculating Expression PIK3CA.pvalues levels and thresholding to find those
# which are below a 10^(-8) cutoff. Then finding the corresponding SNPs
pvalues.PIK3CA <- data.frame(snp=colnames(snps), PIK3CA.pvalues=NA)

for (i in 1:nSNPs) {
  mod <- lm(expression.df$PIK3CA~snps[,i])
  pvalues.PIK3CA[i,2] <- coef(summary(mod))[2,4]
}

length(which(pvalues.PIK3CA$PIK3CA.pvalues <10^(-8)))
Low_p <- (which(pvalues.PIK3CA$PIK3CA.pvalues <10^(-8)))
Low_p
pvalues.PIK3CA[Low_p, ]

# snp PIK3CA.pvalues
# 964 rs4977264   1.063495e-11
# 969 rs113505981   4.725661e-14

# Output: 2 pvalues which are below the threshold
##################################

# Calculating Expression pvalues.CDKN2A
pvalues.CDKN2A <- data.frame(snp=colnames(snps), pvalues.CDKN2A=NA)

for (i in 1:nSNPs) {
  mod <- lm(expression.df$CDKN2A~snps[,i])
  pvalues.CDKN2A[i,2] <- coef(summary(mod))[2,4]
}

length(which(pvalues.CDKN2A$CDKN2A.pvalues <10^(-8)))
Low_c <- (which(pvalues.CDKN2A$CDKN2A.pvalues <10^(-8)))
Low_c
pvalues.CDKN2A[Low_c, ]

# OUTPUT: integer(0) pvalues which are below the threshold
##################################

# Calculating Expression pvalues.TP53
pvalues.TP53 <- data.frame(snp=colnames(snps), TP53.pvalues=NA)

for (i in 1:nSNPs) {
  mod <- lm(expression.df$TP53~snps[,i])
  pvalues.TP53[i,2] <- coef(summary(mod))[2,4]
}

length(which(pvalues.TP53$TP53.pvalues <10^(-8)))
Low_t <- (which(pvalues.TP53$TP53.pvalues <10^(-8)))
Low_t
pvalues.TP53[Low_t, ]

# OUTPUT: integer(0) pvalues which are below the threshold
##################################

# Calculating Expression pvalues.SMAD4
pvalues.SMAD4 <- data.frame(snp=colnames(snps), SMAD4.pvalues=NA)

for (i in 1:nSNPs) {
  mod <- lm(expression.df$SMAD4~snps[,i])
  pvalues.SMAD4[i,2] <- coef(summary(mod))[2,4]
}

length(which(pvalues.SMAD4$SMAD4.pvalues <10^(-8)))
Low_s <- (which(pvalues.SMAD4$SMAD4.pvalues <10^(-8)))
Low_s
pvalues.SMAD4[Low_s, ]

# snp SMAD4.pvalues
# 969 rs113505981  9.040235e-12

# OUTPUT: integer(1) pvalues which are below the threshold
#Same snp, rs113505981, as PIK3CA
##################################

selectedSNP_PIK3CA_1 <- snps[,"rs113505981"]
roundSNP_1 <- round(snps[,"rs113505981"])
AF.SNP_1 <- sum(table(roundSNP_1)*c(0, 1, 2))/(2*length(selectedSNP_PIK3CA_1))
AF.SNP_1
# reversal required as AF greater than 0.5.

MAF.SNP_1 <- min(AF.SNP_1, 1-AF.SNP_1)
MAF.SNP_1
theory1 <- c((1-MAF.SNP_1)^2, 2*MAF.SNP_1*(1-MAF.SNP_1), MAF.SNP_1^2)
chisq.test(table(round(snps[,"rs113505981"])), p=theory1)
# Chi-squared test for given probabilities
# data:  table(round(snps[, "rs113505981"]))
# X-squared = 81.345, df = 2, p-value < 2.2e-16

# This SNP_1 does not comply with HWE (p-value < 2.2e-16)

selectedSNP_PIK3CA_2 <- snps[,"rs4977264"]
roundSNP_2 <- round(snps[,"rs4977264"])
MAF.SNP_2 <- sum(table(roundSNP_2)*c(0, 1, 2))/(2*length(selectedSNP_PIK3CA_2))
MAF.SNP_2
# No reversal required as < 0.5. The MAF for the SNP is <0.5,
# confirming that the allele for which we return the measure is the minor (=
# less common) allele.

theory2 <- c((1-MAF.SNP_2)^2, 2*MAF.SNP_2*(1-MAF.SNP_2), MAF.SNP_2^2)
chisq.test(table(round(snps[,"rs4977264"])), p=theory2)
# Chi-squared test for given probabilities
# data:  table(round(snps[, "rs4977264"]))
# X-squared = 0.96215, df = 2, p-value = 0.6181

# This SNP_2 complies with HWE (p-value = 0.6181)
# Plotted are the distribution of gene expression levels across samples for gene PIK3CA

# Base R hist PIK3CA
hist(expression.df$PIK3CA, main="Gene expression profile: PIK3CA", xlab = "Expression level")

library(ggplot2)
# ggplot hist PIK3CA
ggplot(expression.df, aes(PIK3CA))+
  geom_histogram(aes(y=..density..), colour='black', fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
# Influence of SNP rs4977264 on phenotype
# boxplot PIK3CA ~ "rs4977264"

boxplot(expression.df$PIK3CA ~ roundSNP_2,
        main="Gene expression levels of PIK3CA for SNP rs4977264",
        xlab = "Genotypes",
        ylab = "Expression",
        names = c("CC", "CT", "TT"))

# Boxplot of gene expression levels of PIK3CA for SNP rs4977264
# Searched gwas.als using snp rs4977264 to find the corresponding alleles to
# find C and T.


# Q1 b) -----------------------------------------------------------

# grouped together pancreas, brain and the other tissues
expression.df <- as.data.frame(mat.gtex)

expression.df$Tissue <- "Other"
expression.df[which(expression.df$PIK3CA>10),]$Tissue <-"Brain"
expression.df[which(expression.df$PIK3CA<6),]$Tissue <-"Pancreas"

install.packages("data.table")

library(data.table)
tissue_repeats <- setDT(expression.df)[, .N, Tissue]
nOther <- (tissue_repeats)[2,2]
numOther <- as.numeric(nOther)

set.seed(6537332)
expression.df[ expression.df == "Other"] <- sample(1:4, numOther, replace=TRUE)

expression.df[ expression.df == 1] <- "lung"
expression.df[ expression.df == 2] <- "colon"
expression.df[ expression.df == 3] <- "oesophagus"
expression.df[ expression.df == 4] <- "kidney"

# Count number of repeats of the word "Other" then - pass into function which
# replaces the word and assigns it with a random number from 1 to 4
# Set seed of randomisation so that this is reproducible 
# where Other = lung, colon, oesophagus, kidney

expression.df$Tissue <- factor(expression.df$Tissue)
head(expression.df)
dat.covariate <- data.frame(Expression = expression.df$PIK3CA,
                            Genotype = snps[,"rs4977264"],
                            Site = expression.df$Tissue)

# "rs4977264" = selectedSNP_PIK3CA_2 (the snp which passed HWE previously)

linearm <- lm(Expression~Genotype+Site, data = dat.covariate)
summary(linearm)

# Residual standard error: 1.935 on 46 degrees of freedom
# Multiple R-squared:  0.8675,	Adjusted R-squared:  0.8502 
# F-statistic:  50.2 on 6 and 46 DF,  p-value: < 2.2e-16

# We can see that there is a significant association between the SNP (rs4977264)
# and PIK3CA expression, but also between colon, kidney, lung, oesophagus, and
# pancreas expresion.

library(lme4)
mixedm <- lmer (Expression~Genotype+(1|Site), data = dat.covariate)
summary(mixedm)

# The model does not output a p-value, so to get the p-value for the
# fixed effects, we can perform a likelihood ratio test using anova()

# We compared our model with that of a null model where the random effect is factored in:
reducedm <- lmer(Expression~1+(1|Site), data = dat.covariate)
anova(mixedm, reducedm)

# mixedm: Expression ~ Genotype + (1 | Site)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# reducedm  3 246.64 252.55 -120.32   240.64                            
# mixedm    4 239.79 247.67 -115.89   231.79 8.8477      1   0.002935 **

# The two models are significantly different (0.05 threshold), so there is a
# significant association between the selected SNP and PIK3CA expression after
# taking into account the site where measurements were performed. This result is
# similar in comparison to the results in the previous findings. Mixed effect
# models, take into account various random effects not causal to the trait in
# question.

# The lower p-value when using different covariates may indicate that
# the regression describes the data more accurately.



# Q2 a) & b) ------------------------------------------------------

# BIOMARKER analysis
# Below the threshold for gwas.als

length(which(gwas.als$p <10^(-8)))
Low_p_all <- (which(gwas.als$p <10^(-8)))
Low_p_all
gwas.als[Low_p_all, ]
als_p.threshold_val <- gwas.als[Low_p_all, ]

# Output: 124 pvalues which are below the threshold
# 121 for chr 9; 3 for chr 17

# Chr_17 Manhattan plot Base R Biomarker
Chr17_gwas<- subset(gwas.als, chr == 17)
plot(x = Chr17_gwas$bp,
     y=-log10(Chr17_gwas$p),
     pch='+',
     ylab='-log10(p) for biomarker',
     xlab='Position')

# Chr_9 Manhattan plot Base R Biomarker
Chr9_gwas<- subset(gwas.als, chr == 9)
plot(x = Chr9_gwas$bp,
     y=-log10(Chr9_gwas$p),
     pch='+',
     ylab='-log10(p) for biomarker',
     xlab='Position')

# ALL data: Manhattan plot Base R Biomarker qqman
install.packages("qqman")
library(qqman)
manhattan(gwas.als, chr="chr", bp="bp", p="p", snp="snp")

# The plot supports the above results showing a peak at chr 17 (3points) and a
# larger, more dense peak at chr 9 (121 points). The X axis shows that the
# position on the chromosome where these significant SNPs are located are at
# similar position. The Y axis tells how much it is associated with a trait,
# these vary in range surpassing the threshold. Noticeable also is that there is
# a gap in the positions of chromosome 9 where there are no points at all.


# Q3 a) b) c) -----------------------------------------------------

# Selecting the top most significant SNP; One method is to subset the respective
# chromosome then manually filter the data by size order in the "p" column.
# Hence showing the corresponding snp, using the data from Question two.


# CHROMOSOME 17)  Selecting the top most significant SNP from chrosomosome 17
# equals rs35714695. CHROMOSOME 9)  Selecting the top most significant SNP from
# chromosome 9 equals rs3849943. Used the above snps from chromosome 17 and
# chromosome 9 to search (https://gtexportal.org/home/) for single-tissue eQTLs
# of these SNPs. Then downloaded the results in csv format. Also applied a
# threshold of (10^-8) to identify eQTLs.

Imported17_snp = read.csv("17GTExPortal.csv", header = TRUE)
length(which(Imported17_snp$P.Value <10^(-8)))
Imported17_low <- (which(Imported17_snp$P.Value <10^(-8)))
Imported17_low
Imported17_snp[Imported17_low, ]
Imported17.threshold_val <- Imported17_snp[Imported17_low, ]


Imported9_snp = read.csv("9GTExPortal.csv", header = TRUE)
length(which(Imported9_snp$P.Value <10^(-8)))
Imported9_low <- (which(Imported9_snp$P.Value <10^(-8)))
Imported9_low
Imported9_snp[Imported9_low, ]
Imported9.threshold_val <- Imported9_snp[Imported9_low, ]

# For Chr 9 one gene is affected, C9orf72 and the tissues affected are as seen
# in the table above. For Chr 17 two genes, TMEM97 and POLDIP2 are affected and
# the tissues affected are as seen in the table above.

# Chr 9: All SNP.ID rs3849943, gene C9orf72 have Cis-eQTLs (number = 18)
# Chr 17: All SNP.ID rs35714695, gene TMEM97 have Cis-eQTLs (number = 1)
# All SNP.ID rs35714695, gene POLDIP2 have Cis-eQTLs (number = 2)