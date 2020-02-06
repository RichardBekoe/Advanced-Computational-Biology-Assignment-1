
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
# which are below the 10^(-8) cutoff. Then finding the corresponding SNPs
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

#Same snp, rs113505981, as PIK3CA

# OUTPUT: integer(1) pvalues which are below the threshold

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
# # Chi-squared test for given probabilities
# data:  table(round(snps[, "rs113505981"]))
# X-squared = 81.345, df = 2, p-value < 2.2e-16

# This SNP_1 does not comply with HWE (p-value < 2.2e-16)




selectedSNP_PIK3CA_2 <- snps[,"rs4977264"]
roundSNP_2 <- round(snps[,"rs4977264"])


MAF.SNP_2 <- sum(table(roundSNP_2)*c(0, 1, 2))/(2*length(selectedSNP_PIK3CA_2))
MAF.SNP_2

# No reversal required as < 0.5. The MAF for the SNP is <0.5,
# confirming that the allele for which we return the measure is the minor (=
# less common) allele. In an eQTL study often a minimum MAF is required.

theory2 <- c((1-MAF.SNP_2)^2, 2*MAF.SNP_2*(1-MAF.SNP_2), MAF.SNP_2^2)
chisq.test(table(round(snps[,"rs4977264"])), p=theory2)

# Chi-squared test for given probabilities
# data:  table(round(snps[, "rs4977264"]))
# X-squared = 0.96215, df = 2, p-value = 0.6181

# This SNP_2 complies with HWE (p-value = 0.6181)

# Plotted are the distribution of gene expression levels across samples for gene PIK3CA
# and SMAD4. Hence illustrating identified associations.

# Base R hist PIK3CA
hist(expression.df$PIK3CA, main="Gene expression profile: PIK3CA", xlab = "Expression level")

library(ggplot2)

# ggplot hist PIK3CA
ggplot(expression.df, aes(PIK3CA))+
  geom_histogram(aes(y=..density..), colour='black', fill="white")+
  geom_density(alpha=.2, fill="#FF6666")


# Base R hist SMAD4
hist(expression.df$SMAD4, main="Gene expression profile: SMAD4",
     xlab="Expression level")

# ggplot hist SMAD4
ggplot(expression.df, aes(SMAD4))+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
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

# grouped pancreas, brain and the other tissues
# repeat code for debugging purposes; easy delete/ override
expression.df <- as.data.frame(mat.gtex)

expression.df$Tissue <- "Other"
expression.df[which(expression.df$PIK3CA>10),]$Tissue <-"Brain"
expression.df[which(expression.df$PIK3CA<6),]$Tissue <-"Pancreas"

install.packages("data.table")

# Count number of repeats of the word "Other" then -- pass into function which replaces the
# word and assigns it with a random number from 1 to 4

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

# Set seed of randomisation so that this is reproducible:
# could use loop! Or sample as characters instead of numbers - expression.df[
# expression.df == "Other"] <- sample((lung, colon, oesophagus, kidney),
# numOther, replace=TRUE)
# where Other = lung, colon, oesophagus, kidney

expression.df$Tissue <- factor(expression.df$Tissue)
head(expression.df)

dat.covariate <- data.frame(Expression = expression.df$PIK3CA,
                            Genotype = snps[,"rs4977264"],
                            Site = expression.df$Tissue)

# "rs4977264" = selectedSNP_PIK3CA_2

linearm <- lm(Expression~Genotype+Site, data = dat.covariate)
summary(linearm)

# Residual standard error: 1.935 on 46 degrees of freedom
# Multiple R-squared:  0.8675,	Adjusted R-squared:  0.8502 
# F-statistic:  50.2 on 6 and 46 DF,  p-value: < 2.2e-16

# We can see that there is a significant association between the SNP (rs4977264) and PIK3CA
# expression, but also between colon, kidney, lung, oesophagus, and pancreas expresion.

library(lme4)

mixedm <- lmer (Expression~Genotype+(1|Site), data = dat.covariate)
summary(mixedm)

# The model does not output a p-value, so to get the p-value for the
# fixed effects, we can perform a likelihood ratio test using anova()

# Compare our model with that of a null model where the random effect would be factored in:
reducedm <- lmer(Expression~1+(1|Site), data = dat.covariate)
anova(mixedm, reducedm)

# mixedm: Expression ~ Genotype + (1 | Site)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# reducedm  3 246.64 252.55 -120.32   240.64                            
# mixedm    4 239.79 247.67 -115.89   231.79 8.8477      1   0.002935 **

# The two models are significantly different (0.05 threshold), so there is a
# significant association between the selected SNP and PIK3CA expression after
# taking into account the site where measurements were performed. This result is
# similar in comparison to the results in the previous findings. Mixed effect models, take into account
# various random effects not causal to the trait in question, such as the tissue type.

# As we can see the slope is different when covariates are incorporated in the
# model. This addition modifies the estimated slope and its associated p-value,
# and highlights 


# The lower p-value when using covariates may indicate that
# the regression describes the data more accurately.






# Q2 a) & b) ------------------------------------------------------

# BIOMARKER analysis

# Below the threshold for gwas.als

length(which(gwas.als$p <10^(-8)))

# Output: 124 pvalues which are below the threshold

Low_p_all <- (which(gwas.als$p <10^(-8)))

Low_p_all

gwas.als[Low_p_all, ]

als_p.threshold_val <- gwas.als[Low_p_all, ]

# 124 in total, below the threshold
# 121 chr 9
# 3 chr 17


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

# ALL Manhattan plot Base R Biomarker

plot(x = gwas.als$p,
     y = -log10(gwas.als$p),
     pch = '+',
     ylab = '-log10(p) for gwas alas expression',
     xlab = 'Position')


# MIN all gwas
min.p.biomarker <- min(gwas.als$p)
min.p.biomarker
# [1] 1.70882e-24
subset(gwas.als, p == min.p.biomarker)

#         chr   snp       bp     a1 a2   freq         b         se           p
# 4500536   9 rs3849943 27543382  C  T 0.247955 0.0405097 0.00396593 1.70882e-24

# Minimum gwas.als p value (biomarker) snp rs3849943 found on CHR 9



# ALL data: Manhattan plot Base R Biomarker qqman

install.packages("qqman")

library(qqman)

manhattan(gwas.als, chr="chr", bp="bp", p="p", snp="snp")

# All_manhattan <- manhattan(gwas.als, chr="chr", bp="bp", p="p", snp="snp")

# explanation : What do you observe? Is there anything particularly striking?
# How can you explain these observations? (10 points for interpretation) 

# The plot supports the above results # showing a peak at chr 17 (3points) and a
# larger peak at chr 9 (121 points). The X axis shows that the position on the
# chromosome where these significant SNPs occur are at the same position. The Y
# axis tells how much it is associated with a trait, these vary in range surpassing the threshold.
# Noticeable also is that there is a gap in the position of chromosome 9 where there are no points at all.

#  research histogram -> add to the explanation



# Q3 a) b) c) -----------------------------------------------------

# Selecting the top most significant SNP; One method is to subset the respective
# chromosome then, manually filter the data by size order in the "p" column.
# Hence then to show the corresponding snp using the data from Question two.


# CHROMOSOME 17)  Selecting the top most significant SNP from chrosomosome 17.

#         chr     snp     bp       a1 a2    freq     b         se           p
# 6843269	17	rs35714695	26719788	G	A	0.8239640	0.0295130	0.00455194	8.95546e-11

# CHROMOSOME 9)  Selecting the top most significant SNP from chromosome 9.

#          chr    snp        bp     a1  a2    freq       b              se           p
# 
# 4500536	 9	 rs3849943	27543382	C	  T	  0.2479550	  0.0405097	  0.00396593	  1.70882e-24



# Used the above snps from chromosome 17 and chromosome 9 to search
# https://gtexportal.org/home/ for single-tissue eQTLs of these SNPs. And
# downloaded the results in csv format. Also applied a threshold to identify
# eQTLs

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


# For Chr 9 one gene is affected C9orf72 and the tissues affected are as seen in the table.

# For Chr 17 two genes TMEM97 and POLDIP2 are affected and the tissues affected are as seen in the table


# b)     How many cis and how many trans eQTLs in each case? (5points)

# Chr 9 
# All SNP.ID rs3849943, gene C9orf72 have Cis-eQTLs (number = 18) [	rs3849943]

# Chr 17
# All SNP.ID rs35714695, gene TMEM97 have Cis-eQTLs (number = 1) [ rs35714695]
#
# All SNP.ID rs35714695, gene POLDIP2 have Cis-eQTLs (number = 2) [ref - Systematic identification database… ] [ rs35714695]

# A cis eQTL is a SNP showing association with gene expression located
# close by (usually within a 1-Mb window) of the gene that it influences. A
# trans eQTL is a SNP located elsewhere in the genome, away from the gene whose
# expression it is altering. GTEx IGV Browser was used to explore the location
# of the SNP in comparison to the transcription start sites of the genes, hence
# evaluate the relative distance An online eQTL browser was also used as another
# source of information. [ref - https://gtexportal.org/home/] [ref - systematic identification database… ]

# In addiotion, gene expression signatures are cell-type specific, and therefore
# regulatory control of expression may also be cell-type dependent. Significant
# tissue specificity has been reported for multiple cis eQTLs.


# c)     Comment on whether these genes could be involved in the development of ALS and how
# that might look like at a functional level, providing further support from scientific literature or
# other sources. (10 points)



# POLDIP2
# The brain cis-eQTL effect for SNPs in this locus on POLDIP2 suggests that
# POLDIP2 could be the causal gene in this locus. Another overlap was observed
# in the SARM1 locus where rs35714695 and its proxies had the strongest
# exon-level cis-eQTL effect on POLDIP2 in multiple brain tissues (P = 2.32 ×
# 10−3). [ ref - Genome-wide association ]


# TMEM97
# TMEM97 ligands bind to S2R and has a pharmacologic profile
# the same as that of S2R. Over the past few decades, sigma receptors (SRs),
# including sigma 1 and sigma 2 receptor subtypes (S1R and S2R, respectively)
# have been widely associated with aging- and mitochondria-associated disorders,
# such as Parkinson’s and Alzheimer’s disease, multiple sclerosis and
# amyotrophic lateral sclerosis. However, the specific role played by this orphan
# receptor family in cell biology has yet to be clarified. The 3D structure of
# TMEM97/S2R would provide understanding of the biological functions and
# mechanisms. [ref - Sigma receptors… ] [ref - Role of the sigma-1 receptor…]


# C9orf72
# The main functional disease mechanisms proposed have been the loss of function
# of the C9orf72 protein and toxic gain of function from C9orf72 repeat RNA or
# from dipeptide repeat proteins produced by repeat-associated non-ATG
# translation; more than a few hundred repeats represent a risk for ALS and FTD.
# The discovery that repeat expansions in the C9orf72 gene are a frequent cause
# of amyotrophic lateral sclerosis (ALS) has been important in understanding the
# pathogenic mechanisms of the disease. Several downstream processes across a
# range of cellular functions have also been implicated. [ref - C9orf72-mediated
# ALS and FTD: multiple pathways to disease ] [ref - Molecular Mechanisms of
# Neurodegeneration Related to C9orf72 Hexanucleotide Repeat Expansion]


