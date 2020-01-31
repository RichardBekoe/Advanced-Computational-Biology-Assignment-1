# Question 1

# a)

load("snps.RData")
load("mat.gtex.RData")

load("gwas.als.RData")


#  Quality control checks

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

selectedSNP <- snps[,16]
geno.16 <- round(snps[,16])
# "rs7874974" = snp 16


MAF.16 <- sum(table(geno.16)*c(0, 1, 2))/(2*length(selectedSNP))
MAF.16
theory <- c((1-MAF.16)^2, 2*MAF.16*(1-MAF.16), MAF.16^2)
chisq.test(table(round(snps[,16])), p=theory)

# This SNP complies with HWE.
# Now plotted is the gene expression distribution for selected genes.
# Plotting the distribution of gene expression levels across samples for gene CDKN2A: ?

# Output - Chi-squared test for given probabilities
# data:  table(round(snps[, 16]))
# X-squared = 1.2745, df = 2, p-value = 0.5287




# simple linear regression models to identify potential eQTLs

mod <- lm(expression.df$PIK3CA~ snps[,"rs7874974"])
summary(mod)
coef(summary(mod))
coef(summary(mod))[2,4]

# OUTPUT: > coef(summary(mod))[2,4]
# [1] 0.1476368




# Calculating Expression PIK3CA.pvalues Changed column name from original
# MarkerName in the script to snp (as the orginal had MarkerName as one of the
# columns for by.x part)
pvalues.PIK3CA <- data.frame(snp=colnames(snps), PIK3CA.pvalues=NA)

for (i in 1:nSNPs) {
  mod <- lm(expression.df$PIK3CA~snps[,i])
  pvalues.PIK3CA[i,2] <- coef(summary(mod))[2,4]
}

# Merging Calculated Expression pvalues.PIK3CA to gwas.als data frame
#helps for plotting histograms of pvalues.PIK3CA also

gwas.als <- merge(gwas.als, pvalues.PIK3CA,
                  by.x = , by.y = "snp",
                  all.x = FALSE, all.y = FALSE)
head(gwas.als)


length(which(pvalues.PIK3CA$PIK3CA.pvalues <10^(-8)))

# Output: 2;  pvalues which are below the threshold

Low_p <- (which(pvalues.PIK3CA$PIK3CA.pvalues <10^(-8)))

Low_p

pvalues.PIK3CA[Low_p, ]

# OUTPUT: 
# snp PIK3CA.pvalues
# 964 rs4977264   1.063495e-11
# 
# snp PIK3CA.pvalues
# 969 rs113505981   4.725661e-14

# After do the MAF method and check which ones comply


##################################

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

###################################

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

###################################


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


# OUTPUT: integer(1) pvalues which are below the threshold


pvalues.SMAD4[Low_s, ]

# snp SMAD4.pvalues
# 969 rs113505981  9.040235e-12


###################################

##################################

# MAF?? each one pass?



# Warning gwas.als 500,000 entries can take a long time to process

# BIOMARKER
# Below the threshold for gwas.als

length(which(gwas.als$p <10^(-8)))

# Output: 124;  pvalues which are below the threshold

Low_p_all <- (which(gwas.als$p <10^(-8)))

Low_p_all

gwas.als[Low_p_all, ]

als_p.threshold_val <- gwas.als[Low_p_all, ]

# lm for each gene
# ????? QUESTION 1 Boxplot for pik3ca


# 124 in total, below the threshold
# 121 chr 9
# 3 chr 17


# Graphs of identified associations

# Base R hist PIK3CA
hist(expression.df$PIK3CA, main="Gene expression profile: PIK3CA", xlab = "Expression level")

library(ggplot2)

# ggplot hist PIK3CA
ggplot(expression.df, aes(PIK3CA))+
  geom_histogram(aes(y=..density..), colour='black', fill="white")+
  geom_density(alpha=.2, fill="#FF6666")



# Base R hist CDKN2A
hist(expression.df$CDKN2A, main="Gene expression profile: CDKN2A",
     xlab="Expression level")

# ggplot hist CDKN2A
ggplot(expression.df, aes(CDKN2A))+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")



# boxplot PIK3CA ~ geno.16
boxplot(expression.df$PIK3CA ~ geno.16,
        main="Gene expression levels of PIK3CA for SNP rs7874974",
        xlab = "Genotypes",
        ylab = "Expression",
        names = c("CC", "CT", "TT"))



# Question 1)

# b)



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

# the significant snip

# "rs7874974" = snp 16

linearm <- lm(Expression~Genotype+Site, data = dat.covariate)
summary(linearm)

# Residual standard error: 2.01 on 46 degrees of freedom
# Multiple R-squared:  0.857,	Adjusted R-squared:  0.8384 
# F-statistic: 45.96 on 6 and 46 DF,  p-value: < 2.2e-16

# We can see that there is a significant association between the SNP and PIK3CA
# expression, but also between "colon, kidney, lung, oesophagus, Pancreas" and
# Pancreas and PIK3CA expression.

library(lme4)

mixedm <- lmer (Expression~Genotype+(1|Site), data = dat.covariate)
summary(mixedm)

# The model does not output a p-value, so to get the p-value for the
# fixed effects, we can perform a likelihood ratio test using anova()

# Compare our model with that of a null model where the random effect would be factored in:
reducedm <- lmer(Expression~1+(1|Site), data = dat.covariate)
anova(mixedm, reducedm)

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# reducedm  3 246.78 252.69 -120.39   240.78                         
# mixedm    4 246.26 254.14 -119.13   238.26 2.5145      1     0.1128

0.002935 **

# The two models are not significantly different (0.05 threshold), so there is
# not a significant association between the selected SNP and PIK3CA expression
# after taking into account the site where measurements were performed. This
# result is different in comparison to the results in the previous findings


# Simple regression models or correlation statistics
# to assess association between genotype and phenotype, don’t account for population structure.
# Therefore, more commonly used in practice are mixed effect models, which take into account
# various random effects not causal to the trait in question, such as the tissue type


# As we can see the slope is different when covariates are incorporated in the
# model. This addition modifies the estimated slope and its associated p-value,
# and highlights the number of days of treatment as a significant contributor to
# the association found. The lower p-value when using covariates indicates that
# the regression describes the data more accurately.


# We can see that there is a significant association between the SNP and X?
# expression (as before), but also between "SITE" and "SITE" and "X" expression.




# moodle reply help
# snps[,"rsxx1"] - where "rsxx1" is your actual SNP ID. Another way to select it
# in a more general fashion is snps[,which(colnames(snps) == "rsxx1")]. If you
# just want the column number, then which(colnames(snps) == "rsxx1") should give
# you your answer.


















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
# [1] 1.70882e-24  (~20mins all processing time)
subset(gwas.als, p == min.p.biomarker)






# Question 3 - assessment
# Selecting the top most significant SNP; one method is to subset the respective
# chromosome then, manually filter the data by size order in the "p" column hence to show the corresponding snp.

# Selecting the top most significant SNP from chrosomosome 17 is

# 
#         chr     snp     bp       a1 a2    freq     b         se           p
# 6843269	17	rs35714695	26719788	G	A	0.8239640	0.0295130	0.00455194	8.95546e-11 - from file
# 
#             rs35714695	chr17	28392769	true	G	A	17_26719788_G_A_b37 - from online source
# 
# 
# # top most significant SNP from chromosome 9.
# 
#          chr    snp        bp     a1  a2    freq       b              se           p
# 
# 4500536	 9	 rs3849943	27543382	C	  T	  0.2479550	  0.0405097	  0.00396593	  1.70882e-24







# ALL Manhattan plot Base R Biomarker qqman

install.packages("qqman")

library(qqman)

manhattan(gwas.als, chr="chr", bp="bp", p="p", snp="snp")

# All_manhattan <- manhattan(gwas.als, chr="chr", bp="bp", p="p", snp="snp")

# explanation : What do you observe? Is there anything particularly striking?
# How can you explain these observations? (10 points for interpretation) 

# The plot supports the above results # showing a peak at chr 17 (3points) and
# larger peak at chr 9 (121 points). The X axis shows that the position on a
# chromosome where these significant SNPs occur are at the same position. The Y
# axis tells how much it is associated with a trait, these vary in range beyond the threshold.
# Noticeable also is that there is a gap in the position of chromosome 9 where there are points at all.



  


# EXPRESSION

# ALL Manhattan plot Base R Expression; PIK3CA.pvalues

# plot(x = gwas.als$bp,
#      y = -log10(gwas.als$PIK3CA.pvalues),
#      pch = '+',
#      ylab = '-log10(p) for PIK3CA expression',
#      xlab = 'Position')


# Chr17_gwas Manhattan plot Base R Expression     
plot(x = Chr17_gwas$bp ,
     y = -log10(Chr17_gwas$PIK3CA.pvalues),
     pch = '+',
     ylab = '-log10(p) for PIK3CA expression',
     xlab = 'Position')


# Chr9_gwas Manhattan plot Base R Expression     
plot(x = Chr9_gwas$bp ,
     y = -log10(Chr9_gwas$PIK3CA.pvalues),
     pch = '+',
     ylab = '-log10(p) for PIK3CA expression',
     xlab = 'Position')


# MIN all gwas
min.p.biomarker.PIK3CA <- min(gwas.als$PIK3CA.pvalues)
min.p.biomarker.PIK3CA
subset(gwas.als, PIK3CA.pvalues == min.p.biomarker.sort1)

##################################

##################################










# ESSENTIAL FIRST assignment sheet questions

install.packages("MatrixEQTL")
library(MatrixEQTL)

pvOutputThreshold = 1e-8
errorCovariance = numeric()
useModel = modelLINEAR

ssnps = SlicedData$new()
ssnps$CreateFromMatrix(t(snps))
show(ssnps)

genes = SlicedData$new()

mat.expr <- t(as.matrix(mat.gtex[,-1]))
rownames(mat.expr) <- NULL
colnames(mat.expr) <- rownames(mat.gtex)


# Alternative; one gene column only e.g. PIK3CA
# mat.expr <- t(as.matrix(mat.gtex$PIK3CA[,-1]))
# rownames(mat.expr) <- NULL
# colnames(mat.expr) <- rownames(mat.gtex$PIK3CA)

# OR for just one gene instead of the whole matrix e.g. PIK3CA ??
# VALIDATION model should match up with other method ??

genes$CreateFromMatrix(mat.expr)
show(genes)

# cvrt = SlicedData$new()
# mat.cov<- t(as.matrix(cov[,-c(1,3)]))
# colnames(mat.cov) <- cov$sample
# cvrt$CreateFromMatrix(mat.cov)
# show(cvrt)

meqt11 = Matrix_eQTL_engine(
  snps = ssnps,
  gene = genes,
  output_file_name = NULL,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

plot(meqt11)

meqt11$all$neqtls

# [1] 1

head(meqt11$all$eqtls)

# snps   gene        statistic        pvalue        FDR           beta
# 1 rs113505981 row3  8.777644      9.040235e-12    2.712071e-08    12.04821

# We can see that there is a significant association between the SNP and FRK
# expression (as before), but also between sites 5 and 6 and FRK expression.

# We can see that this output matches the number of eQTLs originally found for mat.gtex (with/without adjusting for covariates).
# Print out the number of eQTLs and their details









# OLD - simple linear regression models to identify potential eQTLs
dat.snps16 <- data.frame(Expression = expression.df$PIK3CA,
                      Genotype = snps[,16])

linearm <- lm(Expression~Genotype, data = dat.snps16)
summary(linearm)

# coef(summary(mod))[2,4] ??? USE MOD!! AND LOOP, ALL THE GENES; SAME AS BEFORE CALCULATION?? WITH CORRELATION??


# OUTPUT: Residual standard error: 4.946 on 51 degrees of freedom
# Multiple R-squared:  0.04066,	Adjusted R-squared:  0.02185 
# F-statistic: 2.162 on 1 and 51 DF,  p-value: 0.1476

# We can see that there is a significant association between the SNP and PIK3CA expression

dat.snps20 <- data.frame(Expression = expression.df$CDKN2A,
                         Genotype = snps[,20])

linearm <- lm(Expression~Genotype, data = dat.snps20)
summary(linearm)

# We can not see that there is a significant association between the SNP and CDKN2A expression

# OUTPUT: Residual standard error: 2.721 on 51 degrees of freedom
# Multiple R-squared:  0.008847,	Adjusted R-squared:  -0.01059 
# F-statistic: 0.4552 on 1 and 51 DF,  p-value: 0.5029

dat.snps35 <- data.frame(Expression = expression.df$SMAD4,
                         Genotype = snps[,35])

linearm <- lm(Expression~Genotype, data = dat.snps35)
summary(linearm)

# OUTPUT: Residual standard error: 9.267 on 51 degrees of freedom
# Multiple R-squared:  0.02209,	Adjusted R-squared:  0.002918
# F-statistic: 1.152 on 1 and 51 DF,  p-value: 0.2882











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



# Which genes are affected and in which tissues  Imported17.threshold_val
# 
# For Chr 9 one gene is affected C9orf72 and the tissues affected are as seen in the table
# 
# For Chr 17 two genes are affected and the tissues affected are as seen in the table

# How many cis and how many trans eQTLs in each case? (5points)

# Chr 9
# All SNP.ID C9orf72 have Cis-eQTLs (number = 18) [	rs3849943]
# 
# Chr 17
# 
# All TMEM97 have Cis-eQTLs (number = 1) [ rs35714695]
# 
# All POLDIP2 have Cis-eQTLs (number = 2) [ref - Systematic identification… ] [ rs35714695]


# Gene expression signatures are cell-type specific, and therefore regulatory
# control of expression may also be cell-type dependent. Significant tissue
# specificity has been reported for multiple cis eQTLs.


# Comment on whether these genes could be involved in the development of ALS and how
# that might look like at a functional level, providing further support from scientific literature or
# other sources. (10 points)



# POLDIP2
# The brain cis-eQTL effect for SNPs in this locus on POLDIP2 suggests that
# POLDIP2 could be the causal gene in this locus. Another overlap was observed
# in the SARM1 locus where rs35714695 and its proxies had the strongest
# exon-level cis-eQTL effect on POLDIP2 in multiple brain tissues (P = 2.32 ×
# 10−3). [ ref - Genome-wide association ]


# TMEM97
# TMEM97 ligands bind to S2R (Alon et al., 2017) and has a pharmacologic profile
# the same as that of S2R. Over the past few decades, sigma receptors (SRs),
# including sigma 1 and sigma 2 receptor subtypes (S1R and S2R, respectively)
# have been widely associated with aging- and mitochondria-associated disorders,
# such as Parkinson’s and Alzheimer’s disease, multiple sclerosis and
# amyotrophic lateral sclerosis However, the specific role played by this orphan
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







# SEE ONENOTE FOR FURTHER NOTES



