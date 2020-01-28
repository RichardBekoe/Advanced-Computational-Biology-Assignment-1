load("snps.RData")
load("mat.gtex.RData")

load("gwas.als.RData")

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
# Plotting the distribution of gene expression levels across samples for gene CDKN2A:

# Output - Chi-squared test for given probabilities
# data:  table(round(snps[, 16]))
# X-squared = 1.2745, df = 2, p-value = 0.5287

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

# simple linear regression models to identify potential eQTLs

mod <- lm(expression.df$PIK3CA~ snps[,"rs7874974"])
summary(mod)
coef(summary(mod))
coef(summary(mod))[2,4]

# OUTPUT: > coef(summary(mod))[2,4]
# [1] 0.1476368








# Excercise 3

# Warning gwas.als 500,000 entries can take a long time to process

# Chr_17 Manhattan plot Base R Biomarker
Chr17_gwas<- subset(gwas.als, chr == 17)

plot(x = Chr17_gwas$bp,
     y=-log10(Chr17_gwas$p),
     pch='+',
     ylab='-log10(p) for biomarker',
     xlab='Position')


# Chr_19 Manhattan plot Base R Biomarker
Chr9_gwas<- subset(gwas.als, chr == 9)

plot(x = Chr9_gwas$bp,
     y=-log10(Chr9_gwas$p),
     pch='+',
     ylab='-log10(p) for biomarker',
     xlab='Position')


# ALL Manhattan plot Base R Biomarker
# plot(x = gwas.als$p,
#      y = -log10(gwas.als$p),
#      pch = '+',
#      ylab = '-log10(p) for PIK3CA expression',
#      xlab = 'Position')


# MIN all gwas
min.p.biomarker <- min(gwas.als$p)
min.p.biomarker
# [1] 1.70882e-24  (~20mins all processing time)
subset(gwas.als, p == min.p.biomarker)




##################################

##################################


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


# pvalues.PIK3CA[ (returned Low_p value(s)) ,  ] Search the row number and pass
# into above expression to return the relevant snp and PIK3CA.pvalues row


##################################

##################################


# Calculating Expression pvalues.CDKN2A Changed column name from original
# MarkerName in the script to snp (as the orginal had MarkerName as one of the
# columns for by.x part)
pvalues.CDKN2A <- data.frame(snp=colnames(snps), pvalues.CDKN2A=NA)

for (i in 1:nSNPs) {
  mod <- lm(expression.df$CDKN2A~snps[,i])
  pvalues.CDKN2A[i,2] <- coef(summary(mod))[2,4]
}

# Merging Calculated Expression pvalues.PIK3CA to gwas.als data frame
#helps for plotting histograms of pvalues.PIK3CA also
# 
# gwas.als <- merge(gwas.als, pvalues.CDKN2A,
#                   by.x = , by.y = "snp",
#                   all.x = FALSE, all.y = FALSE)
# head(gwas.als)


length(which(pvalues.CDKN2A$CDKN2A.pvalues <10^(-8)))

Low_c <- (which(pvalues.CDKN2A$CDKN2A.pvalues <10^(-8)))

Low_c

pvalues.CDKN2A[Low_c, ]


# OUTPUT: integer(0) pvalues which are below the threshold

###################################

##################################


# Calculating Expression pvalues.TP53 Changed column name from original
# MarkerName in the script to snp (as the orginal had MarkerName as one of the
# columns for by.x part)
pvalues.TP53 <- data.frame(snp=colnames(snps), TP53.pvalues=NA)

for (i in 1:nSNPs) {
  mod <- lm(expression.df$TP53~snps[,i])
  pvalues.TP53[i,2] <- coef(summary(mod))[2,4]
}

# Merging Calculated Expression pvalues.PIK3CA to gwas.als data frame
#helps for plotting histograms of pvalues.PIK3CA also
# 
# gwas.als <- merge(gwas.als, pvalues.TP53,
#                   by.x = , by.y = "snp",
#                   all.x = FALSE, all.y = FALSE)
# head(gwas.als)

length(which(pvalues.TP53$TP53.pvalues <10^(-8)))

Low_t <- (which(pvalues.TP53$TP53.pvalues <10^(-8)))

Low_t

pvalues.TP53[Low_t, ]


# OUTPUT: integer(0) pvalues which are below the threshold


###################################


##################################


# Calculating Expression pvalues.SMAD4 Changed column name from original
# MarkerName in the script to snp (as the orginal had MarkerName as one of the
# columns for by.x part)
pvalues.SMAD4 <- data.frame(snp=colnames(snps), SMAD4.pvalues=NA)

for (i in 1:nSNPs) {
  mod <- lm(expression.df$SMAD4~snps[,i])
  pvalues.SMAD4[i,2] <- coef(summary(mod))[2,4]
}

# Merging Calculated Expression pvalues.PIK3CA to gwas.als data frame
#helps for plotting histograms of pvalues.PIK3CA also
# 
# gwas.als <- merge(gwas.als, pvalues.SMAD4,
#                   by.x = , by.y = "snp",
#                   all.x = FALSE, all.y = FALSE)
# head(gwas.als)


length(which(pvalues.SMAD4$SMAD4.pvalues <10^(-8)))

Low_s <- (which(pvalues.SMAD4$SMAD4.pvalues <10^(-8)))

Low_s


# OUTPUT: integer(1) pvalues which are below the threshold


pvalues.SMAD4[ , ]

# snp SMAD4.pvalues
# 969 rs113505981  9.040235e-12

###################################








# ALL Manhattan plot Base R Expression

# plot(x = gwas.als$bp,
#      y = -log10(gwas.als$PIK3CA.pvalues),
#      pch = '+',
#      ylab = '-log10(p) for PIK3CA expression',
#      xlab = 'Position')


# Chr17_gwas Manhattan plot Base R Biomarker
plot(x = Chr17_gwas$bp ,
     y = -log10(Chr17_gwas$PIK3CA.pvalues),
     pch = '+',
     ylab = '-log10(p) for PIK3CA expression',
     xlab = 'Position')


# Chr9_gwas Manhattan plot Base R Biomarker
plot(x = Chr9_gwas$bp ,
     y = -log10(Chr9_gwas$PIK3CA.pvalues),
     pch = '+',
     ylab = '-log10(p) for PIK3CA expression',
     xlab = 'Position')


# MIN all gwas
min.p.biomarker.PIK3CA <- min(gwas.als$PIK3CA.pvalues)
min.p.biomarker.PIK3CA
subset(gwas.als, PIK3CA.pvalues == min.p.biomarker.sort1)





# ESSENTIAL FIRST

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
                      Genotype = snps[,"rs7874974"],
                      Site = expression.df$Tissue)

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

# The two models are not significantly different (0.05 threshold), so there is
# not a significant association between the selected SNP and PIK3CA expression
# after taking into account the site where measurements were performed. This
# result is different in comparison to the results in the previous findings


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





