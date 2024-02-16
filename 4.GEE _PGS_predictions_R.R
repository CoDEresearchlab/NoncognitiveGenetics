rm(list = ls())
ls()

library(gee)
library(tidyverse)
library(psych)


df<- read.csv('', strip.white = T, 
              header = T, na.strings = NA)

#Recode Chip variable
levels(df$Chip)
df$chip_dummy<-ifelse(df$Chip == "affy        ", 0,
                      ifelse(df$Chip == "oee         ", 1,NA))

df$chip_dummy = as.numeric(paste(df$chip_dummy))

# Standardize  PCs
df$PC1 = scale(df$PC1)
df$PC2 = scale(df$PC2)
df$PC3 = scale(df$PC3)
df$PC4 = scale(df$PC4)
df$PC5 = scale(df$PC5)
df$PC6 = scale(df$PC6)
df$PC7 = scale(df$PC7)
df$PC8 = scale(df$PC8)
df$PC9 = scale(df$PC9)
df$PC10 = scale(df$PC10)


#Check variables distribution after preparation 
multi.hist(df[, c()])

# Select only participants with DNA data
df1<- df[c(df$genotyped1==1|df$genotyped2==1), ]

#Correlation bewteen Cog and Noncog PGS in the TEDS sammple 
x<- df1[, c( )]
#plot correlations
library(GGally)
cor(x, use = "complete") 
ggpairs(x, title= "correlaitons PGS Noncog Cog") 
ggcorr(x,label = TRUE)

#PGS correlations with sex (no corrleation as expected as PGS regressed fpr sex (sex1) 
#and standardized)
y<- df1[, c("NonCog", 'NonCogExt',
            "Cog",'CogExt','sex1', 'sex2')]
cor(y, use = "complete" )
ggcorr(y)


# Example of GEE Cognitive PGS prediction (using the cognitive score)including both cog and noncog scores in each regression)
# of achievement over development 
```{r, include =FALSE}

MAT <- matrix(NA,  nrow=12, ncol=5)# N col is based on the GEE output 
pval <- matrix(NA, nrow=12, ncol=1)# Create one extra matrix for p value (not provided by standard GEE output)
R2 <- matrix(NA, nrow=12, ncol=1) # Create an extra matrix for R2 (not provided by standard GEE output)

rownames(MAT) <- colnames( df[,c('gte1' ,  'gtm1' , 'gta1' ,  'ite1' ,  'itm1' ,  'it2a1', 'lte1' , 'ltm1' , 'lt2a1', 'pcexgcseengm1',  'pcexgcsematm1',  'pcexgcsecorem1') ]) # columns in your data frame (‘df’) corresponding to PGS you want to use for prediction : in my case corresponsing to criteria

for (i in 1:12){ #same as above: number of criteria
  
  x = df[,c('gte1' ,  'gtm1' , 'gta1' ,  'ite1' ,  'itm1' ,  'it2a1', 'lte1' , 'ltm1' , 'lt2a1', 'pcexgcseengm1',  'pcexgcsematm1',  'pcexgcsecorem1')][i] # same as above: position of PGS predictors in your data frame
  
  GEEfit <- gee(scale(x)~Cog+NonCog+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+chip_dummy,
                id = randomfamid, # this is the family ID in consecutive order 
                data = df, 
                corstr = "exchangeable")
  
  MAT[i,] <- coef(summary(GEEfit))[2,] # beta coefficient for x ~ PRS
  pval[i,] <- (2 * pnorm(abs(coef(summary(GEEfit))[,5]), lower.tail = FALSE))[2] # calculate tw-sided p-value based on robust SE
  R2[i,] <- coef(summary(GEEfit))[2] ^2 # calculate approximate R2
  
}

#This extract the results for the cog PGS prediciton. 

CogPGS_Cog_tabl <- cbind(MAT, pval[,1], R2[,1])
colnames(CogPGS_Cog_tabl) <- c('coeff', 'naive_SE','naive_Z','Robust_SE', 'Robust_Z','pval', 'R2')
CogPGS_Cog_tabl


#### Adjust p value for multipole testing ####
#Add the vector of p values

p<- c()
is.numeric(p)

p.corrected.fdr<- p.adjust(p, method = "fdr", n = length(p))
p.corrected.fdr = as.data.frame(p.corrected.fdr)

p.corrected.bon<- p.adjust(p, method = "bonferroni", n = length(p))
p.corrected.bon = as.data.frame(p.corrected.bon)
