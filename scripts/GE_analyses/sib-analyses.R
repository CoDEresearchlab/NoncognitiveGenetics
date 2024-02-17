## DZ analyses

#libraries 
library(lme4)
library(foreach)
library(ggplot2)

df <- read.csv('data_x_DZ/2021-04-26_DevNoncog_prepared-v4-latentfactors.csv')

#select only DZs
df <- df %>%
  filter(DZtwinpair == 1) 

### double enter PGS and covariates

newnames = c("CogFRC_1","NonCogFRC_1","c_chol_cog_noncog_FRC_1","nc_chol_cog_noncog_FRC_1","EA3_Lee2018_no23andme_FRCT1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","Chip","Batch")

cols = newnames
df2 <- df %>%
  group_by(randomfamid) %>%
  mutate_at(vars(cols), function(x){
    if (n() == 2){
      rev(x)
    } else {
      NA
    }
  }) %>%
  ungroup() %>%
  dplyr::select(newnames) %>%
  rename_all(funs(paste0(., "b"))) %>%
  cbind(df, .)


#change names for twin 1 (a) vs twin 2 (b)
oldnames = c("CogFRC_1","NonCogFRC_1","c_chol_cog_noncog_FRC_1","nc_chol_cog_noncog_FRC_1","EA3_Lee2018_no23andme_FRCT1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","Chip","Batch")

newnames = c("CogFRC_1a","NonCogFRC_1a","c_chol_cog_noncog_FRC_1a","nc_chol_cog_noncog_FRC_1a","EA3_Lee2018_no23andme_FRCT1a","PC1a","PC2a","PC3a","PC4a","PC5a","PC6a","PC7a","PC8a","PC9a","PC10a","Chipa","Batcha")

df2 <- df2 %>% rename_at(vars(oldnames), ~ newnames)

#check PGS correlations between siblings 
with(df2,cor(NonCogFRC_1b,NonCogFRC_1a,use = 'p'))
#[1] 0.5512216
with(df2,cor(CogFRC_1b,CogFRC_1a,use = 'p'))
#[1] 0.5401406


### create within/mean variables 

df3 <- df2 %>%
  mutate(meanCOG = .5*(CogFRC_1a+CogFRC_1b),
         meanNONCOG = .5*(NonCogFRC_1a+NonCogFRC_1b),
         meanCOGext = .5*(c_chol_cog_noncog_FRC_1a+c_chol_cog_noncog_FRC_1b),
         meanNONCOGext = .5*(nc_chol_cog_noncog_FRC_1a+nc_chol_cog_noncog_FRC_1b),
         
         withinCOG = CogFRC_1a-meanCOG,
         withinNONCOG = NonCogFRC_1a-meanNONCOG,
         withinCOGext = c_chol_cog_noncog_FRC_1a-meanCOGext,
         withinNONCOGext = nc_chol_cog_noncog_FRC_1a-meanNONCOGext,
         
         withinCOGb = CogFRC_1b-meanCOG,
         withinNONCOGb = NonCogFRC_1b-meanNONCOG,
         withinCOGextb = c_chol_cog_noncog_FRC_1b-meanCOGext,
         withinNONCOGextb = nc_chol_cog_noncog_FRC_1b-meanNONCOGext,
         
         meanEA = .5*(EA3_Lee2018_no23andme_FRCT1a+EA3_Lee2018_no23andme_FRCT1b),
         withinEA = EA3_Lee2018_no23andme_FRCT1a-meanEA,
         GCSEdiff = pcexgcsecorem1-pcexgcsecorem2,
         withinEAb = EA3_Lee2018_no23andme_FRCT1b-meanEA
  )


### for loop analyses

depvar = c("gta1","it2a1","lt2a1","pcexgcsecorem1")

reslts <- NULL

foreach::foreach(y=depvar) %do%
  {
covar = c("withinCOG","meanCOG"," withinNONCOG","meanNONCOG","PC1a","PC2a","PC3a","PC4a","PC5a","PC6a","PC7a","PC8a","PC9a","PC10a","Chipa","Batcha","sex1")
f <- as.formula(paste(y ,"~", paste(covar, collapse = "+"), "+ (1 | randomfamid)")) #FORMULA
coefs1 <- data.frame(coef(summary(do.call("lmer", list(f, quote(df3)))))) #RUN MIXZED MODEL AND EXTRACT COEF
coefs1$p.z <- 2*(pnorm(-abs(coefs1$t.value))) #PVAL

reslts <- rbind(reslts,coefs1[2:5,1:4])
  }

reslts$effect <- rownames(reslts) 
reslts$depvar <- rep(depvar, each=length(coefs1))

reslts$PRS <- rep(rep(c("cog","noncog"),each=2),length(depvar))
reslts$type <- rep(rep(c("within","between"),2),length(depvar))

reslts$depvar <- factor(reslts$depvar, levels=c("pcexgcsecorem1","lt2a1","it2a1","gta1"))

saveRDS(reslts,"dz.rds")

# format results
reslts$BETA<- round(reslts$Estimate,3)
reslts$SE<- round(reslts$Std..Error,3)
reslts$P<- round(reslts$p.z,3)
reslts$Z <- round(reslts$t.value,3)

write.table(reslts, "table_results-dz.txt", col.names=T, row.names=F, quote=F, sep=',') 


#plot
p <- reslts %>%
  #filter(type == "within") %>%
ggplot(aes(y=depvar,x=Estimate, color = PRS, group = PRS)) + 
  geom_point() + 
  geom_path() + 
  xlim(0,.5) + 
  facet_grid(~type) 
#+ 
    #geom_hline(aes(yintercept=.8), linetype="dashed", color = "blue")  

ggsave(p, file="withinBetween.png", width = 8, height=5, dpi=350)


# repeat analyses with Cholesky GWAS PGS

### for loop 

depvar = c("gta1","it2a1","lt2a1","pcexgcsecorem1")

reslts <- NULL

foreach::foreach(y=depvar) %do%
  {
    covar = c("withinCOGext","meanCOGext"," withinNONCOGext","meanNONCOGext","PC1a","PC2a","PC3a","PC4a","PC5a","PC6a","PC7a","PC8a","PC9a","PC10a","Chipa","Batcha","sex1")
    f <- as.formula(paste(y ,"~", paste(covar, collapse = "+"), "+ (1 | randomfamid)")) #FORMULA
    coefs1 <- data.frame(coef(summary(do.call("lmer", list(f, quote(df3)))))) #RUN MIXZED MODEL AND EXTRACT COEF
    coefs1$p.z <- 2*(pnorm(-abs(coefs1$t.value))) #PVAL
    
    reslts <- rbind(reslts,coefs1[2:5,1:4])
    
  }

reslts$effect <- rownames(reslts) 
reslts$depvar <- rep(depvar, each=length(coefs1))

reslts$PRS <- rep(rep(c("cog","noncog"),each=2),length(depvar))
reslts$type <- rep(rep(c("within","between"),2),length(depvar))

reslts$depvar <- factor(reslts$depvar, levels=c("pcexgcsecorem1","lt2a1","it2a1","gta1"))

saveRDS(reslts,"dz.ext.rds")

# format results
reslts$BETA<- round(reslts$Estimate,3)
reslts$SE<- round(reslts$Std..Error,3)
reslts$P<- round(reslts$p.z,3)
reslts$Z <- round(reslts$t.value,3)

write.table(reslts, "table_results-dz.ext.txt", col.names=T, row.names=F, quote=F, sep=',') 


#plot 

p <- reslts %>%
  #filter(type == "within") %>%
  ggplot(aes(y=depvar,x=Estimate, color = PRS, group = PRS)) + 
  geom_point() + 
  geom_path() + 
  xlim(0,.5) + 
  facet_grid(~type) #+
#geom_hline(aes(yintercept=.8), linetype="dashed", color = "blue")  

ggsave(p, file="data_x_DZ/withinBetween_EXT.png", width = 8, height=5, dpi=350)


