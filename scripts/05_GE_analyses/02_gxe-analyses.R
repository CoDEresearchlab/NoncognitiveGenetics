#GxE analyses

#PACKAGES
library(ggplot2)
library(tidyverse)
library(lme4)

df <- read.csv("2021-04-26_DevNoncog_prepared-v4-latentfactors.csv")

df$CogExt <- df$c_chol_cog_noncog_FRC_1 # fraction 1 LDpred PGS - COG choleskyGWAS
df$NonCogExt <- df$nc_chol_cog_noncog_FRC_1 # fraction 1 LDpred PGS - NONCOG choleskyGWAS


resultsL <- list() #list of results

#educational achievement 
vars <- c("gta1","it2a1","lt2a1","pcexgcsecorem1") 


for (outcome in vars){ 
  
  #make sure all relevant vars for interaction are scaled
  varscld <- c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10',"ases","NonCogExt","CogExt")
  
  df[,c(varscld)]<- scale(df[,c(varscld)])
  
  indep <- c("ases","NonCogExt","CogExt") #main exposures 
  
  #cavoariates
  covars = c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','Batch','sex1')
  
  
  #test interactions between all exposures and covariates
  int <- apply(expand.grid(indep, covars), 1, function(x) paste(x[2], x[1], sep=":")) 
  
  #formula for lmer including main interactions and random effects
  f <- as.formula(paste(outcome,"~", paste(indep, collapse = '+'),'+', paste(int, collapse = '+'),'+',paste(covars, collapse = '+'),'+  ases:CogExt', '+  ases:NonCogExt', "+ (1 | randomfamid)"))
  
  fit <- summary(do.call("lmer", list(f, quote(df))))
  
  coefs  <- round(data.frame(coef(fit)),3) #RUN mixed model and extract coefficients
  
  resultsL[[paste0("lmer_x_" , outcome)]] <- coefs[c(1:4, 81,82),] #save main effects of exposures and interactions
  
  PVAL = 2*(pnorm(-abs(resultsL[[paste0("lmer_x_",i)]][["t.value"]]))) #PVAL
  
  resultsL[[paste0("lmer_x_", outcome)]][["p"]] <- signif(PVAL, digits=3) 
  
}


saveRDS(resultsL,"gxe.rds")


# the code below is only for plotting 

# categorize SES

KeepVars = c(i,"NonCogExt","CogExt","ases",covars)
dat <- df[complete.cases(df[,KeepVars]), ]

#make coarse SES variable
qSES <- quantile(dat$ases, na.rm = T) 

dat$sesCoars<- ifelse(dat$ases< qSES[2], "LOW SES", #lowest quantile
                      ifelse(dat$ases> qSES[4],"HIGH SES", #highest quantile
                             ifelse(dat$ases<= qSES[4] & dat$ases>= qSES[2], "MEDIUM SES", NA))) #middle 50%
dat$sesCoars = factor(dat$sesCoars, levels=c("LOW SES","MEDIUM SES","HIGH SES"), labels=c("LOW SES","MEDIUM SES","HIGH SES")) 


### std residual loop
# residulize outcome for either COG or NONCOG PGS + covars
for (ind in c("CogExt", "NonCogExt")){
  f <- as.formula(paste0(i, "~", ind,"+",paste(covars, collapse = '+'),"+ (1 | randomfamid)"))
  dat[,paste0("r",i,ind)] <- scale(resid(do.call("lmer", list(f, quote(dat)))))
  
  
  yCog <- paste0("r",i,"CogExt")
  yNCog <- paste0("r",i,"NonCogExt")
  
  #PLOT COG and NON COG adjusted slopes 
  
  x <- ggplot(dat) + 
    geom_jitter(aes_string("NonCogExt",noquote(yCog)), colour="#E69F00",alpha=0.5) + 
    geom_smooth(aes_string("NonCogExt",noquote(yCog)), method=lm, colour = "#D55E00", se=FALSE) +
    geom_jitter(aes_string("CogExt",noquote(yNCog)), colour="#56B4E9",alpha=0.3) + 
    geom_smooth(aes_string("CogExt",noquote(yNCog)), method=lm, colour = '#0072B2',  se=FALSE) +
    facet_wrap(~sesCoars, scales="free_x") +
    labs(x = "PGS", y = i) 
  
  png(paste0("plots/",i,"_lmer_x_adjusted_slopes.png"), res=300, height = 1200, width = 2000)
  plot(x)
  dev.off()
  
}
