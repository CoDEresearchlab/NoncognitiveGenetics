# PACKAGES
library(semPlot)
library(tidySEM)
library(lavaan)
library(ggplot2)
library(dplyr)

##############
#show that increase in prediction is significant
##############

df <- read.csv('2021-04-26_DevNoncog_prepared-v4-latentfactors.csv')

df$NonCog <- df$NonCogFRC_1
df$Cog <- df$CogFRC_1
df$CogExt <- df$c_chol_cog_noncog_FRC_1
df$NonCogExt <- df$nc_chol_cog_noncog_FRC_1

df2 <- subset (df, selectunpaired == 1) #exclude 1 sib at rand among genotyped DZ pairs


#regress out covariates: 
df2$pcexgcsecorem1rs <- resid(lm(pcexgcsecorem1 ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1 + pcexgcseage1, df2, na.action=na.exclude))

#check correlation
#cor(df2$pcexgcsecorem1rs,df2$pcexgcsecorem1, use ='p')

df2$gta1rs <- resid(lm(gta1 ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1 + gtqage1, df2, na.action=na.exclude))

df2$it2a1rs <- resid(lm(it2a1 ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1 + itage1 , df2, na.action=na.exclude))

df2$lt2a1rs <- resid(lm(lt2a1 ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1 + ltqage1 , df2, na.action=na.exclude))

df2$CogExtrs <- resid(lm(CogExt ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1, df2, na.action=na.exclude))

df2$NonCogExtrs <- resid(lm(NonCogExt ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1, df2, na.action=na.exclude))


#SEM models 

#all estiamtes free 

model_1 <- "
gta1rs + it2a1rs + lt2a1rs + pcexgcsecorem1rs ~ CogExt+NonCogExt 
 
gta1rs ~~ it2a1rs + lt2a1rs + pcexgcsecorem1rs
it2a1rs ~~ lt2a1rs + pcexgcsecorem1rs
lt2a1rs ~~ pcexgcsecorem1rs
"

fit1 <- sem(model_1,data = df2)

#COG estimates constrianed to be equal
#NONcog estimates constrianed to be equal

model_2 <- "

 gta1rs  ~ a*CogExt 
  it2a1rs  ~ a*CogExt 
   lt2a1rs  ~ a*CogExt 
 pcexgcsecorem1rs ~ a*CogExt
 
  gta1rs  ~ b*NonCogExt 
  it2a1rs  ~ b*NonCogExt 
   lt2a1rs  ~ b*NonCogExt 
 pcexgcsecorem1rs ~ b*NonCogExt
 
 
gta1rs ~~ it2a1rs + lt2a1rs + pcexgcsecorem1rs
it2a1rs ~~ lt2a1rs + pcexgcsecorem1rs
lt2a1rs ~~ pcexgcsecorem1rs

"

fit2 <- sem(model_2,data = df2)

anova(fit1,fit2)

##only NON COG estimates constrained to be equal

model_3 <- "

 gta1rs  ~ CogExt 
  it2a1rs  ~ CogExt 
   lt2a1rs  ~ CogExt 
 pcexgcsecorem1rs ~ CogExt
 
  gta1rs  ~ b*NonCogExt 
  it2a1rs  ~ b*NonCogExt 
   lt2a1rs  ~ b*NonCogExt 
 pcexgcsecorem1rs ~ b*NonCogExt
 
 
gta1rs ~~ it2a1rs + lt2a1rs + pcexgcsecorem1rs
it2a1rs ~~ lt2a1rs + pcexgcsecorem1rs
lt2a1rs ~~ pcexgcsecorem1rs

"

fit3 <- sem(model_3,data = df2)

##only  COG estiamtes constrianed to be equal

model_4 <- "

gta1rs  ~ a*CogExt 
it2a1rs  ~ a*CogExt 
lt2a1rs  ~ a*CogExt 
pcexgcsecorem1rs ~ a*CogExt
 
gta1rs  ~ NonCogExt 
it2a1rs  ~ NonCogExt 
lt2a1rs  ~ NonCogExt 
pcexgcsecorem1rs ~ NonCogExt
 
 
gta1rs ~~ it2a1rs + lt2a1rs + pcexgcsecorem1rs
it2a1rs ~~ lt2a1rs + pcexgcsecorem1rs
lt2a1rs ~~ pcexgcsecorem1rs


"

fit4 <- sem(model_4,data = df2)

#summary(fit4, standardized = T)


#save results

write.table(round(anova(fit1,fit2),4), quote = F,sep=',')
write.table(round(anova(fit1,fit2,fit3),4), quote = F,sep=',')
write.table(round(anova(fit1,fit2,fit4),4), quote = F,sep=',')

#plot models 

#graph_sem(model = fit1, layout= )

edge_width = 2
node_width = 2
label_cex  = 2


label<- c("ea7", "ea9", "ea12","ea16","C PGS","NC PGS")


q1 <- semPaths(fit1, title=FALSE, 
               curvePivot = F, residuals = FALSE, exoCov = T,
               intercepts = F,
               edge.width = edge_width,
               node.width = node_width,
               label.cex = label_cex,
               style = "lisrel",
               #,
               #layout = "tree",
               layoutSplit = FALSE,
               edge.label.cex = 2,
               nodeLabels = label
               #,edge.label.position = c(.2,.2,.2,.2,.6,.6,.6,.6)
)


q2 <- semPaths(fit2, title=FALSE, 
               curvePivot = F, residuals = FALSE, exoCov = T,
               intercepts = F,
               edge.width = edge_width,
               node.width = node_width,
               label.cex = label_cex,
               style = "lisrel",
               #,
               #layout = "tree",
               layoutSplit = FALSE,
               edge.label.cex = 2
               ,edge.label.position = c(.2),
               nodeLabels = label
)


q3 <- semPaths(fit3, title=FALSE, 
               curvePivot = F, residuals = FALSE, exoCov = T,
               intercepts = F,
               edge.width = edge_width,
               node.width = node_width,
               label.cex = label_cex,
               style = "lisrel",
               #,
               #layout = "tree",
               layoutSplit = FALSE,
               edge.label.cex = 2
               ,edge.label.position = c(.2),
               nodeLabels = label
)

q4 <- semPaths(fit4, title=FALSE, 
               curvePivot = F, residuals = FALSE, exoCov = T,
               intercepts = F,
               edge.width = edge_width,
               node.width = node_width,
               label.cex = label_cex,
               style = "lisrel",
               #,
               #layout = "tree",
               layoutSplit = FALSE,
               edge.label.cex = 2
               ,edge.label.position = c(.2),
               nodeLabels = label
)



SIZE = 2.5

tiff('model-figures.png',  units="in", width=16, height=11, res=250)

par(mfrow=c(2,2)
    ,mar = c(2,2, 2, 2)
)    # set the plotting area into a 1*2 array

plot(q1)
title('a', adj = 0, cex.main=SIZE)
plot(q2)
title('b', adj = 0, cex.main=SIZE)
plot(q3)
title('c', adj = 0, cex.main=SIZE)
plot(q4)
title('d', adj = 0, cex.main=SIZE)
dev.off()


##############
#more fine grained comparisons for NonCog
##############

#only first two are equal
model_5 <- "

 gta1rs  ~ a*CogExt 
  it2a1rs  ~ a*CogExt 
   lt2a1rs  ~ a*CogExt 
 pcexgcsecorem1rs ~ a*CogExt
 
  gta1rs  ~ b*NonCogExt 
  it2a1rs  ~ b*NonCogExt 
   lt2a1rs  ~ NonCogExt 
 pcexgcsecorem1rs ~ NonCogExt
 
 
gta1rs ~~ it2a1rs + lt2a1rs + pcexgcsecorem1rs
it2a1rs ~~ lt2a1rs + pcexgcsecorem1rs
lt2a1rs ~~ pcexgcsecorem1rs

"

fit5 <- sem(model_5,data = df2)

#summary(fit5, standardized = T)

anova(fit4,fit5)


model_6 <- "

 gta1rs  ~ a*CogExt 
  it2a1rs  ~ a*CogExt 
   lt2a1rs  ~ a*CogExt 
 pcexgcsecorem1rs ~ a*CogExt
 
  gta1rs  ~ b*NonCogExt 
  it2a1rs  ~ b*NonCogExt 
   lt2a1rs  ~ b*NonCogExt 
 pcexgcsecorem1rs ~ NonCogExt
 
 
gta1rs ~~ it2a1rs + lt2a1rs + pcexgcsecorem1rs
it2a1rs ~~ lt2a1rs + pcexgcsecorem1rs
lt2a1rs ~~ pcexgcsecorem1rs

"

fit6 <- sem(model_6,data = df2)

#summary(fit6, standardized = T)

anova(fit4,fit5,fit6)

model_7 <- "

gta1rs  ~ a*CogExt 
it2a1rs  ~ a*CogExt 
lt2a1rs  ~ a*CogExt 
pcexgcsecorem1rs ~ a*CogExt

gta1rs  ~ b*NonCogExt 
it2a1rs  ~ NonCogExt 
lt2a1rs  ~ b*NonCogExt 
pcexgcsecorem1rs ~ NonCogExt
 
 
gta1rs ~~ it2a1rs + lt2a1rs + pcexgcsecorem1rs
it2a1rs ~~ lt2a1rs + pcexgcsecorem1rs
lt2a1rs ~~ pcexgcsecorem1rs

"

fit7 <- sem(model_7,data = df2)

#summary(fit6, standardized = T)

anova(fit4,fit6,fit7)

write.table(round(anova(fit4,fit5,fit6),4), quote = F,sep=',')
write.table(round(anova(fit4,fit6,fit7),4), quote = F,sep=',')

write.table(round(anova(fit1,fit2,fit3),4), quote = F,sep=',')
write.table(round(anova(fit1,fit2,fit4),4), quote = F,sep=',')


#graph_sem(model = fit1, layout= )

edge_width = 2
node_width = 2
label_cex  = 2

label<- c("ea7", "ea9", "ea12","ea16","C PGS","NC PGS")


q4 <- semPaths(fit4, title=FALSE, 
               curvePivot = F, residuals = FALSE, exoCov = T,
               intercepts = F,
               edge.width = edge_width,
               node.width = node_width,
               label.cex = label_cex,
               style = "lisrel",
               #,
               #layout = "tree",
               layoutSplit = FALSE,
               edge.label.cex = 2,
               nodeLabels = label
               #,edge.label.position = c(.2,.2,.2,.2,.6,.6,.6,.6)
)


q5 <- semPaths(fit5, title=FALSE, 
               curvePivot = F, residuals = FALSE, exoCov = T,
               intercepts = F,
               edge.width = edge_width,
               node.width = node_width,
               label.cex = label_cex,
               style = "lisrel",
               #,
               #layout = "tree",
               layoutSplit = FALSE,
               edge.label.cex = 2
               ,edge.label.position = c(.2),
               nodeLabels = label
)


q6 <- semPaths(fit6, title=FALSE, 
               curvePivot = F, residuals = FALSE, exoCov = T,
               intercepts = F,
               edge.width = edge_width,
               node.width = node_width,
               label.cex = label_cex,
               style = "lisrel",
               #,
               #layout = "tree",
               layoutSplit = FALSE,
               edge.label.cex = 2
               ,edge.label.position = c(.2),
               nodeLabels = label
)

q7 <- semPaths(fit7, title=FALSE, 
               curvePivot = F, residuals = FALSE, exoCov = T,
               intercepts = F,
               edge.width = edge_width,
               node.width = node_width,
               label.cex = label_cex,
               style = "lisrel",
               #,
               #layout = "tree",
               layoutSplit = FALSE,
               edge.label.cex = 2
               ,edge.label.position = c(.2),
               nodeLabels = label
)


SIZE = 2.5

tiff('model-figures-b.png',  units="in", width=16, height=11, res=250)

par(mfrow=c(2,2)
    ,mar = c(2,2, 2, 2)
)    # set the plotting area into a 1*2 array

plot(q4)
title('d', adj = 0, cex.main=SIZE)
plot(q5)
title('e', adj = 0, cex.main=SIZE)
plot(q6)
title('f', adj = 0, cex.main=SIZE)
plot(q7)
title('g', adj = 0, cex.main=SIZE)
dev.off()


###############
#show that increase in prediciton is significant after adjsuting for family SES
##############

df <- read.csv('2021-04-26_DevNoncog_prepared-v4-latentfactors.csv')

df$NonCog <- df$NonCogFRC_1
df$Cog <- df$CogFRC_1
df$CogExt <- df$c_chol_cog_noncog_FRC_1
df$NonCogExt <- df$nc_chol_cog_noncog_FRC_1

df2 <- subset(df, selectunpaired == 1) 

#regress out covariates: 

summary(lm(pcexgcsecorem1 ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1 + pcexgcseage1, df2, na.action=na.exclude))

df2$pcexgcsecorem1rs <- resid(lm(pcexgcsecorem1 ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1 + pcexgcseage1 + ases, df2, na.action=na.exclude))

#cor(df2$pcexgcsecorem1rs,df2$pcexgcsecorem1, use ='p')

df2$gta1rs <- resid(lm(gta1 ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1 + gtqage1 + ases, df2, na.action=na.exclude))

#cor(df2$gta1rs,df2$gta1, use ='p')

df2$it2a1rs <- resid(lm(it2a1 ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1 + itage1 + ases, df2, na.action=na.exclude))

#cor(df2$it2a1rs,df2$it2a1, use ='p')

df2$lt2a1rs <- resid(lm(lt2a1 ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1 + ltqage1 + ases, df2, na.action=na.exclude))

#cor(df2$lt2a1rs,df2$lt2a1, use ='p')

df2$CogExtrs <- resid(lm(CogExt ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1 + ases, df2, na.action=na.exclude))

#cor(df2$CogExtrs,df2$CogExt, use ='p')

df2$NonCogExtrs <- resid(lm(NonCogExt ~  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ Batch + sex1+ ases, df2, na.action=na.exclude))

#cor(df2$NonCogExtrs,df2$NonCogExt, use ='p')

fit1 <- sem(model_1,data = df2)
fit2 <- sem(model_2,data = df2)

anova(fit1,fit2)

fit3 <- sem(model_3,data = df2)

anova(fit1,fit2,fit3)

fit4 <- sem(model_4,data = df2)

#summary(fit4, standardized = T)

anova(fit1,fit2,fit4)

fit5 <- sem(model_5,data = df2)

#summary(fit5, standardized = T)

anova(fit1,fit2,fit4,fit5)

fit6 <- sem(model_6,data = df2)

#summary(fit6, standardized = T)

anova(fit1,fit2,fit4,fit5,fit6)

fit7 <- sem(model_7,data = df2)

#summary(fit6, standardized = T)

anova(fit1,fit2,fit4,fit6,fit7)

write.table(round(anova(fit1,fit2),4), quote = F,sep=',')
write.table(round(anova(fit1,fit2,fit3),4), quote = F,sep=',')
write.table(round(anova(fit1,fit2,fit4),4), quote = F,sep=',')

