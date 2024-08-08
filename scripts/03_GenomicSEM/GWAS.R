#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

GWASchunk = as.character(args[1])

library(data.table)
library(GenomicSEM)

sumstats <- fread(paste0("../out/",GWASchunk), data.table = F)

mldsc <- readRDS("../out/mldsc_cogNonCog30092020.rds")

model <- "
NC =~ NA*EDUC03 + INCO03 + HINCOME  
C =~ NA*INTE03 + MEMO + SymbolDigit + TMTB  + RT + EDUC03 + INCO03 + HINCOME  

NC~~1*NC
C~~1*C

NC ~~ 0*C

C ~ SNP
NC ~ SNP

EDUC03 ~~ 0*INTE03 + 0*MEMO + 0*SymbolDigit + 0*TMTB  + 0*RT  
INCO03 ~~  0*INTE03 + 0*MEMO + 0*SymbolDigit + 0*TMTB  + 0*RT  
HINCOME ~~ 0*INTE03 + 0*MEMO + 0*SymbolDigit + 0*TMTB  + 0*RT  

a> .001
b> .001
c> .001
d> .001
e> .001
f> .001
g> .001
h> .001

TMTB ~~ a*TMTB
EDUC03 ~~ b*EDUC03
INCO03 ~~ c*INCO03
INTE03 ~~ d*INTE03
HINCOME ~~ e*HINCOME
MEMO ~~ f*MEMO
RT ~~ g*RT
SymbolDigit ~~ h*SymbolDigit
"

#run the multivariate GWAS using parallel processing
directedFactors <- userGWAS(covstruc = mldsc, 
                            SNPs = sumstats, 
                            estimation = "DWLS", 
                            model = model, 
                            modelchi = TRUE,
                            printwarn = FALSE,
                            sub=c("C~SNP", "NC~SNP"), 
                            cores = NULL,
                            toler = FALSE,
                            SNPSE = 0.0005,
                            parallel = TRUE,
                            Output = NULL)

trait = c("C","NC")

for (sumstat in seq(2)){
  
  nm  <- trait[sumstat]
  
  write.table(directedFactors[[sumstat]], paste0("../out/",nm,".",GWASchunk,".txt"), col.names=T, row.names=F, quote=F)

  }