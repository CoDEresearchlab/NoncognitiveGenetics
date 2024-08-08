# format summary statistics 

library(GenomicSEM)

#sumstats names
trait.names<-c("EDUC03",
"INCO03",
"INTE03",
"HINCOME",
"MEMO",
"RT",
"SymbolDigit",
"TMTB")

#SE
se.logit <-c(F,F,F,F,F,F,F,F)

OLS <- c(T,T,T,T,T,T,T,T)

N=c(NA,NA,NA,NA,331679,330024,87741,78547)

files <- list("../sumstats/EA3.txt.gz", # Lee et al., 2018
 "../sumstats/Townsend.txt.gz", # Hill et al., 2019
  "../sumstats/IQ3_noTEDS.txt.gz", # Savage et al., 2018 - No TEDS
  "../sumstats/Income.txt.gz", # Hill et al., 2019
 "../sumstats/Memory.txt.gz", # de la Fuente et al., 2020
 "../sumstats/RT.txt.gz", # de la Fuente et al., 2020
 "../sumstats/SymbolDigit.txt.gz", # de la Fuente et al., 2020
 "../sumstats/TMTB.txt.gz") # de la Fuente et al., 2020

ref <- "../ref/reference.1000G.maf.0.005.txt"

formatted_sumstats <- sumstats(files=files,
        ref=ref,
        trait.names=trait.names,
        se.logit=se.logit,
        OLS = OLS, 
        N=N,
        maf.filter=0.01)


write.table(formatted_sumstats, "../out/formatted_sumstats_07102020.txt", col.names = T, row.names = F, quote = F)

