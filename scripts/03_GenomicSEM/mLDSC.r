#multivariable LDSC

library(GenomicSEM)

ld <- "../ref/eur_w_ld_chr/"
wld <- "../ref/eur_w_ld_chr/"

traits<-c("../sumstats/munged/EA3.munged.gz", # Lee et al., 2018
          "../sumstats/munged/Townsend.munged.gz", # Hill et al., 2019
          "../sumstats/munged/IQ3_noTEDS.munged.gz", # Savage et al., 2018 - No TEDS
          "../sumstats/munged/Income.munged.gz", # Hill et al., 2019
          "../sumstats/munged/Memory.munged.gz", # de la Fuente et al., 2020
          "../sumstats/munged/RT.munged.gz", # de la Fuente et al., 2020
          "../sumstats/munged/SymbolDigit.munged.gz", # de la Fuente et al., 2020
          "../sumstats/munged/TMTB.munged.gz") # de la Fuente et al., 2020


sample.prev <- c(
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA)

population.prev <- c(
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA)

trait.names<-c("EDUC03",
               "INCO03",
               "INTE03",
               "HINCOME",
               "MEMO",
               "RT",
               "SymbolDigit",
               "TMTB")

mldsc <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)

saveRDS(mldsc, file="../out/mldsc_cogNonCog30092020.rds")