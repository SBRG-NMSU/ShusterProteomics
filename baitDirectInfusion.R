########### Prereqs ###########
options(stringsAsFactors = FALSE, scipen = 900)
oldPar <- par()
os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/ShusterProteomics/", 
                  "~/gdrive/ShusterProteomics/")
setwd(baseDir)
library(MSGFplus)
library(tidyverse)
library(mzID)

########### Setup parameters and run ###########
par <- msgfPar()

dbFile <- paste0(baseDir, "UP000007110.fasta")
db(par) <- dbFile
tolerance(par) <- '10 ppm'
chargeRange(par) <- c(2, 6)
lengthRange(par) <- c(6, 25)
instrument(par) <- 'HighRes'
enzyme(par) <- 'Trypsin'
fragmentation(par) <- 0
protocol(par) <- 0
isotopeError(par) <- c(0, 1)
matches(par) <- 6
ntt(par) <- 1
tda(par) <- TRUE
mods(par)[[1]] <- msgfParModification(name = 'Carbamidomethyl', 
                                      composition = 'C2H3N1O1', 
                                      residues = 'C', 
                                      type = 'opt', 
                                      position = 'any')
mods(par)[[2]] <- msgfParModification(name = 'Oxidation', 
                                      mass = 15.994915, 
                                      residues = 'M', 
                                      type = 'opt', 
                                      position = 'any')
nMod(par) <- 2

# Run!
res <- runMSGF(par, 'Bait_DirectInfusion2.mzML')

########### Process results ###########
# Get peptides result:
resPep <- res@peptides@peptides
names(resPep)[names(resPep)=="id"] <- "peptide_ref"

# Get database result:
resDB <- res@database@database
names(resDB)[names(resDB) == "id"] <- "dBSequence_ref"

# Get evidence:
resEv <- res@evidence@evidence

# Get PSM result:
resPSM <- res@psm@id
names(resPSM)[names(resPSM) == "id"] <- "scanID"

# Joins:
resPSM <- resPSM %>% left_join(resPep)
resPSM <- resPSM %>% left_join(resEv, by = c("peptide_ref" = "peptide_ref"))
resPSM <- resPSM %>% left_join(resDB)

resPSMq10 <- resPSM %>% filter(`MS-GF:QValue` < .1)
resPSMq10 <- resPSMq10 %>% filter(!isDecoy)
writexl::write_xlsx(resPSMq10, path = "resPSMq10.xlsx")

