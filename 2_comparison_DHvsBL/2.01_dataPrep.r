###################################################################
###################################################################
####
#### merge genotype files of DH lines and breeding line
####
#### Manfred Mayer (Technical University of Munich, Plant Breeding)
#### manfred.mayer@tum.de
####
#### date: 17.04.2020
###################################################################
###################################################################

# general settings
options(stringsAsFactors=FALSE)
options(scipen=99)
options(warn = 1)
set.seed(212)

# load packages
library(synbreed)

###################################################################
###
### merge genotype files of DH lines and breeding line
###

# load genotypic data of DH lines
load("gpDH.RData")

# load gpData of breeding lines (already filtered for the same marker set and imputed)
load("gpBL.RData")

# merge the two data files
all.equal(colnames(gpDH$geno), colnames(gpBL$geno))
geno <- rbind(gpDH$geno, gpBL$geno)
map <- gpDH$map
fam <- c(gpDH$covar$family, gpBL$covar$family)
fam <- as.data.frame(fam, stringsAsFactors = FALSE)
rownames(fam) <- rownames(geno)

gpDH.BL <- create.gpData(geno = geno, map = map, family = fam, map.unit = "bp")	
gpDH.BL$info$codeGeno <- TRUE
save(gpDH.BL, file = "gpDH.BL.RData")
rm(geno, map, fam)

###
###################################################################
