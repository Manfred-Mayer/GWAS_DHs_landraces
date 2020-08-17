###################################################################
###################################################################
####
#### 1) convert geno file of DH lines into gpData object (synbreed package)
####
#### Manfred Mayer (Technical University of Munich, Plant Breeding)
#### manfred.mayer@tum.de
####
#### date: 25.05.2020
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
### convert geno file of DH lines into gpData object (synbreed package)
###

# load genotypic data of DH lines (filtered and imputed)
genoData <- read.table("Mayer_et_al_genotypes_DHlines_600k_filteredImputed.txt.gz", header = TRUE, stringsAsFactors = FALSE)

# separate into map and geno
map <- genoData[, 2:3]
rownames(map) <- genoData$marker
str(map)
head(map)
geno <- genoData[, -(1:5)]
geno <- as.matrix(geno)
geno <- t(geno)
colnames(geno) <- rownames(map)
str(geno)

# generate gpData
fam <- substr(rownames(geno), 4, 5)
fam <- as.data.frame(fam, stringsAsFactors = FALSE)
rownames(fam) <- rownames(geno)
str(fam)

gpDH <- create.gpData(geno = geno, map = map, family = fam, map.unit = "bp")	
gpDH$info$codeGeno <- TRUE
save(gpDH, file = "gpDH.RData")
rm(geno, map, fam)

###
###################################################################
