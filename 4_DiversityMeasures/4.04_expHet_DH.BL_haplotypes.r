###################################################################
###################################################################
####
#### calculate gene diversity (expected heterozygosity) according to Nei 1973
#### for DHs of KE, LL, PE separately as well as overall
#### for BLs
####
#### use haplotypes
####
#### Manfred Mayer (Technical University of Munich, Plant Breeding)
#### manfred.mayer@tum.de
####
#### date: 30.06.2020
###################################################################
###################################################################

# general settings
options(stringsAsFactors=FALSE)
options(scipen=99)
options(warn = 1)
set.seed(212)

# arguments
nSNPs <- 10
steps <- 10

# load packages
library(synbreed)

# function (expHet)
# x = genotype matrix, with individuals in rows and markers (or haplotype windows) in columns
expHet_f <- function(x){
	expHet_per_marker <- apply(x, 2, expHet_per_matker_f)
	expHet_avg <- mean(expHet_per_marker)
	return(list(expHet_avg = expHet_avg,
				expHet_per_marker = expHet_per_marker))
}
# x = vector of alleles
expHet_per_matker_f <- function(x){
	allele_counts <- table(x)
	allele_freqs <- allele_counts/sum(allele_counts)
	expHet <- 1 - sum(allele_freqs^2)
	return(expHet)
}

###################################################################
###

# output folder
dir.create(paste("diversityParams", sep = ""))
dir.create(paste("diversityParams/expHet", sep = ""))
outfolder <- paste("diversityParams/expHet", sep = "")
outfolder

info_expHet_avg <- matrix(NA, nrow = 5, ncol = 1)
rownames(info_expHet_avg) <- c("KE", "LL", "PE", "DHall", "BL")
colnames(info_expHet_avg) <- "expHet"
info_expHet_avg

# load genotypic data of DH and BL lines
load(paste("haplo_DH.BL_m", nSNPs, "s", steps, ".RData", sep = ""))
str(haplo)
haplo[1:10,1:3]

# calculate expHet for each landrace separately
for(LR in c("KE", "LL", "PE")){
print(LR)
geno <- haplo[which(substr(rownames(haplo), 1, 5) == paste0("DH_", LR)), ]
str(geno)
expHet <- expHet_f(geno)
save(expHet, file = paste(outfolder, "/expHet_", LR, "_haplo.RData", sep = ""))
info_expHet_avg[LR, "expHet"] <- expHet$expHet_avg
write.table(info_expHet_avg, file = paste0(outfolder, "/info_expHet_avg_haplo.csv"), sep = ";", dec = ".", quote = FALSE, col.names = FALSE, row.names = TRUE)
print(info_expHet_avg)
rm(geno, expHet)
}

# calculate expHet for all landraces together
geno <- haplo[which(substr(rownames(haplo), 1, 3) == paste0("DH_")), ]
str(geno)
expHet <- expHet_f(geno)
save(expHet, file = paste(outfolder, "/expHet_DHall_haplo.RData", sep = ""))
info_expHet_avg["DHall", "expHet"] <- expHet$expHet_avg
write.table(info_expHet_avg, file = paste0(outfolder, "/info_expHet_avg_haplo.csv"), sep = ";", dec = ".", quote = FALSE, col.names = FALSE, row.names = TRUE)
print(info_expHet_avg)
rm(geno, expHet)

# calculate expHet for BLs
geno <- haplo[which(substr(rownames(haplo), 1, 3) == paste0("EL_")), ]
str(geno)
expHet <- expHet_f(geno)
save(expHet, file = paste(outfolder, "/expHet_BLs_haplo.RData", sep = ""))
info_expHet_avg["BL", "expHet"] <- expHet$expHet_avg
write.table(info_expHet_avg, file = paste0(outfolder, "/info_expHet_avg_haplo.csv"), sep = ";", dec = ".", quote = FALSE, col.names = FALSE, row.names = TRUE)
print(info_expHet_avg)
rm(geno, expHet)

###
###################################################################
