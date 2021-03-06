###################################################################
###################################################################
####
#### calculate Polymorphism Information Content (PIC) according to Botstein et al. (1980)
#### for DHs of KE, LL, PE separately as well as overall
#### for BLs
####
#### use haplotype windows
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

# function (PIC)
# x = genotype matrix, with individuals in rows and markers (or haplotype windows) in columns
PIC_f <- function(x){
	PIC_per_marker <- apply(x, 2, PIC_per_matker_f)
	PIC_avg <- mean(PIC_per_marker)
	return(list(PIC_avg = PIC_avg,
				PIC_per_marker = PIC_per_marker))
}
# x = vector of alleles
PIC_per_matker_f <- function(x){
	allele_counts <- table(x)
	allele_freqs <- allele_counts/sum(allele_counts)
	alleles <- names(allele_freqs)
	rows <- cbind(alleles[sequence(1:length(alleles))],alleles[rep(1:length(alleles),1:length(alleles))])
	rows <- data.frame(	allele1 = rows[,1],
						allele2 = rows[,2])
	rows <- rows[which(!(rows[,1]==rows[,2])),]
	rows$p_square_allele1 <- allele_freqs[rows$allele1]^2
	rows$p_square_allele2 <- allele_freqs[rows$allele2]^2
	last_term_vec <- apply(rows[ ,3:4], 1, function(x){2 * x[1] * x[2]})
	PIC <- 1 - sum(allele_freqs^2) - sum(last_term_vec)
	return(PIC)
}

###################################################################
###
### load haplo data and merge chromosomes
###

haplo_all <- NULL

for(CHR in 1:10){
print(CHR)

load(paste("haplo_DH.BL_chr", CHR, "_m", nSNPs, "s", steps, ".RData", sep=""))
str(haplo)
haplo_all <- cbind(haplo_all, haplo)
rm(haplo)

}
str(haplo_all)
haplo <- haplo_all
rm(haplo_all)
save(haplo, file = paste("haplo_DH.BL_m", nSNPs, "s", steps, ".RData", sep = ""))
str(haplo)

###
###################################################################


###################################################################
###
### PIC calculation
###

# output folder
dir.create(paste("diversityParams", sep = ""))
dir.create(paste("diversityParams/PIC", sep = ""))
outfolder <- paste("diversityParams/PIC", sep = "")
outfolder

info_PIC_avg <- matrix(NA, nrow = 5, ncol = 1)
rownames(info_PIC_avg) <- c("KE", "LL", "PE", "DHall", "BL")
colnames(info_PIC_avg) <- "PIC"
info_PIC_avg

# calculate PIC for each landrace separately
for(LR in c("KE", "LL", "PE")){
print(LR)
geno <- haplo[which(substr(rownames(haplo), 1, 5) == paste0("DH_", LR)), ]
str(geno)
PIC <- PIC_f(geno)
save(PIC, file = paste(outfolder, "/PIC_", LR, "_haplo_m", nSNPs, "s", steps, ".RData", sep = ""))
info_PIC_avg[LR, "PIC"] <- PIC$PIC_avg
write.table(info_PIC_avg, file = paste0(outfolder, "/info_PIC_avg_haplo_m", nSNPs, "s", steps, ".csv"), sep = ";", dec = ".", quote = FALSE, col.names = FALSE, row.names = TRUE)
print(info_PIC_avg)
rm(geno, PIC)
}

# calculate PIC for all landraces together
geno <- haplo[which(substr(rownames(haplo), 1, 3) == paste0("DH_")), ]
str(geno)
PIC <- PIC_f(geno)
save(PIC, file = paste(outfolder, "/PIC_DHall_haplo_m", nSNPs, "s", steps, ".RData", sep = ""))
info_PIC_avg["DHall", "PIC"] <- PIC$PIC_avg
write.table(info_PIC_avg, file = paste0(outfolder, "/info_PIC_avg_haplo_m", nSNPs, "s", steps, ".csv"), sep = ";", dec = ".", quote = FALSE, col.names = FALSE, row.names = TRUE)
print(info_PIC_avg)
rm(geno, PIC)

# calculate PIC for BLs
geno <- haplo[which(substr(rownames(haplo), 1, 3) == paste0("EL_")), ]
str(geno)
PIC <- PIC_f(geno)
save(PIC, file = paste(outfolder, "/PIC_BLs_haplo_m", nSNPs, "s", steps, ".RData", sep = ""))
info_PIC_avg["BL", "PIC"] <- PIC$PIC_avg
write.table(info_PIC_avg, file = paste0(outfolder, "/info_PIC_avg_haplo_m", nSNPs, "s", steps, ".csv"), sep = ";", dec = ".", quote = FALSE, col.names = FALSE, row.names = TRUE)
print(info_PIC_avg)
rm(geno, PIC)

###
###################################################################
