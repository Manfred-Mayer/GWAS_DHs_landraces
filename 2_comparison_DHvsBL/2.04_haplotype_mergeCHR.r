###################################################################
###################################################################
####
#### prepare haplotype data
#### merge chromosomes to single file
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

# arguments
nSNPs <- 10
steps <- 10

###################################################################
###
### merge chromosomes
###

geno <- NULL
InfoList_all <- NULL

for(CHR in 1:10){
print(CHR)

load(paste("haplo_binary_DH.BL_chr", CHR, "_m", nSNPs, "s", steps, ".RData", sep = ""))
geno <- cbind(geno, haplo_binary)
load(paste("InfoList_DH.BL_chr", CHR, "_m", nSNPs, "s", steps, ".RData", sep = ""))
InfoList_all <- rbind(InfoList_all, InfoList)
rm(haplo_binary, InfoList)

}
InfoList <- InfoList_all
rm(InfoList_all)

str(InfoList)
str(geno)

# recode into numeric
recode_num_f <- function(x){
	as.numeric(x)
}
names_geno <- rownames(geno)
geno <- apply(geno, 2, recode_num_f)
rownames(geno) <- names_geno
str(geno)
save(geno, InfoList, file = paste("geno_InfoList_DH.BL_m", nSNPs, "s", steps, ".RData", sep = ""))

###
###################################################################
