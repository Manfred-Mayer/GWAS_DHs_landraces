###################################################################
###################################################################
####
#### generate GWAS input files
#### 1) merge chromosomes of haplotype dataset into one object
#### 2) calculate kinship matrix based on SNPs
#### 3) prepare BLUEs within environments for each trait
#### 4) prepare input files for gemma
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

# arguments
nSNPs <- 10
steps <- 10
minCount <- 3

###################################################################
###
### 1) merge chromosomes of haplotype dataset into one object
###

geno <- NULL
InfoList_all <- NULL

for(CHR in 1:10){
print(CHR)

load(paste("haplo_binary_DH_chr", CHR, "_m", nSNPs, "s", steps, ".RData", sep = ""))
geno <- cbind(geno, haplo_binary)
load(paste("InfoList_DH_chr", CHR, "_m", nSNPs, "s", steps, ".RData", sep = ""))
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
save(geno, InfoList, file = paste("geno_InfoList_DH_m", nSNPs, "s", steps, ".RData", sep = ""))

###
###################################################################


###################################################################
###
### 2) calculate kinship matrix based on SNPs
###

# load marker data
load("gpDH.RData")
str(gpDH)

# filter for SNPs with the minor allele being present at least minCount times
count_f <- function(x){
        y <- length(which(x==2))
        return(y)
}
count <- apply(gpDH$geno, 2, count_f)
summary(count)
# for reference as well as for alternative allele
ex_m <- names(count)[which(count <= minCount | count >= (nrow(gpDH$geno) - minCount))]
length(ex_m)

# generate filtered gp object(count)
gp <- discard.markers(gpDH, ex_m)
save(gp, file = paste("gp_DH_SNPs_count", minCount, ".RData", sep = ""))
str(gp)

###
### calculate kinship matrix
###

# function for recoding geno matrix according to (0 = homozygous for the major allele, 1 = heterozygous, 2 = homozygous for the minor allele)
recodeMAF_f <- function(x) {
	if((length(which(x == 0))*2 + length(which(x == 1))) < length(x)){
	  x[which(x == 0)] <- 99
	  x[which(x == 2)] <- 0
	  x[which(x == 99)] <- 2
	}
	return(x)
}

# edit allele coding in gp$geno
gp_kin <- gp
gp_kin$geno <- apply(gp_kin$geno, 2, recodeMAF_f)

# Astle and Balding (2009)
Kinship <- kin(gp_kin, ret="realizedAB")
save(Kinship, file = paste("Kin_AB_DHs_SNPs.RData", sep = ""))

###
###################################################################



###################################################################
###
### 3) prepare BLUEs within environment for each trait
###

# load phenotypic data
# BLUEs across environment
pheno_accEnv <- read.table("Mayer_et_al_phenotypes_BLUEs_acrossAllEnvironments.txt.gz", header = TRUE, stringsAsFactors = FALSE)
str(pheno_accEnv)
# remove all non-DHs from the pheno file
pheno_accEnv <- pheno_accEnv[which(substr(pheno_accEnv$Genotype, 1, 2) == "DH"), ]
rownames(pheno_accEnv) <- pheno_accEnv$Genotype

# BLUEs within each environment
pheno_perEnv <- read.table("Mayer_et_al_phenotypes_BLUEs_perEnvironment.txt.gz", header = TRUE, stringsAsFactors = FALSE)
str(pheno_perEnv)
# also here, remove all non-DHs from the pheno file
pheno_perEnv <- pheno_perEnv[which(substr(pheno_perEnv$Genotype, 1, 2) == "DH"), ]

# convert it into a list with tables for each trait
traits <- colnames(pheno_perEnv)[3:ncol(pheno_perEnv)]
traits

# generate template table for GWAS input
genotypes <- unique(pheno_perEnv$Genotype)
genotypes <- sort(genotypes)

environments <- unique(pheno_perEnv$Environment)
environments
environments_matrix <- matrix(as.numeric(NA), ncol = length(environments), nrow = length(genotypes))
colnames(environments_matrix) <- environments

GWAS_table_template <- data.frame(Genotype = genotypes)
GWAS_table_template <- cbind(GWAS_table_template, environments_matrix)
str(GWAS_table_template)
head(GWAS_table_template)

pheno_perEnv_list <- list()
for(TRAIT in traits){
  BLUEs_trait <- pheno_perEnv[, c("Genotype", "Environment", TRAIT)]
  GWAS_table <- GWAS_table_template
  rownames(GWAS_table) <- GWAS_table$Genotype
  for(ENV in environments){
    BLUEs_trait_env <- BLUEs_trait[which(BLUEs_trait$Env == ENV),]
    rownames(BLUEs_trait_env) <- BLUEs_trait_env$Genotype
    for(GENO in BLUEs_trait_env$Genotype){
      GWAS_table[GENO, ENV] <- BLUEs_trait_env[GENO, 3]
    }
  }
  colnames(GWAS_table) <- c("Genotype", colnames(GWAS_table)[-1])
  # remove columns with only NAs (trait not measured in this environment)
  rm_na_f <- function(x){
    any(is.na(x) == FALSE)
  }
  GWAS_table <- GWAS_table[, c(1, which(apply(GWAS_table[ , -1], 2, rm_na_f)) + 1)]
  
  pheno_perEnv_list[[TRAIT]] <- GWAS_table
}
str(pheno_perEnv_list)
save(pheno_perEnv_list, file = "pheno_perEnv_list_DHs.RData")

###
###################################################################



###################################################################
###
### 4) prepare input files for gemma
###

# while writing out input files also write gemma command lines
# command file (sh)
i <- 0
Script <- list()
## function for writing start script
write_file_f <- function(x, name) {
  write(x, append = T, file = name, ncolumns = length(x))
}

for(TRAIT in traits){
print(TRAIT)
pheno_trait <- pheno_perEnv_list[[TRAIT]]

# make directory for each trait
dir.create(TRAIT)
outDir <- TRAIT

# filter pheno_trait for individuals
# single environment
pheno_trait <- pheno_trait[which(pheno_trait$Genotype %in% rownames(geno)), ]
rownames(pheno_trait) <- pheno_trait$Genotype
# across environment
BLUEsAcross_trait <- pheno_accEnv[rownames(pheno_trait), ]
pheno_trait_across <- BLUEsAcross_trait[, TRAIT]
pheno_trait$Across <- pheno_trait_across

# filter geno for overlapping individuals
geno_trait <- geno[which(rownames(geno) %in% pheno_trait$Genotype), ]
print(all.equal(rownames(pheno_trait), rownames(geno_trait)))

# filter kinship for overlapping individuals
load(paste("Kin_AB_DHs_SNPs.RData", sep = ""))
Kinship <- Kinship[rownames(geno_trait), rownames(geno_trait)]
print(all.equal(rownames(Kinship), rownames(geno_trait)))
print(all.equal(colnames(Kinship), rownames(geno_trait)))

##
## write GEMMA input files
##
#
# marker file
#
M <- ncol(geno_trait)
chr <- InfoList[,1]
snpID <- rownames(InfoList)
bpPos <- apply(InfoList[ , 2:3], 1, function(x) {round(mean(x))})

#bimbam genotype file
prefix <- "DHs"
bimbam_geno <- data.frame(snpID, rep("A", length(snpID)), rep("G", length(snpID)), t(geno_trait))
write.table(bimbam_geno, file = file.path(outDir, paste(prefix, ".geno.txt", sep = "")), quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
system(paste("gzip ", file.path(outDir, paste(prefix, ".geno.txt", sep = "")), sep =""))
write.table(rownames(geno_trait), file = file.path(outDir, paste(prefix, ".geno_order.txt", sep = "")), quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)

#bimbam annotation file
all.equal(snpID, colnames(geno_trait))
bimbam_anno <- data.frame(snpID, bpPos, chr)
write.table(bimbam_anno, file = file.path(outDir, paste(prefix, ".anno.txt", sep = "")), quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)

#
# phenotype file
#

pheno <- pheno_trait[, -1]
print(all.equal(rownames(pheno), rownames(geno_trait)))

#bimbam phenotype file
bimbam_pheno <- data.frame(pheno)
write.table(bimbam_pheno, file = file.path(outDir, paste(prefix, ".pheno.txt", sep = "")), quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
write.table(colnames(pheno_trait)[-1], file = file.path(outDir, paste(prefix, ".phenotypes_order.txt", sep = "")), quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)

# prepare genomic relationship matrix
Kinship <- as.data.frame(Kinship)
write.table(Kinship, file = file.path(outDir, paste(prefix, ".Kinship.txt", sep = "")), quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
system(paste("gzip ", file.path(outDir, paste(prefix, ".Kinship.txt", sep = "")), sep =""))

# prepare landrace covariates
intercept <- rep(1, nrow(geno_trait))
KE <- rep(0, nrow(geno_trait))
KE[which(substr(rownames(geno_trait), 1, 5) == "DH_KE")] <- 1
LL <- rep(0, nrow(geno_trait))
LL[which(substr(rownames(geno_trait), 1, 5) == "DH_LL")] <- 1
LR_cov <- data.frame(intercept = intercept, KE = KE, LL = LL)
print(str(LR_cov))
write.table(LR_cov, file = file.path(outDir, paste(prefix, ".covariables.txt", sep = "")), quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)

	# each environment has different numbers of missing values for the phenotypes,
	# thus calculate maf threshold per environment according to minCount
	# also calculate set of "non-collinear" markers (r2 < 1 per chr)
	maf_thresh <- NULL
	nonColl_list <- list()
	for(env_i in colnames(pheno_trait)[-1]){
		geno_i <- pheno_trait$Geno[which(is.na(pheno_trait[, env_i]) == FALSE)]
		thresh_i <- minCount / (length(geno_i)-0.1)
		maf_thresh <- c(maf_thresh, thresh_i)
		duplicated_markers <- NULL
		for(CHR in 1:10){
			geno_chr <- geno[geno_i, which(InfoList[,1] == CHR)]
			duplicated_markers <- c(duplicated_markers, colnames(geno_chr)[which(duplicated(geno_chr, MARGIN = 2))])
		}
		nonColl_list[[env_i]] <- colnames(geno)[which((colnames(geno) %in% duplicated_markers) == FALSE)]		
	}
	names(maf_thresh) <- colnames(pheno_trait)[-1]
	print(maf_thresh)
	print(str(nonColl_list))
	save(nonColl_list, file = file.path(outDir, paste(prefix, ".nonColl_list.RData", sep = "")))
	
##
# write command lines, bash script
##

##
## for starting GWAS
##
	i <- i + 1
	Script[[i]] <- paste("", sep = "")
	i <- i + 1
	Script[[i]] <- paste("### ", TRAIT, " ###", sep = "")
	i <- i + 1
	Script[[i]] <- paste("cd ", outDir, sep = "")

for(env_i in colnames(pheno_trait)[-1]){
	index <- which(colnames(pheno_trait)[-1] == env_i)
	i <- i + 1
	Script[[i]] <- paste("gemma -g ", prefix, ".geno.txt.gz -p ", prefix, ".pheno.txt -a ", prefix, ".anno.txt -k ", prefix, ".Kinship.txt.gz -c ", prefix, ".covariables.txt -n ", index, " -lmm 4 -o ", TRAIT, "_", env_i, " -maf ", maf_thresh[env_i], " & disown -h", sep = "")
}


}
lapply(Script, write_file_f, name = paste("GEMMA_commands.sh", sep =""))


###
###################################################################
