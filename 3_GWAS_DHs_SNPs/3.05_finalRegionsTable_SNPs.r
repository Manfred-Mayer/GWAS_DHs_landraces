###################################################################
###################################################################
####
#### after the final set of trait-associated SNPs is defined,
#### now generate table of corresponding associated genomic regions
#### to be used for candidate gene search etc.
####
#### we go back to our already defined regions according to our GWAS results (based on LD; script 3.03)
#### and merge for each SNP the final list the corresponding regions of the environment-specific gwas
####
#### Manfred Mayer (Technical University of Munich, Plant Breeding)
#### manfred.mayer@tum.de
####
#### date: 03.07.2020
###################################################################
###################################################################

# general settings
options(stringsAsFactors=FALSE)
options(scipen=99)
options(warn = 1)
set.seed(212)

# load packages
library(asreml)

# arguments
TRAIT <- "PH_final"
FDR <- "0.15"
p_thresh <- 0.01

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
TRAIT

###################################################################
###
### 
###

# input folder
infolder <- paste("GWAS_SNPs/", TRAIT, "/QTLregs", sep = "")
infolder

# output folder
outfolder <- paste("GWAS_SNPs/", TRAIT, "/QTLregs/finalQTL", sep = "")
outfolder

# load data
# geno and map
load(paste("geno_InfoList_DH_SNPs.RData", sep = ""))
str(InfoList)
str(geno)

# load annotated gene information for B73v4
# it's a dataframe with columns: name, chr, start_bp, end_bp, strand
# and one row per gene
load(paste("gene_pos_Zea_mays.B73_RefGen_v4.43_FGS.RData", sep = ""))
str(gene_pos)
# consider only genes mapped on one of the 10 chromosomes
gene_pos <- gene_pos[which(gene_pos$chr %in% c(1,2,3,4,5,6,7,8,9,10)), ]
str(gene_pos)
# also include regions ...kb upstream of genes
extGenes <- 5
for(i in 1: nrow(gene_pos)){
	if(gene_pos[i, 5] == "+"){
		gene_pos[i, 3] <- gene_pos[i, 3] - (extGenes * 1000)
	}
	if(gene_pos[i, 5] == "-"){
		gene_pos[i, 4] <- gene_pos[i, 4] + (extGenes * 1000)
	}
}
gene_pos$chr <- as.numeric(gene_pos$chr)
names(gene_pos)[1] <- "Gene_ID"
str(gene_pos)
head(gene_pos)

# environments
Envs <- read.table(paste("GWAS_SNPs/", TRAIT, "/", list.files(path = paste0("GWAS_SNPs/", TRAIT), pattern = ".phenotypes_order.txt"), sep = ""), stringsAsFactors = FALSE)
Envs <- Envs[,1]
print(Envs)

load(paste(outfolder, "/OutFinalQTLset_", TRAIT , "_FDR", FDR, "_p", p_thresh, ".RData", sep =""))
str(info_qtls_sig_95)
# how many haplotypes from the initial candidate set made it into the final set
str(candQTLset)
str(finalQTLset)
length(finalQTLset) / length(candQTLset)

# load regions defined within envs (candQTLset)
reg_all <- NULL
	for(env_i in Envs){
		if(file.exists(paste(infolder, "/reg_list_merge_", TRAIT , ".", env_i, "_FDR", FDR, ".txt", sep =""))){
			res <- read.table(paste(infolder, "/reg_list_merge_", TRAIT , ".", env_i, "_FDR", FDR, ".txt", sep =""), header = TRUE)
			env <- rep(env_i, nrow(res))
			res <- cbind(res, env)
			reg_all <- rbind(reg_all, res)
			rm(res, env)
		}
	}
reg_all <- reg_all[order(reg_all$Chr, reg_all$Start_bp, reg_all$End_bp), ]
str(reg_all)
head(reg_all)

# filter for finalQTLset
reg_finalQTL <- reg_all[which(reg_all$mostSigMarker %in% finalQTLset), ]
str(reg_finalQTL)
head(reg_finalQTL)

#
# generate QTL reg list, using originally defined regions for the QTLs of the finalQTLset
# if QTL is significant in multiple environments, take outer borders to define region
#
finalQTLregs <- NULL
gene_list <- list()
for(qtl_i in finalQTLset){
	# qtl reg
	chr_qtl <- reg_finalQTL[which(reg_finalQTL$mostSigMarker == qtl_i), "Chr"][1]
	start_qtl <- min(reg_finalQTL[which(reg_finalQTL$mostSigMarker == qtl_i), "Start_bp"])
	end_qtl <- max(reg_finalQTL[which(reg_finalQTL$mostSigMarker == qtl_i), "End_bp"])
	size_qtl_kb <- (end_qtl - start_qtl + 1) / 1000
	
	# effects
	info_qtl <- info_qtls_sig_95[which(info_qtls_sig_95$Marker == qtl_i), ]
	info_qtl <- info_qtl[which((info_qtl$Qeff != 0) & (info_qtl$seQeff != 0)), ]
	# generate Qeff_sd vector
	Qeff_sd <- NULL
		for(env_i in unique(info_qtl$Environment)){
			Y_temp <- read.table(paste0("GWAS_SNPs/", TRAIT, "/", list.files(path = paste0("GWAS_SNPs/", TRAIT), pattern = ".pheno.txt")), stringsAsFactors = FALSE)
			colnames(Y_temp) <- Envs
			Y_temp <- Y_temp[ , env_i]			
			Qeff_sd <- c(Qeff_sd, info_qtl$Qeff[which(info_qtl$Environment == env_i)] / sd(Y_temp, na.rm = TRUE))
			names(Qeff_sd)[length(Qeff_sd)] <- env_i
		}
	
	# just in very rare cases, the haplotype might not be significant for CI95% in any environment, then go down to CI90%
	if(is.null(Qeff_sd)){
		info_qtl <- info_qtls[which(info_qtls$Marker == qtl_i), ]
		lowerCI_90 <- info_qtl$Qeff - qnorm(0.95, mean = 0, sd = 1) * info_qtl$seQeff
		upperCI_90 <- info_qtl$Qeff + qnorm(0.95, mean = 0, sd = 1) * info_qtl$seQeff
		sig_f <- function(x){
			(findInterval(0, x) == 1) == FALSE
		}		
		info_qtl <- info_qtl[apply(cbind(lowerCI_90, upperCI_90), 1, sig_f), ]
		for(env_i in unique(info_qtl$Environment)){
			Y_temp <- read.table(paste0("GWAS_SNPs/", TRAIT, "/", list.files(path = paste0("GWAS_SNPs/", TRAIT), ".pheno.txt")), stringsAsFactors = FALSE)
			colnames(Y_temp) <- Envs
			Y_temp <- Y_temp[ , env_i]			
			Qeff_sd <- c(Qeff_sd, info_qtl$Qeff[which(info_qtl$Environment == env_i)] / sd(Y_temp, na.rm = TRUE))
			names(Qeff_sd)[length(Qeff_sd)] <- env_i
		}

		if(is.null(Qeff_sd) == FALSE){
			max_eff_sd <- Qeff_sd[which(Qeff_sd == max(Qeff_sd))[1]]
			min_eff_sd <- Qeff_sd[which(Qeff_sd == min(Qeff_sd))[1]]
			mean_eff_sd <- mean(Qeff_sd)
			max_env_sd <- paste(names(Qeff_sd)[which(Qeff_sd == max_eff_sd)], "_CI90", sep = "")
			min_env_sd <- paste(names(Qeff_sd)[which(Qeff_sd == min_eff_sd)], "_CI90", sep = "")	
			# absolute effects
			max_eff_abs <-info_qtl$Qeff[which(info_qtl$Qeff == max(info_qtl$Qeff))[1]]
			min_eff_abs <- info_qtl$Qeff[which(info_qtl$Qeff == min(info_qtl$Qeff))[1]]
			mean_eff_abs <- mean(info_qtl$Qeff)
			max_env_abs <- info_qtl$Environment[which(info_qtl$Qeff == max_eff_abs)]
			min_env_abs <- info_qtl$Environment[which(info_qtl$Qeff == min_eff_abs)]
			# number of environments with significant positive/negative effect
			n_env <- length(Qeff_sd)		
		}

	} else {	
	max_eff_sd <- Qeff_sd[which(Qeff_sd == max(Qeff_sd))[1]]
	min_eff_sd <- Qeff_sd[which(Qeff_sd == min(Qeff_sd))[1]]
	mean_eff_sd <- mean(Qeff_sd)
	max_env_sd <- names(Qeff_sd)[which(Qeff_sd == max_eff_sd)]
	min_env_sd <- names(Qeff_sd)[which(Qeff_sd == min_eff_sd)]	
	# absolute effects
	max_eff_abs <- info_qtl$Qeff[which(info_qtl$Qeff == max(info_qtl$Qeff))[1]]
	min_eff_abs <- info_qtl$Qeff[which(info_qtl$Qeff == min(info_qtl$Qeff))[1]]
	mean_eff_abs <- mean(info_qtl$Qeff)
	max_env_abs <- info_qtl$Environment[which(info_qtl$Qeff == max_eff_abs)]
	min_env_abs <- info_qtl$Environment[which(info_qtl$Qeff == min_eff_abs)]
	# number of environments with significant positive/negative effect
	n_env <- length(Qeff_sd)
	}

	if(is.null(Qeff_sd) == FALSE){
	# how many annotated genes lie within that region
	geneIDs <- gene_pos$Gene_ID[which((gene_pos$chr == chr_qtl) & (gene_pos$start_bp <= end_qtl) & (gene_pos$end_bp >= start_qtl))]
	if(length(geneIDs) > 0){
	gene_list[[qtl_i]] <- geneIDs
	} else {
	gene_list[[qtl_i]] <- NA
	}
	n_genes <- length(geneIDs)
	
	# region of marker
	start_hap <- InfoList[qtl_i, 2]
	end_hap <- InfoList[qtl_i, 3]
	
	info_temp <- data.frame(qtl_i, chr_qtl, start_hap, end_hap, start_qtl, end_qtl, size_qtl_kb, min_eff_sd, max_eff_sd, mean_eff_sd, min_env_sd, max_env_sd, n_env, n_genes, min_eff_abs, max_eff_abs, mean_eff_abs)
	
	} else {
	# how many annotated genes lie within that region
	geneIDs <- gene_pos$Gene_ID[which((gene_pos$chr == chr_qtl) & (gene_pos$start_bp <= end_qtl) & (gene_pos$end_bp >= start_qtl))]
	if(length(geneIDs) > 0){
	gene_list[[qtl_i]] <- geneIDs
	} else {
	gene_list[[qtl_i]] <- NA
	}
	n_genes <- length(geneIDs)

	# region of marker
	start_hap <- InfoList[qtl_i, 2]
	end_hap <- InfoList[qtl_i, 3]

	info_temp <- data.frame(qtl_i, chr_qtl, start_hap, end_hap, start_qtl, end_qtl, size_qtl_kb, NA, NA, NA, NA, NA, 0, n_genes, NA, NA, NA)
	colnames(info_temp) <- c("qtl_i", "chr_qtl", "start_hap", "end_hap", "start_qtl", "end_qtl", "size_qtl_kb", "min_eff_sd", "max_eff_sd", "mean_eff_sd",
							 "min_env_sd", "max_env_sd", "n_env", "n_genes", "min_eff_abs", "max_eff_abs", "mean_eff_abs")
	
	}
	finalQTLregs <- rbind(finalQTLregs, info_temp)
}
write.table(finalQTLregs, file = paste(outfolder, "/finalQTLregs_", TRAIT , ".csv", sep = ""), sep = ";", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)
finalQTLregs
save(gene_list, file = paste(outfolder, "/gene_list_", TRAIT , ".RData", sep = ""))
gene_list

###
###################################################################
