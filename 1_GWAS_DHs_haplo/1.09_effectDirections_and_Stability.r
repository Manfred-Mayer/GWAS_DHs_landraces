###################################################################
###################################################################
####
#### for the final set of focus haplotypes,
#### distinguish into favorable, unfavorable, and interacting haplotypes
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

# arguments
nSNPs <- 10
steps <- 10
FDR <- "0.15"
p_thresh <- 0.01

###################################################################
###
### 
###

dir.create("FavUnfav_Stability")

# for all traits analyzed in the study
traits <- c("EV_V4", "EV_V6", "PH_V4", "PH_V6",
			"PH_final",
			"FF", "MF",
			"LO", "TILL")

n_fav_unfav_inter <- data.frame(TRAIT = traits,
								n_fav = rep(NA, length(traits)),
								n_unfav = rep(NA, length(traits)),
								n_inter = rep(NA, length(traits)))
rownames(n_fav_unfav_inter) <- traits
n_fav_unfav_inter

for(TRAIT in traits){
	# input folder
	infolder <- paste(TRAIT, "/QTLregs/finalQTL", sep = "")
	infolder

	# geno and map
	load(paste("geno_InfoList_DH_m", nSNPs, "s", steps, ".RData", sep = ""))
	str(InfoList)
	str(geno)

	dir.create(paste("FavUnfav_Stability/", TRAIT, sep = ""))
	outfolder <- paste("FavUnfav_Stability/", TRAIT, sep = "")
	outfolder

	Envs <- read.table(paste(TRAIT, "/", list.files(path = TRAIT, pattern = ".phenotypes_order.txt"), sep = ""), stringsAsFactors = FALSE)
	Envs <- Envs[,1]
	print(Envs)

	# load hap set final
	if(file.exists(paste(infolder, "/final_HapSet_", TRAIT, ".RData", sep =""))){
		load(paste(infolder, "//final_HapSet_", TRAIT, ".RData", sep =""))
		str(info_qtls_sig_95_final)

		HapSet_info_sig_95 <- info_qtls_sig_95_final
		HapSet_info_all <- info_qtls_all_final

		n_env_all <- length(Envs[which(Envs != "Across")])
		n_env_all

		# if in biallelic case, haplotypes were added, we have to do it here also for geno and InfoList
		if(any(unique(info_qtls_sig_95_final$Marker) %in% colnames(geno) == FALSE)){
			temp_x <- unique(info_qtls_sig_95_final$Marker)
			temp_x <- temp_x[which(temp_x %in% colnames(geno) == FALSE)]
			for(temp_xi in temp_x){
				temp_xix <- geno[ , which(substr(colnames(geno), 1, 13) == substr(temp_xi, 1, 13))]
				temp_xix[which(temp_xix == 0)] <- 99
				temp_xix[which(temp_xix == 2)] <- 0
				temp_xix[which(temp_xix == 99)] <- 2
				geno <- cbind(geno, temp_xix)
				colnames(geno)[ncol(geno)] <- temp_xi

				InfoList_xi <- InfoList[which(substr(rownames(InfoList), 1, 13) == substr(temp_xi, 1, 13)), ]
				InfoList <- rbind(InfoList, InfoList_xi)
				rownames(InfoList)[nrow(InfoList)] <- temp_xi
				print(tail(InfoList, n = 10))
			}
		}


		#
		# calculate for each haplotype
		# number of environments with
		# 	positive
		# 	negative
		# 	non-significant
		# effect
		#
		n_env_PosNegInfo <- NULL
		for(qtl_i in final_HapSet){
			
			# effects
			info_qtl <- HapSet_info_sig_95[which(HapSet_info_sig_95$Marker == qtl_i), ]
			info_qtl <- info_qtl[which((info_qtl$Qeff != 0) & (info_qtl$seQeff != 0)), ]
			if(nrow(info_qtl) > 0){
			# generate Qeff_sd vector
			Qeff_sd <- NULL
				for(env_i in unique(info_qtl$Environment)){
					Y_temp <- read.table(paste(TRAIT, "/", list.files(path = TRAIT, pattern = ".pheno.txt"), sep = ""), stringsAsFactors = FALSE)
					colnames(Y_temp) <- Envs
					Y_temp <- Y_temp[ , env_i]			
					Qeff_sd <- c(Qeff_sd, info_qtl$Qeff[which(info_qtl$Environment == env_i)] / sd(Y_temp, na.rm = TRUE))
					names(Qeff_sd)[length(Qeff_sd)] <- env_i
				}
			max_eff_sd <- Qeff_sd[which(Qeff_sd == max(Qeff_sd))[1]]
			min_eff_sd <- Qeff_sd[which(Qeff_sd == min(Qeff_sd))[1]]
			mean_eff_sd <- mean(Qeff_sd)
			max_diff <- abs(max_eff_sd - min_eff_sd)
			n_env_pos <- length(which(Qeff_sd > 0))
			n_env_neg <- length(which(Qeff_sd < 0))
			n_env_nonsig <- n_env_all - n_env_pos - n_env_neg
			
			info_temp <- data.frame(qtl_i, n_env_pos, n_env_neg, n_env_nonsig, max_eff_sd, min_eff_sd, mean_eff_sd, max_diff)
			
			} else {
			
			info_temp <- data.frame(qtl_i, 0, 0, n_env_all, NA, NA, NA, NA)
			colnames(info_temp) <- c("qtl_i", "n_env_pos", "n_env_neg", "n_env_nonsig", "max_eff_sd", "min_eff_sd", "mean_eff_sd", "max_diff")
			
			}
			
			n_env_PosNegInfo <- rbind(n_env_PosNegInfo, info_temp)

		}

		# define what is favorable/unfavorable per trait
			if(TRAIT %in% c('EV_V4', 'EV_V6','PH_V4', 'PH_V6',
							'PH_final')){
			
				colnames(n_env_PosNegInfo)[2] <- "n_env_fav"
				colnames(n_env_PosNegInfo)[3] <- "n_env_unfav"
				colnames(n_env_PosNegInfo)[4] <- "n_env_neut"
			
			} else {
			#FF, MF, LO, TILL

				colnames(n_env_PosNegInfo)[3] <- "n_env_fav"
				colnames(n_env_PosNegInfo)[2] <- "n_env_unfav"
				colnames(n_env_PosNegInfo)[4] <- "n_env_neut"
			
			}

		n_env_FavUnfavInfo <- n_env_PosNegInfo[ , c("qtl_i", "n_env_fav", "n_env_unfav", "n_env_neut", "max_eff_sd", "min_eff_sd", "mean_eff_sd", "max_diff")]
		n_env_FavUnfavInfo

		haps_Fav_info <- n_env_FavUnfavInfo[which((n_env_FavUnfavInfo[ , "n_env_unfav"] == 0) & (n_env_FavUnfavInfo[ , "n_env_fav"] != 0)), ]
		haps_Unfav_info <- n_env_FavUnfavInfo[which((n_env_FavUnfavInfo[ , "n_env_fav"] == 0) & (n_env_FavUnfavInfo[ , "n_env_unfav"] != 0)), ]
		haps_Inter_info <- n_env_FavUnfavInfo[which((n_env_FavUnfavInfo[ , "n_env_unfav"] != 0) & (n_env_FavUnfavInfo[ , "n_env_fav"] != 0)), ]

		info_nEnv_fav_unfav <-rbind(summary(haps_Fav_info[ , "n_env_fav"]),
									summary(haps_Unfav_info[ , "n_env_unfav"])
									)
		rownames(info_nEnv_fav_unfav) <- c("fav_haps", "unfav_haps")
		info_nEnv_fav_unfav

		write.table(haps_Fav_info, file = paste(outfolder, "/haps_Fav_info_", TRAIT , ".csv", sep = ""), sep = ";", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)
		write.table(haps_Unfav_info, file = paste(outfolder, "/haps_Unfav_info_", TRAIT , ".csv", sep = ""), sep = ";", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)
		write.table(haps_Inter_info, file = paste(outfolder, "/haps_Inter_info_", TRAIT , ".csv", sep = ""), sep = ";", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)
		write.table(info_nEnv_fav_unfav, file = paste(outfolder, "/info_nEnv_fav_unfav_", TRAIT , ".csv", sep = ""), sep = ";", dec = ".", col.names = NA, row.names = TRUE, quote = FALSE)
		save(haps_Fav_info, haps_Unfav_info, haps_Inter_info, file = paste(outfolder, "/haps_all_", TRAIT , ".RData", sep = ""))

		n_fav_unfav_inter[TRAIT, "n_fav"] <- nrow(haps_Fav_info)
		n_fav_unfav_inter[TRAIT, "n_unfav"] <- nrow(haps_Unfav_info)
		n_fav_unfav_inter[TRAIT, "n_inter"] <- nrow(haps_Inter_info)

	}
}

write.table(n_fav_unfav_inter, file = paste("FavUnfav_Stability/n_fav_unfav_Inter_allTraits.csv", sep = ""), sep = ";", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)

###
###################################################################
