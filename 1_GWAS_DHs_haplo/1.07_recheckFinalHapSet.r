###################################################################
###################################################################
####
#### in the backward selection (script 1.06), we already defined a
#### (final) set of focus haplotypes. However, theoretically it would be
#### possible that within the same window of such a haplotype, also an additional
#### haplotype was significant in GWAS but did not pass the backward elimination as
#### it is collinear with the selected focus haplotype
####
#### here, we check for such situations and if applicable add such alternative haplotypes
#### to the final list of haplotypes
#### (which then are checked for being favorable/unfavorable/interacting and compared with the breeding lines)
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
library(asreml)

# arguments
nSNPs <- 10
steps <- 10
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
infolder <- paste(TRAIT, "/QTLregs/finalQTL", sep = "")
infolder

# output folder
outfolder <- paste(TRAIT, "/QTLregs/finalQTL", sep = "")
outfolder

# load data
# geno and map
load(paste("geno_InfoList_DH_m", nSNPs, "s", steps, ".RData", sep = ""))
str(InfoList)
str(geno)

# load qtls of TRAIT
load(paste(infolder, "/OutFinalQTLset_", TRAIT , "_FDR", FDR, "_p", p_thresh, ".RData", sep =""))

#
# check initial GWAS results
#

# generate list of haplotypes which passed FDR15% in any GWAS scenario

	Envs <- read.table(paste(TRAIT, "/", list.files(path = TRAIT, pattern = ".phenotypes_order.txt"), sep = ""), stringsAsFactors = FALSE)
	Envs <- Envs[,1]
	print(Envs)

	reg_all <- NULL
		for(env_i in Envs){
			if(file.exists(paste(TRAIT, "/summary_GWAS/Hits_", TRAIT , "_", env_i, "_FDR", FDR, ".csv", sep =""))){
				res <- read.table(paste(TRAIT, "/summary_GWAS/Hits_", TRAIT , "_", env_i, "_FDR", FDR, ".csv", sep =""), header = TRUE, sep = ";", row.names = 1)
				if(is.na(res[1,1]) == FALSE){
					res$mk <- rownames(res)
					res$env <- rep(env_i, nrow(res))
					reg_all <- rbind(reg_all, res)
				}
				rm(res)
			}
		}
	print(str(reg_all))
	# vector of markers which were at least one time significant in any of the environment-specific GWAS runs
	all_sig_mk <- unique(reg_all$mk)
	n_haps <- length(all_sig_mk)
	# the corresponding genomic windows
	n_windows <- length(unique(substr(all_sig_mk, 1, 13)))
	print(paste("n_haps", n_haps))
	print(paste("n_windows", n_windows))

	# the haplotypes of the initial Q should all be included in this GWAS FDR15% list
	Q_markers <- unique(info_qtls_sig_95$Marker)
	print(paste("haps in Q passing FDR15% in GWAS:", length(which(Q_markers %in% all_sig_mk)), "/", length(Q_markers)))

# now, check if there are additional haplotypes in windows where a focus haplotype was selected
# windows containing a focus haplotype
windows_with_focus_haplotype <- substr(finalQTLset, 1, 13)
# haplotypes which passed FDR15% in GWAS and are within such a window
sig_haps <- all_sig_mk[which(substr(all_sig_mk, 1, 13) %in% windows_with_focus_haplotype)]
# potentially additional haplotypes to consider
add_haps <- sig_haps[which(sig_haps %in% finalQTLset == FALSE)]

add_haps_sig <- NULL
info_add_haps <- NULL
info_add_haps_sig95 <- NULL
if(length(add_haps) > 0){

	# test the potentially new haplotypes in the multi-locus, multi-environment model
	# by substituting the focus haplotype of the respective window with the alternative haplotype (the rest stays the same)
	
	# regenerate the input dataframe
	QTL_cand <- c(finalQTLset, add_haps)
	Y <- NULL
	ENV <- NULL
	GROUP <- NULL
	GENOTYPE <- NULL
	QTLs <- NULL
		# combine all environments (within-environment BLUEs), but exclude across-environment BLUEs
		for(env_i in Envs[which(Envs != "Across")]){
			Y_temp <- read.table(paste0(TRAIT, "/", list.files(path = TRAIT, pattern = ".pheno.txt")), stringsAsFactors = FALSE)
			colnames(Y_temp) <- Envs
			Y_temp <- Y_temp[ , env_i]

			GENOTYPE_temp <- read.table(paste0(TRAIT, "/", list.files(path = TRAIT, pattern = "geno_order.txt")), stringsAsFactors = FALSE)
			GENOTYPE_temp <- GENOTYPE_temp[,1]
			GENOTYPE_temp <- GENOTYPE_temp[which(is.na(Y_temp) == FALSE)]

			geno_gwas <- geno[GENOTYPE_temp, ]
			
			GROUP_temp <- read.table(paste0(TRAIT, "/", list.files(path = TRAIT, pattern = ".covariables.txt")), stringsAsFactors = FALSE)
			GROUP_temp <- GROUP_temp[,-1]
			GROUP_temp <- apply(GROUP_temp, 1, function(x){
				if(x[1] == 1){
					res <- "KE"
				}
				if(x[2] == 1){
					res <- "LL"
				}
				if((x[1] == 0) & (x[2] == 0)){
					res <- "PE"
				}
				return(res)
				}
				)
			GROUP_temp <- GROUP_temp[which(is.na(Y_temp) == FALSE)]
			
			Y_temp <- Y_temp[which(is.na(Y_temp) == FALSE)]
			
			ENV_temp <- rep(env_i, length(Y_temp))
			
			QTLs_temp <- data.frame(geno_gwas[, QTL_cand])
			
			Y <- c(Y, Y_temp)
			ENV <- c(ENV, ENV_temp)
			GENOTYPE <- c(GENOTYPE, GENOTYPE_temp)
			GROUP <- c(GROUP, GROUP_temp)
			QTLs <- rbind(QTLs, QTLs_temp)
		}
	input_data <- data.frame(Y, ENV, GROUP, GENOTYPE, QTLs)
	input_data$ENV <- as.factor(input_data$ENV)
	input_data$GROUP <- as.factor(input_data$GROUP)
	input_data$GENOTYPE <- as.factor(input_data$GENOTYPE)
	print(str(input_data))

	# bring focus haplotypes into the right order
	focusHaps <- final_model$fixed.formula
	focusHaps <- strsplit(as.character(focusHaps)[3], split = "+", fixed = TRUE)
	focusHaps <- gsub(" ", "", focusHaps[[1]])
	focusHaps <- gsub("\n", "", focusHaps)
	focusHaps <- gsub(")", "", focusHaps)
	focusHaps <- focusHaps[-(1:2)]
	focusHaps[1] <- substr(focusHaps[1], 6, nchar(focusHaps[1]))
	
	# start values
		test_model <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
							ENV:(", paste0(focusHaps, collapse = " + "),")", collapse = "")),
			   random = ~ GENOTYPE,
			   rcov   = ~ at(ENV):units,
			   data   = input_data,
			   maxit  = 500,
			   workspace   = 3e+08,
			   start.values=TRUE)
		iv <- test_model$gammas.table

	# now test each of the add_haps for significance
	for(m_i in add_haps){
			# swap haplotypes and put tested haplotype at the end
			wind.m_i <- substr(m_i, 1, 13)
			focusHaps.m_i <- focusHaps[-which(substr(focusHaps, 1, 13) == wind.m_i)]
			focusHaps.m_i <- c(focusHaps.m_i, m_i)
			
			# run model
			test_model <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
							ENV:(", paste0(focusHaps.m_i, collapse = " + "),")", collapse = "")),
			   random = ~ GENOTYPE,
			   rcov   = ~ at(ENV):units,
			   data   = input_data,
			   maxit  = 500,
			   workspace   = 3e+08,
			   G.param=iv, 
			   R.param=iv)
			
			# test significance of m_i
			wald.obj <- wald(test_model, denDF = "default")
			p_wald <- wald.obj$Wald$Pr
			names(p_wald) <- rownames(wald.obj$Wald)
			p_wald <- p_wald[which(substr(names(p_wald), 1, 4) == "ENV:")]
			names(p_wald) <- gsub("ENV:", "", names(p_wald))
			p_value.m_i <- p_wald[m_i]
		
			if(p_value.m_i < p_thresh){				
				# extract effect estimates
				Qeff <- test_model$coefficients$fixed[grep(":wind",names(test_model$coefficients$fixed))]
				seQeff <- sqrt(test_model$vcoeff$fixed)[1:length(Qeff)]
				
				# "Confidence intervals" (95%)
				lowerCI95 <- Qeff - qnorm(0.975, mean = 0, sd = 1) * seQeff
				upperCI95 <- Qeff + qnorm(0.975, mean = 0, sd = 1) * seQeff

				# summary table
				EnvMarker <- strsplit(names(Qeff), split = ":", fixed = TRUE)
				Marker <- unlist(lapply(EnvMarker, function(x){x[2]}))
				Environment <- unlist(lapply(EnvMarker, function(x){x[1]}))
				Environment <- gsub("ENV_", "", Environment)
				p_wald <- p_wald[Marker]
				info_qtls_m_i <- data.frame(Marker, Environment, p_wald, Qeff, seQeff, lowerCI95, upperCI95)
				info_qtls_m_i <- info_qtls_m_i[which(info_qtls_m_i$Marker == m_i), ]

				# add position information
				info_qtls_m_i$Chr <- InfoList[info_qtls_m_i$Marker, 1]
				info_qtls_m_i$Start_bp <- InfoList[info_qtls_m_i$Marker, 2]
				info_qtls_m_i$End_bp <- InfoList[info_qtls_m_i$Marker, 3]
				info_qtls_m_i$Pos_bp <- apply(InfoList[info_qtls_m_i$Marker, 2:3], 1, mean)
				info_qtls_m_i <- info_qtls_m_i[order(info_qtls_m_i$Chr, info_qtls_m_i$Start_bp, info_qtls_m_i$End_bp, info_qtls_m_i$Environment), ]

				# select significant effects (environment-specific)
				sig_f <- function(x){
					(findInterval(0, x) == 1) == FALSE
				}		
				sig_m <- apply(info_qtls_m_i[, c("lowerCI95", "upperCI95")], 1, sig_f)
				info_qtls_m_i_sig95 <- info_qtls_m_i[which(sig_m), ]

				# put in add object
				add_haps_sig <- c(add_haps_sig, m_i)
				info_add_haps <- rbind(info_add_haps, info_qtls_m_i)
				info_add_haps_sig95 <- rbind(info_add_haps_sig95, info_qtls_m_i)
			}
	}
	print(paste(add_haps_sig, " added"))
	
	# update final haplotype set
	final_HapSet <- c(focusHaps, add_haps_sig)

} else {

	final_HapSet <- finalQTLset

}


#
# for bi-allelic windows, only one haplotype was tested in GWAS, but of course the other haplotype would be significant as well with opposite effect sign
# thus, for considering favorable/unfavorable haplotypes, the respectiv alternativ haplotype has to be included as well
# (this happens only in very rare cases)
#

# screen for windows with only one hap in the gWAS results
biallelicWindows <- NULL
		for(env_i in Envs){
			if(file.exists(paste(TRAIT, "/output/", TRAIT , "_", env_i, ".assoc.txt.gz", sep =""))){
				res <- read.table(paste(TRAIT, "/output/", TRAIT , "_", env_i, ".assoc.txt.gz", sep =""), header = TRUE, stringsAsFactors = FALSE)
				winds <- substr(res$rs, 1, 13)
				n_alleles_per_wind <- table(winds)
				bi_allelic <- names(n_alleles_per_wind)[which(n_alleles_per_wind == 1)]
				biallelicWindows <- c(biallelicWindows, bi_allelic)
				rm(res)
			}
		}		
	biallelicWindows <- unique(biallelicWindows)
	print(str(biallelicWindows))

# are some of these bi-allelic windows in the final_HapSet? if yes, add the alternative haplotype
if(any(substr(final_HapSet, 1, 13) %in% biallelicWindows)){
	add_bi <- final_HapSet[which(substr(final_HapSet, 1, 13) %in% biallelicWindows)]
	for(add_bi_i in add_bi){
		add_table <- info_qtls_sig_95[which(info_qtls_sig_95$Marker == add_bi_i), ]
		add_table$Marker <- substr(add_table$Marker, 1, 14)
		add_table$Marker <- paste(add_table$Marker, "2", sep = "")
		
		add_table$Qeff <- -add_table$Qeff
		add_table$lowerCI95 <- -add_table$upperCI95
		add_table$upperCI95 <- -add_table$lowerCI95
		
		info_qtls_sig_95 <- rbind(info_qtls_sig_95, add_table)
		final_HapSet <- c(final_HapSet, add_table$Marker[1])
		rm(add_table)
		
		add_table <- info_qtls[which(info_qtls$Marker == add_bi_i), ]
		add_table$Marker <- substr(add_table$Marker, 1, 14)
		add_table$Marker <- paste(add_table$Marker, "2", sep = "")
		
		add_table$Qeff <- -add_table$Qeff
		add_table$lowerCI95 <- -add_table$upperCI95
		add_table$upperCI95 <- -add_table$lowerCI95
		
		info_qtls <- rbind(info_qtls, add_table)
		rm(add_table)
	}
}

#
# write out final set of haplotypes and the corresponding effect tables
#

info_qtls_sig_95_final <- rbind(info_qtls_sig_95, info_add_haps_sig95)
info_qtls_sig_95_final <- info_qtls_sig_95_final[order(info_qtls_sig_95_final$Marker, info_qtls_sig_95_final$Environment), ]
str(info_qtls_sig_95_final)

info_qtls_all_final <- rbind(info_qtls, info_add_haps)
info_qtls_all_final <- info_qtls_all_final[order(info_qtls_all_final$Marker, info_qtls_all_final$Environment), ]
str(info_qtls_all_final)

final_HapSet <- sort(final_HapSet)
final_HapSet

save(info_qtls_sig_95_final, info_qtls_all_final, final_HapSet, file = paste(outfolder, "/final_HapSet_", TRAIT, ".RData", sep = ""))
write.table(info_qtls_sig_95_final, file = paste(outfolder, "/final_HapSet_", TRAIT, ".csv", sep = ""), sep = ";", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)

###
###################################################################
