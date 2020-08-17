###################################################################
###################################################################
####
#### for focus haplotypes present in two landraces,
#### compare effect directions and significance of the haplotype between landraces
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

###################################################################
###
### 
###

dir.create("QTL_LRspecificEff")
outfolder <- paste("QTL_LRspecificEff", sep = "")
outfolder

# geno and map
load(paste("geno_InfoList_DH_m", nSNPs, "s", steps, ".RData", sep = ""))
str(InfoList)
str(geno)

traits <- c("EV_V4", "EV_V6", "PH_V4", "PH_V6",
			"PH_final",
			"FF", "MF",
			"LO", "TILL")

for(TRAIT in traits){

Envs <- read.table(paste0(TRAIT, "/", list.files(path = TRAIT, pattern = ".phenotypes_order.txt")), stringsAsFactors = FALSE)
Envs <- Envs[,1]
Envs

load(paste(TRAIT, "/QTLregs/finalQTL/OutFinalQTLset_", TRAIT, "_FDR0.15_p0.01.RData", sep = ""))	
QTL_cand <- finalQTLset

	# generate input dataframe
	Y <- NULL
	ENV <- NULL
	GROUP <- NULL
	GENOTYPE <- NULL
	QTLs <- NULL
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
	if(length(QTL_cand) == 1){
		colnames(input_data)[5] <- QTL_cand
	}
	print(str(input_data))

# sorting of markers
x <- final_model$fixed.formula
x <- strsplit(as.character(x)[3], split = "+", fixed = TRUE)
x <- gsub(" ", "", x[[1]])
x <- gsub("\n", "", x)
x <- gsub(")", "", x)
x <- x[-(1:2)]
x[1] <- substr(x[1], 6, nchar(x[1]))
x

	###
	### 1. first a dummy call to get starting values:
	###
	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:GROUP:(", paste0(x, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08,
					   start.values=TRUE)	
	iv     <- reml.obj$gammas.table

	###
	### 2. then run asreml model
	###
	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:GROUP:(", paste0(x, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08
					   )
	iv$Value  <- summary(reml.obj)$varcomp$component

	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:GROUP:(", paste0(x, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace = 3e+08,
					   G.param=iv, 
					   R.param=iv)

			## get the pvalue of the wald test for the snp
			wald.obj <- wald(reml.obj, denDF = "default")
			p_wald <- wald.obj$Wald$Pr
			names(p_wald) <- rownames(wald.obj$Wald)
			p_wald
			
	# extract effect estimates
	Qeff <- reml.obj$coefficients$fixed[grep(":wind",names(reml.obj$coefficients$fixed))]
	seQeff <- sqrt(reml.obj$vcoeff$fixed)[1:length(Qeff)]
		
	# "Confidence intervals" (95%)
	lowerCI95 <- Qeff - qnorm(0.975, mean = 0, sd = 1) * seQeff
	upperCI95 <- Qeff + qnorm(0.975, mean = 0, sd = 1) * seQeff

	# summary table
		EnvGroupMarker <- strsplit(names(Qeff), split = ":", fixed = TRUE)
		Marker <- unlist(lapply(EnvGroupMarker, function(x){x[3]}))
		Landrace <- unlist(lapply(EnvGroupMarker, function(x){x[2]}))
		Landrace <- gsub("GROUP_", "", Landrace)
		Environment <- unlist(lapply(EnvGroupMarker, function(x){x[1]}))
		Environment <- gsub("ENV_", "", Environment)
		names(p_wald) <- gsub("ENV:GROUP:", "", names(p_wald))
		p_wald <- p_wald[Marker]
		info_qtls_LReff <- data.frame(Marker, Environment, Landrace, p_wald, Qeff, seQeff, lowerCI95, upperCI95)

	# for each haplotype, write out for each landrace
	# and for each environment where it was significant in the combined analysis
	# - is it significant for the particular landrace
	# - the effect within the landrace
	# - the number of individuals from that landrace in the particular environment
	info_res <- NULL
	
	for(hap_i in unique(info_qtls_LReff$Marker)){
		info_qtls_LReff_i <- info_qtls_LReff[which(info_qtls_LReff$Marker == hap_i), ]
		info_qtls_sig95_i <- info_qtls_sig_95[which(info_qtls_sig_95$Marker == hap_i), ]
	
		info_qtls_LReff_i <- info_qtls_LReff_i[which(info_qtls_LReff_i$Environment %in% info_qtls_sig95_i$Environment), ]
		
		n_env_sig <- nrow(info_qtls_sig95_i)
		
		for(env_i in unique(info_qtls_LReff_i$Environment)){
			info_qtls_LReff_i_env <- info_qtls_LReff_i[which(info_qtls_LReff_i$Environment == env_i), ]
			Eff_KE <- info_qtls_LReff_i_env$Qeff[which(info_qtls_LReff_i_env$Landrace == "KE")]
			Eff_LL <- info_qtls_LReff_i_env$Qeff[which(info_qtls_LReff_i_env$Landrace == "LL")]
			Eff_PE <- info_qtls_LReff_i_env$Qeff[which(info_qtls_LReff_i_env$Landrace == "PE")]
			
			sig_Eff_KE <- ifelse(	(info_qtls_LReff_i_env$lowerCI95[which(info_qtls_LReff_i_env$Landrace == "KE")] <= 0) &
									(info_qtls_LReff_i_env$upperCI95[which(info_qtls_LReff_i_env$Landrace == "KE")] >= 0),
									FALSE,
									TRUE)
			sig_Eff_LL <- ifelse(	(info_qtls_LReff_i_env$lowerCI95[which(info_qtls_LReff_i_env$Landrace == "LL")] <= 0) &
									(info_qtls_LReff_i_env$upperCI95[which(info_qtls_LReff_i_env$Landrace == "LL")] >= 0),
									FALSE,
									TRUE)
			sig_Eff_PE <- ifelse(	(info_qtls_LReff_i_env$lowerCI95[which(info_qtls_LReff_i_env$Landrace == "PE")] <= 0) &
									(info_qtls_LReff_i_env$upperCI95[which(info_qtls_LReff_i_env$Landrace == "PE")] >= 0),
									FALSE,
									TRUE)
			
			pheno_temp <- read.table(paste0(TRAIT, "/", list.files(path = TRAIT, pattern = ".pheno.txt")), stringsAsFactors = FALSE)
			colnames(pheno_temp) <- Envs
			geno_names <- read.table(paste0(TRAIT, "/", list.files(path = TRAIT, pattern = "geno_order.txt")), stringsAsFactors = FALSE)
			pheno_env_i <- pheno_temp[, env_i]
			names(pheno_env_i) <- geno_names[,1]

			pheno_env_i_KE <- pheno_env_i[which(substr(names(pheno_env_i), 1, 5) == "DH_KE")]
			pheno_env_i_KE <- pheno_env_i_KE[which(is.na(pheno_env_i_KE) == FALSE)]
			pheno_env_i_KE <- pheno_env_i_KE[which(names(pheno_env_i_KE) %in% rownames(geno)[which(geno[, hap_i] == 2)])]
			n_KE <- length(pheno_env_i_KE)
			
			pheno_env_i_LL <- pheno_env_i[which(substr(names(pheno_env_i), 1, 5) == "DH_LL")]
			pheno_env_i_LL <- pheno_env_i_LL[which(is.na(pheno_env_i_LL) == FALSE)]
			pheno_env_i_LL <- pheno_env_i_LL[which(names(pheno_env_i_LL) %in% rownames(geno)[which(geno[, hap_i] == 2)])]
			n_LL <- length(pheno_env_i_LL)

			pheno_env_i_PE <- pheno_env_i[which(substr(names(pheno_env_i), 1, 5) == "DH_PE")]
			pheno_env_i_PE <- pheno_env_i_PE[which(is.na(pheno_env_i_PE) == FALSE)]
			pheno_env_i_PE <- pheno_env_i_PE[which(names(pheno_env_i_PE) %in% rownames(geno)[which(geno[, hap_i] == 2)])]
			n_PE <- length(pheno_env_i_PE)
			
			if(n_KE == 0){
				Eff_KE <- NA
				sig_Eff_KE <- NA
			}
			
			if(n_LL == 0){
				Eff_LL <- NA
				sig_Eff_LL <- NA
			}
			
			if(n_PE == 0){
				Eff_PE <- NA
				sig_Eff_PE <- NA
			}
			
			info_i <- data.frame(Haplotype = hap_i,
								 Environment = env_i,
								 Eff_KE = Eff_KE,
								 Eff_LL = Eff_LL,
								 Eff_PE = Eff_PE,
								 sig_Eff_KE = sig_Eff_KE,
								 sig_Eff_LL = sig_Eff_LL,
								 sig_Eff_PE = sig_Eff_PE,
								 n_KE = n_KE,
								 n_LL = n_LL,
								 n_PE = n_PE
								)
			
			info_res <- rbind(info_res, info_i)
		}
		
	}

# for some haplotypes, effects could not be estimated due to collinearity (within landraces)
# for those ones, test them separately in univariate models
hap_uni <- NULL
print(info_res[which(info_res$Eff_KE == 0), ])
hap_uni <- c(hap_uni, unique(info_res[which(info_res$Eff_KE == 0), "Haplotype"]))
print(info_res[which(info_res$Eff_LL == 0), ])
hap_uni <- c(hap_uni, unique(info_res[which(info_res$Eff_LL == 0), "Haplotype"]))
print(info_res[which(info_res$Eff_PE == 0), ])
hap_uni <- c(hap_uni, unique(info_res[which(info_res$Eff_PE == 0), "Haplotype"]))
hap_uni <- sort(unique(hap_uni))
print(hap_uni)

if(length(hap_uni) > 0){

	for (hap_i in hap_uni){
		
	###
	### 1. first a dummy call to get starting values:
	###
	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:GROUP:(", paste0(hap_i, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08,
					   start.values=TRUE)	
	iv     <- reml.obj$gammas.table

	###
	### 2. then run asreml model
	###
	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:GROUP:(", paste0(hap_i, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08
					   )
	iv$Value  <- summary(reml.obj)$varcomp$component

	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:GROUP:(", paste0(hap_i, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace = 3e+08,
					   G.param=iv, 
					   R.param=iv)

			## get the pvalue of the wald test for the snp
			wald.obj <- wald(reml.obj, denDF = "default")
			p_wald <- wald.obj$Wald$Pr
			names(p_wald) <- rownames(wald.obj$Wald)
			p_wald
			
	# extract effect estimates
	Qeff <- reml.obj$coefficients$fixed[grep(":wind",names(reml.obj$coefficients$fixed))]
	seQeff <- sqrt(reml.obj$vcoeff$fixed)[1:length(Qeff)]
		
	# "Confidence intervals" (95%)
	lowerCI95 <- Qeff - qnorm(0.975, mean = 0, sd = 1) * seQeff
	upperCI95 <- Qeff + qnorm(0.975, mean = 0, sd = 1) * seQeff

	# summary table
		EnvGroupMarker <- strsplit(names(Qeff), split = ":", fixed = TRUE)
		Marker <- unlist(lapply(EnvGroupMarker, function(x){x[3]}))
		Landrace <- unlist(lapply(EnvGroupMarker, function(x){x[2]}))
		Landrace <- gsub("GROUP_", "", Landrace)
		Environment <- unlist(lapply(EnvGroupMarker, function(x){x[1]}))
		Environment <- gsub("ENV_", "", Environment)
		names(p_wald) <- gsub("ENV:GROUP:", "", names(p_wald))
		p_wald <- p_wald[Marker]
		info_qtls_LReff_2 <- data.frame(Marker, Environment, Landrace, p_wald, Qeff, seQeff, lowerCI95, upperCI95)

		info_qtls_LReff_i <- info_qtls_LReff_2[which(info_qtls_LReff_2$Marker == hap_i), ]
		info_qtls_sig95_i <- info_qtls_sig_95[which(info_qtls_sig_95$Marker == hap_i), ]
	
		info_qtls_LReff_i <- info_qtls_LReff_i[which(info_qtls_LReff_i$Environment %in% info_qtls_sig95_i$Environment), ]
		
		n_env_sig <- nrow(info_qtls_sig95_i)
		
		for(env_i in unique(info_qtls_LReff_i$Environment)){
			info_qtls_LReff_i_env <- info_qtls_LReff_i[which(info_qtls_LReff_i$Environment == env_i), ]
			Eff_KE <- info_qtls_LReff_i_env$Qeff[which(info_qtls_LReff_i_env$Landrace == "KE")]
			Eff_LL <- info_qtls_LReff_i_env$Qeff[which(info_qtls_LReff_i_env$Landrace == "LL")]
			Eff_PE <- info_qtls_LReff_i_env$Qeff[which(info_qtls_LReff_i_env$Landrace == "PE")]
			
			sig_Eff_KE <- ifelse(	(info_qtls_LReff_i_env$lowerCI95[which(info_qtls_LReff_i_env$Landrace == "KE")] <= 0) &
									(info_qtls_LReff_i_env$upperCI95[which(info_qtls_LReff_i_env$Landrace == "KE")] >= 0),
									FALSE,
									TRUE)
			sig_Eff_LL <- ifelse(	(info_qtls_LReff_i_env$lowerCI95[which(info_qtls_LReff_i_env$Landrace == "LL")] <= 0) &
									(info_qtls_LReff_i_env$upperCI95[which(info_qtls_LReff_i_env$Landrace == "LL")] >= 0),
									FALSE,
									TRUE)
			sig_Eff_PE <- ifelse(	(info_qtls_LReff_i_env$lowerCI95[which(info_qtls_LReff_i_env$Landrace == "PE")] <= 0) &
									(info_qtls_LReff_i_env$upperCI95[which(info_qtls_LReff_i_env$Landrace == "PE")] >= 0),
									FALSE,
									TRUE)
			
			pheno_temp <- read.table(paste0(TRAIT, "/", list.files(path = TRAIT, pattern = ".pheno.txt")), stringsAsFactors = FALSE)
			colnames(pheno_temp) <- Envs
			geno_names <- read.table(paste0(TRAIT, "/", list.files(path = TRAIT, pattern = "geno_order.txt")), stringsAsFactors = FALSE)
			pheno_env_i <- pheno_temp[, env_i]
			names(pheno_env_i) <- geno_names[,1]

			pheno_env_i_KE <- pheno_env_i[which(substr(names(pheno_env_i), 1, 5) == "DH_KE")]
			pheno_env_i_KE <- pheno_env_i_KE[which(is.na(pheno_env_i_KE) == FALSE)]
			pheno_env_i_KE <- pheno_env_i_KE[which(names(pheno_env_i_KE) %in% rownames(geno)[which(geno[, hap_i] == 2)])]
			n_KE <- length(pheno_env_i_KE)
			
			pheno_env_i_LL <- pheno_env_i[which(substr(names(pheno_env_i), 1, 5) == "DH_LL")]
			pheno_env_i_LL <- pheno_env_i_LL[which(is.na(pheno_env_i_LL) == FALSE)]
			pheno_env_i_LL <- pheno_env_i_LL[which(names(pheno_env_i_LL) %in% rownames(geno)[which(geno[, hap_i] == 2)])]
			n_LL <- length(pheno_env_i_LL)

			pheno_env_i_PE <- pheno_env_i[which(substr(names(pheno_env_i), 1, 5) == "DH_PE")]
			pheno_env_i_PE <- pheno_env_i_PE[which(is.na(pheno_env_i_PE) == FALSE)]
			pheno_env_i_PE <- pheno_env_i_PE[which(names(pheno_env_i_PE) %in% rownames(geno)[which(geno[, hap_i] == 2)])]
			n_PE <- length(pheno_env_i_PE)
			
			if(n_KE == 0){
				Eff_KE <- NA
				sig_Eff_KE <- NA
			}
			
			if(n_LL == 0){
				Eff_LL <- NA
				sig_Eff_LL <- NA
			}
			
			if(n_PE == 0){
				Eff_PE <- NA
				sig_Eff_PE <- NA
			}
			
			info_i <- data.frame(Haplotype = paste0(hap_i, "_uni"),
								 Environment = env_i,
								 Eff_KE = Eff_KE,
								 Eff_LL = Eff_LL,
								 Eff_PE = Eff_PE,
								 sig_Eff_KE = sig_Eff_KE,
								 sig_Eff_LL = sig_Eff_LL,
								 sig_Eff_PE = sig_Eff_PE,
								 n_KE = n_KE,
								 n_LL = n_LL,
								 n_PE = n_PE
								)
			
			info_res <- rbind(info_res, info_i)
		}
		
	}

}
	
save(info_res, file = paste(outfolder, "/QTL_LRspecificEff_", TRAIT, "_onlySigEnv.RData", sep = ""))
write.table(info_res, file = paste(outfolder, "/QTL_LRspecificEff_", TRAIT, "_onlySigEnv.csv", sep = ""), sep = ";", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)


#
# as LL has a low number of individuals,
# only compare effects between KE and PE
#

  hapEnv_1LR <- NULL
  hapEnv_1LR_KE <- NULL
  hapEnv_1LR_PE <- NULL
  hapEnv_2LR <- NULL

  significant <- NULL
  not_significant <- NULL
  
  equal_sign <- NULL
  opposite_sign <- NULL
  
  equal_sign_significant <- NULL
  opposite_sign_significant <- NULL
  
  info_res <- info_res[ , -grep("LL", colnames(info_res))]
  uni_i <- info_res[grep("_uni", info_res$Haplotype), ]
  multi_i <- info_res[which(info_res$Haplotype %in% uni_i$Haplotype == FALSE), ]
  
  # check each haplotype
  for(hap_i in unique(multi_i$Haplotype)){
    res_i <- multi_i[which(multi_i$Haplotype == hap_i), ]
    
    if(any(res_i$Eff_KE == 0, na.rm = TRUE)){
      env_i <- res_i$Environment[which(res_i$Eff_KE == 0)]
      uni_replace <- uni_i[which((uni_i$Haplotype == paste0(hap_i, "_uni")) & (uni_i$Environment %in% env_i)), "Eff_KE"]
      res_i$Eff_KE[which(res_i$Environment %in% env_i)] <- uni_replace
    }
    
    if(any(res_i$Eff_PE == 0, na.rm = TRUE)){
      env_i <- res_i$Environment[which(res_i$Eff_PE == 0)]
      uni_replace <- uni_i[which((uni_i$Haplotype == paste0(hap_i, "_uni")) & (uni_i$Environment %in% env_i)), "Eff_PE"]
      res_i$Eff_PE[which(res_i$Environment %in% env_i)] <- uni_replace
    }
  
   # take haplotype x environment combinations as units
    for(env_i in res_i$Environment){
      res_i_env_i <- res_i[which(res_i$Environment == env_i), ]
      
   # in how many landraces does this haplotype occur
      n_LR <- length(which(res_i_env_i[ , c("n_KE", "n_PE")] != 0))
      if(n_LR == 1){
        hapEnv_1LR <- c(hapEnv_1LR, paste(hap_i, env_i, sep = "."))
        if(res_i_env_i[ , c("n_KE")] == 0){
          hapEnv_1LR_PE <- c(hapEnv_1LR_PE, paste(hap_i, env_i, sep = "."))
        }
        if(res_i_env_i[ , c("n_PE")] == 0){
          hapEnv_1LR_KE <- c(hapEnv_1LR_KE, paste(hap_i, env_i, sep = "."))
        }
      } 
      if(n_LR == 2){
        hapEnv_2LR <- c(hapEnv_2LR, paste(hap_i, env_i, sep = "."))
      } 
      
    if(n_LR == 2){
    
    # compare effect sign between landraces
      eff_i <- res_i_env_i[ , c("Eff_KE", "Eff_PE")]
      eff_i <- eff_i[which(is.na(eff_i) == FALSE)]
      sign_i <- as.numeric(sign(eff_i))
        if(length(unique(sign_i)) == 1){
          equal_sign <- c(equal_sign, paste(hap_i, env_i, sep = "."))
          
          # check in addition if the effect is significant in the two landraces
          sig_i <- res_i_env_i[ , c("sig_Eff_KE", "sig_Eff_PE")]
          sig_i <- sig_i[which(is.na(sig_i) == FALSE)]
          sig_i <- as.logical(sig_i)
          n_sig_i <- sum(sig_i)
          if(n_sig_i > 1){
            equal_sign_significant <- c(equal_sign_significant, paste(hap_i, env_i, sep = "."))
            significant <- c(significant, paste(hap_i, env_i, sep = "."))
          } else {
            not_significant <- c(not_significant, paste(hap_i, env_i, sep = "."))
          }
        
        } else {
          
          opposite_sign <- c(opposite_sign, paste(hap_i, env_i, sep = "."))
          
          # check in addition if the effect is significant the two landraces
          sig_i <- res_i_env_i[ , c("sig_Eff_KE", "sig_Eff_PE")]
          sig_i <- sig_i[which(is.na(sig_i) == FALSE)]
          sig_i <- as.logical(sig_i)
          n_sig_i <- sum(sig_i)
          if(n_sig_i > 1){
            opposite_sign_significant <- c(opposite_sign_significant, paste(hap_i, env_i, sep = "."))
            significant <- c(significant, paste(hap_i, env_i, sep = "."))
          } else {
            not_significant <- c(not_significant, paste(hap_i, env_i, sep = "."))
          }
          
        }
      }

  }
  }

  hapEnv_1LR
  hapEnv_1LR_KE
  hapEnv_1LR_PE
  hapEnv_2LR

  equal_sign
  opposite_sign
  
  significant
  not_significant
  
  equal_sign_significant
  opposite_sign_significant

  # plot results
  png(paste(outfolder, "/summary_PEKE_effects_", TRAIT, ".png", sep = ""), height = 1500, width = 1400, res = 200)
  par(mai = c(0,0,0.8,0))
  plot(-10,-10, xlim = c(2, 48), ylim = c(3,57), main = TRAIT, xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", cex.main = 1.6, bty="n")
  
  points(x = c(45, 52) - 10, y = c(20, 5), type = "l", lwd = 3, col = "gray75")
  points(x = c(45, 40) - 10, y = c(20, 5), type = "l", lwd = 3, col = "gray75")
  points(x = c(25, 30) - 10, y = c(20, 5), type = "l", lwd = 3, col = "gray75")
  points(x = c(25, 18) - 10, y = c(20, 5), type = "l", lwd = 3, col = "gray75")
  points(x = c(35, 45) - 10, y = c(35, 20), type = "l", lwd = 3, col = "gray75")
  points(x = c(35, 25) - 10, y = c(35, 20), type = "l", lwd = 3, col = "gray75")
  points(x = c(25, 35) - 10, y = c(50, 35), type = "l", lwd = 3, col = "gray75")
  points(x = c(25, 15) - 10, y = c(50, 35), type = "l", lwd = 3, col = "gray75")
  
  points(x = 25 - 10, y = 50, cex = 0.5 + 25, col = "gray70", pch = 16)
  text(x = 25 - 10, y = 50, cex = 2, labels = length(hapEnv_1LR) + length(hapEnv_2LR))
  text(x = 25 - 10, y = 53, cex = 1, labels = "Hap*Env units")
  
  points(x = 15 - 10, y = 35, cex = 0.5 + 25*(length(hapEnv_1LR) / (length(hapEnv_1LR) + length(hapEnv_2LR))), col = "gray70", pch = 16)
  text(x = 15 - 10, y = 35, cex = 2, labels = length(hapEnv_1LR))
  text(x = 15 - 10, y = 38, cex = 1, labels = "In one LR")
  
  points(x = 35 - 10, y = 35, cex = 0.5 + 25*(length(hapEnv_2LR) / (length(hapEnv_1LR) + length(hapEnv_2LR))), col = "gray70", pch = 16)
  text(x = 35 - 10, y = 35, cex = 2, labels = length(hapEnv_2LR))
  text(x = 35 - 10, y = 38, cex = 1, labels = "In both LR")
  
  points(x = 25 - 10, y = 20, cex = 0.5 + 25*((length(equal_sign_significant) + length(opposite_sign_significant)) / (length(hapEnv_1LR) + length(hapEnv_2LR))), col = "gray70", pch = 16)
  text(x = 25 - 10, y = 20, cex = 2, labels = (length(equal_sign_significant) + length(opposite_sign_significant)))
  text(x = 25 - 10, y = 24.5, cex = 1, labels = "Significant")
  text(x = 25 - 10, y = 23, cex = 1, labels = "for both LRs")
  
  points(x = 45 - 10, y = 20, cex = 0.5 + 25*((length(hapEnv_2LR) - (length(equal_sign_significant) + length(opposite_sign_significant))) / (length(hapEnv_1LR) + length(hapEnv_2LR))), col = "gray70", pch = 16)
  text(x = 45 - 10, y = 20, cex = 2, labels = (length(hapEnv_2LR) - (length(equal_sign_significant) + length(opposite_sign_significant))))
  text(x = 45 - 10, y = 24.5, cex = 1, labels = "Not significant")
  text(x = 45 - 10, y = 23, cex = 1, labels = "for both LRs")
  
  points(x = 18 - 10, y = 5, cex = 0.5 + 25*(length(equal_sign_significant) / (length(hapEnv_1LR) + length(hapEnv_2LR))), col = "gray70", pch = 16)
  text(x = 18 - 10, y = 5, cex = 2, labels = length(equal_sign_significant))
  text(x = 18 - 10, y = 8, cex = 1, labels = "Equal sign")
  
  points(x = 30 - 10, y = 5, cex = 0.5 + 25*(length(opposite_sign_significant) / (length(hapEnv_1LR) + length(hapEnv_2LR))), col = "gray70", pch = 16)
  text(x = 30 - 10, y = 5, cex = 2, labels = length(opposite_sign_significant))
  text(x = 30 - 10, y = 8, cex = 1, labels = "Opposite sign")

  points(x = 40 - 10, y = 5, cex = 0.5 + 25*((length(hapEnv_2LR) - (length(equal_sign_significant) + length(opposite_sign_significant)) -
                                                (length(opposite_sign) - length(opposite_sign_significant))) / (length(hapEnv_1LR) + length(hapEnv_2LR))), col = "gray70", pch = 16)
  text(x = 40 - 10, y = 5, cex = 2, labels = (length(hapEnv_2LR) - (length(equal_sign_significant) + length(opposite_sign_significant)) -
                                                (length(opposite_sign) - length(opposite_sign_significant))))
  text(x = 40 - 10, y = 8, cex = 1, labels = "Equal sign")
  
  points(x = 52 - 10, y = 5, cex = 0.5 + 25*((length(hapEnv_2LR) - (length(equal_sign_significant) + length(opposite_sign_significant)) -
                                                (length(equal_sign) - length(equal_sign_significant))) / (length(hapEnv_1LR) + length(hapEnv_2LR))), col = "gray70", pch = 16)
  text(x = 52 - 10, y = 5, cex = 2, labels = (length(hapEnv_2LR) - (length(equal_sign_significant) + length(opposite_sign_significant)) -
                                                (length(equal_sign) - length(equal_sign_significant))))
  text(x = 52 - 10, y = 8, cex = 1, labels = "Opposite sign")
  
  dev.off()

}

# check how many haplotypes are involved in each number (HxE combination) plotted above
info <- NULL

for(TRAIT in traits){
  print(TRAIT)
  
  load(paste(outfolder, "/QTL_LRspecificEff_", TRAIT, "_onlySigEnv.RData", sep = ""))
  info_res <- info_res[ , -grep("LL", colnames(info_res))]
  info_res <- info_res[-grep("_uni", info_res$Haplotype), ]
  info_res <- info_res[which((info_res$n_KE > 0) | (info_res$n_PE > 0)), ]
  
  n_hap_all <- length(unique(info_res$Haplotype))
  
  info_res_both <- info_res[which((info_res$n_KE > 0) & (info_res$n_PE > 0)), ]
  n_hap_both <- length(unique(info_res_both$Haplotype))
  
  info_res_both_sig <- info_res_both[which(info_res_both$sig_Eff_KE & info_res_both$sig_Eff_PE), ]
  n_hap_both_sig <- length(unique(info_res_both_sig$Haplotype))
  
  info_res_both_notsig <- info_res_both[which((info_res_both$sig_Eff_KE & info_res_both$sig_Eff_PE) == FALSE), ]
  n_hap_both_notsig <- length(unique(info_res_both_notsig$Haplotype))
  
  info_res_both_sameDir <- info_res[which(sign(info_res$Eff_KE) == sign(info_res$Eff_PE)), ]
  n_hap_both_sameDir <- length(unique(info_res_both_sameDir$Haplotype))
  info_res_both_sameDir_sig <- info_res_both_sameDir[which(info_res_both_sameDir$sig_Eff_KE & info_res_both_sameDir$sig_Eff_PE), ]
  n_hap_both_sameDir_sig <- length(unique(info_res_both_sameDir_sig$Haplotype))
  
  # check for potentially additional haplotypes for which KE and PE are significant but in different environments
  info_res_both_sameDir_notsig <- info_res_both_sameDir[which(info_res_both_sameDir$sig_Eff_KE != info_res_both_sameDir$sig_Eff_PE), ]
  haps_temp <- NULL
  for(hap_i in unique(info_res_both_sameDir_notsig$Haplotype)){
    info_temp <- info_res_both_sameDir_notsig[which(info_res_both_sameDir_notsig$Haplotype == hap_i), ]
    any_sig_KE <- any(info_temp$sig_Eff_KE)
    any_sig_PE <- any(info_temp$sig_Eff_PE)
    if(any_sig_KE & any_sig_PE){
      haps_temp <- c(haps_temp, hap_i)
    }
  }
  add_haps <- length(haps_temp[which(haps_temp %in% info_res_both_sameDir_sig$Haplotype == FALSE)])
  
  info_res_both_oppoDir <- info_res[which(sign(info_res$Eff_KE) != sign(info_res$Eff_PE)), ]
  info_res_both_oppoDir <- info_res_both_oppoDir[which((info_res_both_oppoDir$Eff_KE != 0) & (info_res_both_oppoDir$Eff_PE != 0)), ]
  n_hap_both_oppoDir <- length(unique(info_res_both_oppoDir$Haplotype))
  info_res_both_oppoDir_sig <- info_res_both_sameDir[which(info_res_both_oppoDir$sig_Eff_KE & info_res_both_oppoDir$sig_Eff_PE), ]
  n_hap_both_oppoDir_sig <- length(unique(info_res_both_oppoDir_sig$Haplotype))
  
  info_i <- c(n_hap_all, n_hap_both, n_hap_both_sig, n_hap_both_notsig, n_hap_both_sameDir, n_hap_both_oppoDir, n_hap_both_sameDir_sig, add_haps, n_hap_both_oppoDir_sig)
  info <- rbind(info, info_i)
}
rownames(info) <- traits
colnames(info) <- c("n_hap_all", "n_hap_both", "n_hap_both_sig", "n_hap_both_notsig", "n_hap_both_sameDir", "n_hap_both_oppoDir", "n_hap_both_sameDir_sig", "add_haps_sig", "n_hap_both_oppoDir_sig")
info

###
###################################################################
