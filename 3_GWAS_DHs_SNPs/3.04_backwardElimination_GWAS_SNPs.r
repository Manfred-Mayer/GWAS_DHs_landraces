###################################################################
###################################################################
####
#### in the previous scripts, we defined trait-associated genomic regions for each environment 
#### and within each region one focus SNP (the one with the lowest p-value) was selected
#### these are "candidate SNPs"
####
#### now, we perform the backward elimination to obtain the final set of SNPs
#### also, based on the multi-locus multi-environment model, we extract the corresponding
#### environment-specific allelic effects
####
#### Manfred Mayer (Technical University of Munich, Plant Breeding)
#### manfred.mayer@tum.de
####
#### date: 29.06.2020
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
dir.create(paste("GWAS_SNPs/", TRAIT, "/QTLregs/finalQTL", sep = ""))
outfolder <- paste("GWAS_SNPs/", TRAIT, "/QTLregs/finalQTL", sep = "")
outfolder

# load data
# geno and map
load(paste("geno_InfoList_DH_SNPs.RData", sep = ""))
str(InfoList)
str(geno)

# environments
Envs <- read.table(paste0("GWAS_SNPs/", TRAIT, "/", list.files(path = paste0("GWAS_SNPs/", TRAIT), pattern = ".phenotypes_order.txt")), stringsAsFactors = FALSE)
Envs <- Envs[,1]
Envs

# load regions defined within environments
reg_all <- NULL
	for(env_i in Envs){
		if(file.exists(paste(infolder, "/reg_list_merge_", TRAIT , ".", env_i, "_FDR", FDR, ".txt", sep =""))){
			res <- read.table(paste(infolder, "/reg_list_merge_", TRAIT , ".", env_i, "_FDR", FDR, ".txt", sep =""), header = TRUE)
			reg_all <- rbind(reg_all, res)
			rm(res)
		}
	}
str(reg_all)

if(length(reg_all) > 0){
	
# generate input dataframe for model
	# order regions according to pvalue
	reg_all <- reg_all[order(reg_all$pvalue_log10, decreasing = TRUE), ]
	# remove redundant markers
	uni <- seq(1, nrow(reg_all), 1)
	names(uni) <- reg_all$mostSigMarker
	uni <- uni[unique(reg_all$mostSigMarker)]
	reg_all <- reg_all[uni , ]
	print(str(reg_all))
	
	# define candidate set of haplotypes
	QTL_cand <- reg_all$mostSigMarker
	candQTLset <- QTL_cand

	# generate input dataframe
	Y <- NULL
	ENV <- NULL
	GROUP <- NULL
	GENOTYPE <- NULL
	QTLs <- NULL
		# combine all environments (within-environment BLUEs), but exclude across-environment BLUEs
		for(env_i in Envs[which(Envs != "Across")]){
			Y_temp <- read.table(paste0("GWAS_SNPs/", TRAIT, "/", list.files(path = paste0("GWAS_SNPs/", TRAIT), pattern = ".pheno.txt")), stringsAsFactors = FALSE)
			colnames(Y_temp) <- Envs
			Y_temp <- Y_temp[ , env_i]

			GENOTYPE_temp <- read.table(paste0("GWAS_SNPs/", TRAIT, "/", list.files(path = paste0("GWAS_SNPs/", TRAIT), pattern = "geno_order.txt")), stringsAsFactors = FALSE)
			GENOTYPE_temp <- GENOTYPE_temp[,1]
			GENOTYPE_temp <- GENOTYPE_temp[which(is.na(Y_temp) == FALSE)]

			geno_gwas <- geno[GENOTYPE_temp, ]
			
			GROUP_temp <- read.table(paste0("GWAS_SNPs/", TRAIT, "/", list.files(path = paste0("GWAS_SNPs/", TRAIT), pattern = ".covariables.txt")), stringsAsFactors = FALSE)
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

##
## check if all chosen markers are polymorph within all environments
##
	mono_f <- function(x){
		sd(x, na.rm = TRUE) == 0
	}

	mono <- NULL		
	for(env_i in Envs[which(Envs != "Across")]){
		geno_i <- as.matrix(input_data[which(input_data$ENV == env_i) , -which(colnames(input_data) %in% c("Y", "ENV", "GROUP", "GENOTYPE"))])
		mono_i <- which(apply(geno_i, 2, mono_f))
		if(length(mono_i) > 0){
			mono_i <- names(mono_i)
			env <- rep(env_i, length(mono_i))
			res <- data.frame(ENV = env, Monomorph = mono_i)
			mono <- rbind(mono, res)
		}
	}
	print(mono)
	# exclude markers which are monomorphic within particular environments, because that would cause singularity in ASreml analysis
	# (that can happen in rare cases, as not all individuals were phenotyped in all environments (compare Hölker et al. 2019))
	print(dim(input_data))
	if(length(mono > 0)){
	ex <- unique(mono$Monomorph)
	input_data <- input_data[ , which((colnames(input_data)) %in% ex == FALSE)]
	QTLs <- input_data[ , -which(colnames(input_data) %in% c("Y", "ENV", "GROUP", "GENOTYPE"))]
	candQTLset <- colnames(QTLs)
	}
	print(dim(input_data))

##
## check how many lines per landrace are present within each environment
##

n_per_group <- matrix(NA, nrow = length(which(Envs != "Across")), ncol = 3)
rownames(n_per_group) <- Envs[which(Envs != "Across")]
colnames(n_per_group) <- c("KE", "LL", "PE")
for(env_i in Envs[which(Envs != "Across")]){
	geno_i <- as.matrix(input_data[which(input_data$ENV == env_i) , -which(colnames(input_data) %in% c("Y", "ENV", "GROUP", "GENOTYPE"))])
	res <- table(substr(rownames(geno_i), 1, 5))
	n_per_group[env_i, ] <- res
}
print(n_per_group)


###
### multi-locus, multi-environment model
###

	### 1) initial sorting of QTLs

	##
	## 1.1 first a dummy call to get starting values:
	##
	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(colnames(QTLs), collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08,
					   start.values=TRUE)
	iv     <- reml.obj$gammas.table
	
	###
	### 1.2 then run asreml for the initial model
	###
	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(colnames(QTLs), collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08)
	iv$Value  <- summary(reml.obj)$varcomp$component

	###
	### 1.3 sort markers based on pvalue
	###
	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(colnames(QTLs), collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08,
					   G.param=iv, 
					   R.param=iv)
	iv$Value  <- summary(reml.obj)$varcomp$component
	# get the pvalue of the wald test for the markers and sort them based on the pvalue
	wald.obj <- wald(reml.obj, denDF = "default")
	p_wald <- wald.obj$Wald$Pr
	names(p_wald) <- rownames(wald.obj$Wald)
	p_wald <- p_wald[which(substr(names(p_wald), 1, 4) == "ENV:")]
	p_wald <- sort(p_wald)
	input_data <- input_data[ , c("Y", "ENV", "GROUP", "GENOTYPE", gsub("ENV:", "", names(p_wald)))]

	
	
	### 2) backward selection
	
	###
	### 2.1 run asreml again for the initial model (with already ordered markers)
	###
	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(names(p_wald), collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08)
	# save variance components
	iv$Value  <- summary(reml.obj)$varcomp$component

	# generate initial vector of p-values (as template)
	p_values_new <- p_wald
	names(p_values_new) <- gsub("ENV:", "", names(p_values_new))

	###
	### 2.2 iteratively perform backward elimination (this can take a while ...)
	###     test all markers, always putting the tested marker as last entry in the model
	###     record p-values for each tested marker
	###     remove the one with the highest p-value
	###     until all p-values are < 0.01
	###
	
	removed <- NULL
	iter_n <- 1
	repeat{
	print(paste("###################", sep = ""))
	print(paste("###################", sep = ""))
	print(paste("###################", sep = ""))
	print(paste("Iteration ", iter_n, sep = ""))
	print(paste("", sep = ""))
	iter_n <- iter_n + 1
	
	# the significance of every QTL can only be interpreted if it is the last fixed effect in the model,
	# thus, circle through the QTLs and always record the p-value of the last one
	p_values_old <- p_values_new
	p_values_old <- sort(p_values_old)
	QTL_order_list <- list()
	for(qtl_i in names(p_values_old)){
		temp <- names(p_values_old)
		temp <- temp[-which(temp == qtl_i)]
		temp <- c(temp, qtl_i)
		QTL_order_list[[qtl_i]] <- paste0(temp, collapse = " + ")
	}
	p_values_new[] <- NA
		
	for(qtl_i in names(p_values_new)){
		print(paste("QTL ", which(names(p_values_new) == qtl_i), " of ", length(p_values_new), sep = ""))
		
			reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
										ENV:(", QTL_order_list[[qtl_i]],")", collapse = "")),
						   random = ~ GENOTYPE,
						   rcov   = ~ at(ENV):units,
						   data   = input_data,
						   maxit  = 500,
						   workspace   = 3e+08,
						   G.param=iv, 
						   R.param=iv)

			## get the pvalue of the wald test for the marker
			wald.obj <- wald(reml.obj, denDF = "default")
			p_wald <- wald.obj$Wald$Pr
			names(p_wald) <- rownames(wald.obj$Wald)
			p_wald <- p_wald[which(substr(names(p_wald), 1, 4) == "ENV:")]
			names(p_wald) <- gsub("ENV:", "", names(p_wald))
			p_values_new[qtl_i] <- p_wald[qtl_i]
		
			if(qtl_i == names(p_values_new)[length(p_values_new)]){
				# store the estimated variance components; 
				# to be used as starting values in the next asreml call
				iv$Value  <- summary(reml.obj)$varcomp$component
			}
	}
	print(summary(p_values_new))
	p_values_new <- sort(p_values_new)
	
		## remove the last snp if the threshold has not been reached yet
		if(p_values_new[length(p_values_new)] > p_thresh){
			removed <- c(removed, names(p_values_new)[(length(p_values_new))])
			input_data <- input_data[ , c("Y", "ENV", "GROUP", "GENOTYPE", names(p_values_new)[1:(length(p_values_new)-1)])]
			print(paste(names(p_values_new)[(length(p_values_new))], "removed, p-value:", p_values_new[length(p_values_new)]))
			p_values_new <- p_values_new[-length(p_values_new)]
			save(removed, file = paste(outfolder, "/removed_qtls.RData", sep =""))
		} else {
		## if all p-values are below threshold, stop
			break
		}
	}

	### 3) run final model and extract haplotype effects and variance components

	## 3.5 order final set of QTLs again according to their p-value and run the final model
		input_data <- input_data[ , c("Y", "ENV", "GROUP", "GENOTYPE", names(p_values_new)[1:(length(p_values_new))])]
			final_model <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
										ENV:(", paste0(names(p_values_new), collapse = " + "),")", collapse = "")),
						   random = ~ GENOTYPE,
						   rcov   = ~ at(ENV):units,
						   data   = input_data,
						   maxit  = 500,
						   workspace   = 3e+08,
						   G.param=iv, 
						   R.param=iv)
		p_wald <- p_values_new

		# extract effect estimates
		Qeff <- final_model$coefficients$fixed[grep(":AX.",names(final_model$coefficients$fixed))]
		seQeff <- sqrt(final_model$vcoeff$fixed)[1:length(Qeff)]
		
		# "Confidence intervals" (95%)
		lowerCI95 <- Qeff - qnorm(0.975, mean = 0, sd = 1) * seQeff
		upperCI95 <- Qeff + qnorm(0.975, mean = 0, sd = 1) * seQeff

		# summary table
		EnvMarker <- strsplit(names(Qeff), split = ":", fixed = TRUE)
		Marker <- unlist(lapply(EnvMarker, function(x){x[2]}))
		Marker <- gsub("X.", "X-", Marker)
		Environment <- unlist(lapply(EnvMarker, function(x){x[1]}))
		Environment <- gsub("ENV_", "", Environment)
		p_wald <- p_wald[Marker]
		info_qtls <- data.frame(Marker, Environment, p_wald, Qeff, seQeff, lowerCI95, upperCI95)

		# add position information
		info_qtls$Chr <- InfoList[info_qtls$Marker, 1]
		info_qtls$Start_bp <- InfoList[info_qtls$Marker, 2]
		info_qtls$End_bp <- InfoList[info_qtls$Marker, 3]
		info_qtls$Pos_bp <- apply(InfoList[info_qtls$Marker, 2:3], 1, mean)
		info_qtls <- info_qtls[order(info_qtls$Chr, info_qtls$Start_bp, info_qtls$End_bp, info_qtls$Environment), ]
		str(info_qtls)
		
		# write out final haplotype set
		finalQTLset <- unique(info_qtls$Marker)
		save(finalQTLset, file = paste(outfolder, "/finalQTLset_", TRAIT , "_FDR", FDR, "_p", p_thresh, ".RData", sep =""))

		# select significant effects (environment-specific)
		# do CIs include 0 or not
		sig_f <- function(x){
			(findInterval(0, x) == 1) == FALSE
		}		
		sig_m <- apply(info_qtls[, c("lowerCI95", "upperCI95")], 1, sig_f)
		info_qtls_sig_95 <- info_qtls[which(sig_m), ]
		str(info_qtls_sig_95)

	
	### compare model with and without QTLs	
	
	# calculate null-model without any haplotype as fixed effects
		# first get template
		null_model <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP", collapse = "")),
						   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08,
					   start.values=TRUE)
		iv <- null_model$gammas.table
		# update starting values
		null_model <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP", collapse = "")),
						   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08)
		iv$Value  <- summary(null_model)$varcomp$component
		# now, use starting values from preveous fit
		null_model <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP", collapse = "")),
						   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08)	
	
	# calculate the proportion of genotypic variance explained by the final set of haplotypes	
		varcomp_null <- summary(null_model)$varcomp$component
		varcomp_full <- summary(final_model)$varcomp$component
		diff_varcomp <- varcomp_null - varcomp_full
		
		print(paste("prop explained G: ", diff_varcomp[1] / varcomp_null[1], sep = ""))
		prop_explained <- diff_varcomp / varcomp_null
		varcomps <- cbind(varcomp_null, varcomp_full, diff_varcomp, prop_explained)
		rownames(varcomps) <- rownames(summary(null_model)$varcomp)

	# write output
		save(info_qtls_sig_95, info_qtls, null_model, final_model, finalQTLset, candQTLset, file = paste(outfolder, "/OutFinalQTLset_", TRAIT , "_FDR", FDR, "_p", p_thresh, ".RData", sep =""))
		write.table(info_qtls_sig_95, file = paste(outfolder, "/info_qtls_sig95_", TRAIT , "_FDR", FDR, "_p", p_thresh, ".csv", sep =""), na = "", sep = ";", col.names = TRUE, row.names = FALSE, quote = FALSE)
		write.table(table(info_qtls_sig_95$Environment), file = paste(outfolder, "/nEnv_sig95_", TRAIT , "_FDR", FDR, "_p", p_thresh, ".csv", sep =""), sep = ";", col.names = FALSE, row.names = FALSE, quote = FALSE)
		write.table(varcomps, file = paste(outfolder, "/varcompsFinalQTLset_", TRAIT , "_FDR", FDR, "_p", p_thresh, ".csv", sep =""), sep = ";", col.names = NA, row.names = TRUE, quote = FALSE)
		
}

###
###################################################################
