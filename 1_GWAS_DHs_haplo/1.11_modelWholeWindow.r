###################################################################
###################################################################
####
#### the selected focus haplotypes always represent the haplotype with the
#### smallest p-value when tested against the remaining haplotypes within a given window.
#### Now, we look at all haplotypes within a given window simultaneously and:
#### - calculate the proportion of genetic variance explained by the window (including all alternative haplotypes)
#### - the effect of each alternative haplotype within the respective window relative to the focus haplotype
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
library(asreml)

# graphical parameters
cex <- 0.6
cex.axis <- 1.6
cex.lab <- 1.6

# generate letter vector for naming haplotype variants
letters_MM <- c(letters[1:26], paste(rep(letters[1:26], rep(26,26)), rep(letters[1:26], 26), sep = ""))
letters_MM <- letters_MM[-1]
letters_MM

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

# geno and map
load(paste("geno_InfoList_DH_m", nSNPs, "s", steps, ".RData", sep = ""))
str(InfoList)
str(geno)

# environments for that trait
Envs <- read.table(paste0(TRAIT, "/", list.files(path = TRAIT, pattern = ".phenotypes_order.txt")), stringsAsFactors = FALSE)
Envs <- Envs[,1]
# in model we only need the within environment blues, so exclude the across environment blues
Envs <- Envs[-which(Envs == "Across")]
Envs

# load focus haplotypes
load(paste(TRAIT, "/QTLregs/finalQTL/OutFinalQTLset_", TRAIT, "_FDR0.15_p0.01.RData", sep = ""))	
QTL_cand <- finalQTLset

# sorting of markers in the model
x <- final_model$fixed.formula
x <- strsplit(as.character(x)[3], split = "+", fixed = TRUE)
x <- gsub(" ", "", x[[1]])
x <- gsub("\n", "", x)
x <- gsub(")", "", x)
x <- x[-(1:2)]
x[1] <- substr(x[1], 6, nchar(x[1]))
x

# for biallelic case add alternative haplotype
for(x_i in 1:length(x)){
wind <- substr(x[x_i], 1, 13)
hap_focus <- x[x_i]
if(length(which(substr(rownames(InfoList), 1, 13) == wind)) == 1){
temp_xix <- geno[ , which(substr(colnames(geno), 1, 13) == substr(hap_focus, 1, 13))]
temp_xix[which(temp_xix == 0)] <- 99
temp_xix[which(temp_xix == 2)] <- 0
temp_xix[which(temp_xix == 99)] <- 2
geno <- cbind(geno, temp_xix)
colnames(geno)[ncol(geno)] <- paste0(wind, "_99")
InfoList_xi <- InfoList[which(substr(rownames(InfoList), 1, 13) == substr(hap_focus, 1, 13)), ]
InfoList <- rbind(InfoList, InfoList_xi)
rownames(InfoList)[nrow(InfoList)] <- paste0(wind, "_99")
print(tail(InfoList, n = 10))
}
}
rm(wind)
rm(x_i)
rm(hap_focus)

###################################################################################################

dummy <- FALSE
if(length(x) == 1){
x <- c(x,x)
dummy <- TRUE
}

# make output folders
dir.create(paste("wholeWindModel", sep = ""))
dir.create(paste("wholeWindModel/", TRAIT, sep = ""))

# generate haplotype vector for the respective windows (coded as factor)
hap_vec_matrix <- matrix(NA, nrow = nrow(geno), ncol = length(x))
rownames(hap_vec_matrix) <- rownames(geno)
colnames(hap_vec_matrix) <- substr(x, 1, 13)

for(x_i in 1:length(x)){
wind <- colnames(hap_vec_matrix)[x_i]
hap_focus <- x[x_i]
geno_wind <- geno[,which(substr(colnames(geno), 1, 13) == wind)]

recode_f <- function(x) {
 y <- which(x == 2)
 if(length(y) == 0){
 y <- NA
 } else {
 y <- letters_MM[y]
 }
 return(y)
}
hap_vec <- apply(geno_wind, 1, recode_f)
# mark focus haplotype (always as "a")
hap_vec[which(geno_wind[, hap_focus] == 2)] <- "a"
# rename haplotypes, so that "b" is the most frequent, "c" the next most frequent, ...
hap_vec_new <- hap_vec
freqs <- sort(table(hap_vec), decreasing = TRUE)
freqs <- freqs[-which(names(freqs) == "a")]
for(i in 1:length(freqs)){
	hap_vec_new[which(hap_vec == names(freqs)[i])] <- letters_MM[i]
}
sort(table(hap_vec), decreasing = TRUE)
sort(table(hap_vec_new), decreasing = TRUE)
hap_vec <- hap_vec_new
rm(hap_vec_new)

hap_vec_matrix[ , wind] <- hap_vec
rm(hap_vec)
}

	# generate input dataframe (having the focus haplotypes number coded (as before) and the respective windows factor coded
	Y <- NULL
	ENV <- NULL
	GROUP <- NULL
	GENOTYPE <- NULL
	QTL <- NULL
	QTL2 <- NULL
		for(env_i in Envs){
			print(env_i)
			Y_temp <- read.table(paste(TRAIT, "/", list.files(path = TRAIT, pattern = ".pheno.txt"), sep = ""), stringsAsFactors = FALSE)
			colnames(Y_temp) <- Envs
			Y_temp <- Y_temp[ , env_i]

			GENOTYPE_temp <- read.table(paste(TRAIT, "/", list.files(path = TRAIT, pattern = "geno_order.txt"), sep = ""), stringsAsFactors = FALSE)
			GENOTYPE_temp <- GENOTYPE_temp[,1]
			names(Y_temp) <- GENOTYPE_temp

			GROUP_temp <- read.table(paste(TRAIT, "/", list.files(path = TRAIT, pattern = ".covariables.txt"), sep = ""), stringsAsFactors = FALSE)
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
			names(GROUP_temp) <- GENOTYPE_temp				
						
			GENOTYPE_temp <- GENOTYPE_temp[which(is.na(Y_temp) == FALSE)]
			
			QTLs_temp <- hap_vec_matrix[GENOTYPE_temp, ]
			QTLs2_temp <- geno[GENOTYPE_temp, x]
			GROUP_temp <- GROUP_temp[GENOTYPE_temp]
			Y_temp <- Y_temp[GENOTYPE_temp]					
			ENV_temp <- rep(env_i, length(Y_temp))
					
			Y <- c(Y, Y_temp)
			ENV <- c(ENV, ENV_temp)
			GENOTYPE <- c(GENOTYPE, GENOTYPE_temp)
			GROUP <- c(GROUP, GROUP_temp)
			QTL <- rbind(QTL, QTLs_temp)
			QTL2 <- rbind(QTL2, QTLs2_temp)
		}
	input_data <- data.frame(Y, ENV, GROUP, GENOTYPE, QTL, QTL2)
	input_data$ENV <- as.factor(input_data$ENV)
	input_data$GROUP <- as.factor(input_data$GROUP)
	input_data$GENOTYPE <- as.factor(input_data$GENOTYPE)	
	if(dummy){
		input_data[, paste0(colnames(QTL)[2], ".1")] <- NULL
		input_data[, paste0(colnames(QTL2)[2], ".1")] <- NULL
		x <- x[1]
	}
	for(QTL_i in 5:(ncol(hap_vec_matrix) + 4)){
	input_data[, QTL_i] <- factor(input_data[, QTL_i], levels = unique(sort(input_data[, QTL_i])))	
	}	
	print(str(input_data))
	print(head(input_data, n = 2))

#
# model calculations
#

# NULL model
	reml.obj.null <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08,
					   start.values=TRUE)
	
	iv.null <- reml.obj.null$gammas.table

	reml.obj.null <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08)

	iv.null$Value  <- summary(reml.obj.null)$varcomp$component

	reml.obj.null <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08,
					   G.param=iv.null, 
					   R.param=iv.null)

	iv.null$Value <- summary(reml.obj.null)$varcomp$component

# model for each window
for(hap_focus in x){

# in which window lies the focus haplotype
wind <- substr(hap_focus, 1, 13)
wind
region <- paste0(TRAIT, ".", hap_focus)
region

# exchange the number coded focus haplotype by the respective factor coded window in the model
	hap_focus_x <- x
	hap_focus_x[which(hap_focus_x == hap_focus)] <- wind

# specify outputfolder
dir.create(paste("wholeWindModel/", TRAIT, "/", region, sep = ""))
outfolder <- paste("wholeWindModel/", TRAIT, "/", region, sep = "")
outfolder

# Full model
	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(hap_focus_x, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08,
					   na.method.Y = "omit",
					   na.method.X = "omit",
					   start.values=TRUE)
	
	iv.full <- reml.obj$gammas.table

	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(hap_focus_x, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   na.method.Y = "omit",
					   na.method.X = "omit",
					   workspace   = 3e+08)

	iv.full$Value  <- summary(reml.obj)$varcomp$component

	reml.obj <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(hap_focus_x, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   na.method.Y = "omit",
					   na.method.X = "omit",
					   workspace   = 3e+08,
					   G.param=iv.full, 
					   R.param=iv.full)

	iv.full$Value  <- summary(reml.obj)$varcomp$component

if(dummy == FALSE){
# minus1 model
# model only removing the window from the full model
hap_focus_x
hap_focus_x <- hap_focus_x[-which(hap_focus_x == wind)]
hap_focus_x

	reml.obj.minus1 <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(hap_focus_x, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08,
					   na.method.Y = "omit",
					   na.method.X = "omit",
					   start.values=TRUE)
	
	iv.minus1 <- reml.obj.minus1$gammas.table

	reml.obj.minus1 <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(hap_focus_x, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   na.method.Y = "omit",
					   na.method.X = "omit",
					   workspace   = 3e+08)

	iv.minus1$Value  <- summary(reml.obj.minus1)$varcomp$component

	reml.obj.minus1 <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(hap_focus_x, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   na.method.Y = "omit",
					   na.method.X = "omit",
					   workspace   = 3e+08,
					   G.param=iv.minus1, 
					   R.param=iv.minus1)

	iv.minus1$Value  <- summary(reml.obj.minus1)$varcomp$component
}

# onlyWindow model
# model only including the window but no other focus haplotype
	reml.obj.onlyWind <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(wind, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   workspace   = 3e+08,
					   na.method.Y = "omit",
					   na.method.X = "omit",
					   start.values=TRUE)
	
	iv.onlyWind <- reml.obj.onlyWind$gammas.table

	reml.obj.onlyWind <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(wind, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   na.method.Y = "omit",
					   na.method.X = "omit",
					   workspace   = 3e+08)

	iv.onlyWind$Value  <- summary(reml.obj.onlyWind)$varcomp$component

	reml.obj.onlyWind <- asreml(fixed  = formula(paste0("Y ~ ENV + GROUP + 
								ENV:(", paste0(wind, collapse = " + "),")", collapse = "")),
					   random = ~ GENOTYPE,
					   rcov   = ~ at(ENV):units,
					   data   = input_data,
					   maxit  = 500,
					   na.method.Y = "omit",
					   na.method.X = "omit",
					   workspace   = 3e+08,
					   G.param=iv.onlyWind, 
					   R.param=iv.onlyWind)

	iv.onlyWind$Value  <- summary(reml.obj.onlyWind)$varcomp$component

# explained variance
full_vs_null <- (iv.null$Value[1] - iv.full$Value[1]) / iv.null$Value[1]
if(dummy){
full_vs_minusWindow <- NA
} else {
full_vs_minusWindow <- (iv.minus1$Value[1] - iv.full$Value[1]) / iv.minus1$Value[1]
}
onlyWind_vs_null <- (iv.null$Value[1] - iv.onlyWind$Value[1]) / iv.null$Value[1]
expl_var <- c(full_vs_null, full_vs_minusWindow, onlyWind_vs_null)
names(expl_var) <- c("full_vs_null", "full_vs_minusWindow", "onlyWindow_vs_null")
write.table(expl_var, file = paste(outfolder, "/explVar_", region, ".csv", sep = ""), sep = ";", dec = ".", quote = FALSE, row.names = TRUE, col.names = FALSE)

	## get the pvalues of the wald test
	wald.obj <- wald(reml.obj, denDF = "default")
	print(wald.obj)

	# extract effect estimates
	HAPeff <- reml.obj$coefficients$fixed[grep(":wind_",names(reml.obj$coefficients$fixed))]
	seHAPeff <- sqrt(reml.obj$vcoeff$fixed)[1:length(HAPeff)]
	
	# "Confidence intervals" (95%)
	lowerCI95 <- HAPeff - qnorm(0.975, mean = 0, sd = 1) * seHAPeff
	upperCI95 <- HAPeff + qnorm(0.975, mean = 0, sd = 1) * seHAPeff	
	sig <- apply(cbind(lowerCI95, upperCI95), 1, function(x) {findInterval(0, x) != 1})
	sig[which(sig)] <- "*"
	sig[which(sig != "*")] <- "ns"
	
	# effects in phenotypic standard deviations
	HAPeff_sd <- HAPeff
	for(env_i in Envs){
			Y_temp <- read.table(paste(TRAIT, "/", list.files(path = TRAIT, pattern = ".pheno.txt"), sep = ""), stringsAsFactors = FALSE)
			colnames(Y_temp) <- Envs
			Y_temp <- Y_temp[ , env_i]
			SD_temp <- sd(Y_temp, na.rm = TRUE)
			HAPeff_sd[grep(env_i, names(HAPeff_sd))] <- HAPeff_sd[grep(env_i, names(HAPeff_sd))] / SD_temp
			}
	info <- data.frame(HAPeff, lowerCI95, upperCI95, sig, HAPeff_sd)
	info$haplotype <- substr(rownames(info), nchar(rownames(info)), nchar(rownames(info)))
	info$ENV <- substr(rownames(info), 5, 12)
	str(info)
	head(info, n = 20)
	tail(info, n = 20)
	info[grep(wind, rownames(info)),]
	summary(info[grep(wind, rownames(info)),"HAPeff_sd"])
	save(info, file = paste(outfolder, "/info_", region, "_v3.RData", sep = ""))
	info <- info[grep(wind, rownames(info)),]
	write.table(info, file = paste(outfolder, "/info_", region, ".csv", sep =""),
            sep = ";", dec = ".", quote = FALSE, col.names = TRUE, row.names = FALSE)
}

###
###################################################################
