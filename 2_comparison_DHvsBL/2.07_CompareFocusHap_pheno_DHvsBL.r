###################################################################
###################################################################
####
#### for a given focus haplotype, compare the phenotypic performance between
#### DH lines carrying the focus haplotype and breeding lines not carrying the focus haplotype
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

# packages
library(ggplot2)

# arguments
nSNPs <- 10
steps <- 10
p_thresh <- 0.01
FDR <- "15"
minCount <- 3

# choose according to the focus haplotype you want to plot
TRAIT <- "PH_V6"
direction <- "fav"  #"fav" or "unfav"
hap_focus <- "wind_03_01592_6"
leg_pos <- "topleft"

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
TRAIT
direction
hap_focus
leg_pos

# graphical parameters
cex <- 0.6
cex.axis <- 1.6
cex.lab <- 1.6

# permut est function (two.sided)
permut.test_f <- function(x, y, n = 10000, alternative = "two.sided"){
  # remove NAs
  x <- na.omit(x)
  y <- na.omit(y)
  # count entries per vector
  nx <- length(x)
  ny <- length(y)
  
  if((nx != 0) & (ny != 0)){
  # calculate parameter of interest (here difference in means)
  stat <- mean(x) - mean(y)
  
  # generate permuted samples (in columns)
  perm_matrix <- replicate(n, sample(c(x, y)))
  stat_f <- function(x, nx){
    mean(x[1:nx]) - mean(x[(nx+1):length(x)])
  }
  # apply test to permuted samples to generate null distribution
  null_distr <- apply(perm_matrix, 2, stat_f, nx = nx)
  
  # calculate p-value
  n_NULL.larger <- length(which(null_distr > stat))
  n_NULL.lower <- length(which(null_distr < stat))
  n_NULL.equal <- length(which(null_distr == stat))
  if (alternative == "two.sided") 
    p_value <- (2 * (min(n_NULL.larger, n_NULL.lower) + 0.5 * n_NULL.equal)) / (n+1)
  if (alternative == "less") 
    p_value <- (n_NULL.lower + 0.5 * n_NULL.equal) / (n+1)
  if (alternative == "greater") 
    p_value <- (n_NULL.larger + 0.5 * n_NULL.equal) / (n+1)
  p_value <- min(p_value, 1)
  
  } else {
  stat <- NA
  null_distr <- NA
  p_value <- NA
  
  }
  # generate output object
  perm_out <- list(stat = stat,
                   null_distr = null_distr,
                   p_value = p_value,
                   original_x = x,
                   original_y = y,
                   nx = nx,
                   ny = ny
  )
  return(perm_out)
}

# generate letter vector for naming haplotype variants
letters_MM <- c(letters[1:26], paste(rep(letters[1:26], rep(26,26)), rep(letters[1:26], 26), sep = ""))
letters_MM <- letters_MM[-1]
letters_MM

#
# define colors
#
col2rgb_MM <- function(x, alpha){
  rgb(x[1],x[2],x[3], alpha = alpha)
}

temp <- (col2rgb("lightblue")/255)[,1]
temp
temp <- c(0.6,0.6,0.9)
BL.2_edge <- col2rgb_MM(temp, alpha = 1)
BL.2_fill <- col2rgb_MM(temp, alpha = 0.25)
BL.2_fill2 <- col2rgb_MM(temp, alpha = 0.65)

temp <- (col2rgb("darkblue")/255)[,1]
temp
temp <- c(0,0,0.6)
BL.1_edge <- col2rgb_MM(temp, alpha = 1)
BL.1_fill <- col2rgb_MM(temp, alpha = 0.25)
BL.1_fill2 <- col2rgb_MM(temp, alpha = 0.65)

temp <- (col2rgb("green2")/255)[,1]
temp
temp <- c(0,0.85,0)
LR.2_edge <- col2rgb_MM(temp, alpha = 1)
LR.2_fill <- col2rgb_MM(temp, alpha = 0.25)
LR.2_fill2 <- col2rgb_MM(temp, alpha = 0.65)

temp <- (col2rgb("darkgreen")/255)[,1]
temp
temp <- c(0,0.4,0)
LR.1_edge <- col2rgb_MM(temp, alpha = 1)
LR.1_fill <- col2rgb_MM(temp, alpha = 0.25)
LR.1_fill2 <- col2rgb_MM(temp, alpha = 0.65)


############################################################################################################################################################################################################
############################################################################################################################################################################################################
##############
##############

wind <- substr(hap_focus, 1, 13)
wind
region <- paste0(TRAIT, ".", hap_focus)
region

# input folder
infolder <- paste(TRAIT, sep = "")
infolder

# output folder
outfolder <- paste("wholeWindModel/", TRAIT, "/", region, sep = "")
outfolder

#################################################################################################################
# load DH data
load(paste("geno_InfoList_DH_m", nSNPs, "s", steps, ".RData", sep = ""))
str(InfoList)
str(geno)

#
# check phenotypic distribution
#
  # BLUEs for 2017
  phenoData <- read.table("Mayer_et_al_phenotypes_BLUEs_acrossSixLocations2017.txt.gz", header = TRUE, stringsAsFactors = FALSE)
  str(phenoData)
  
  pheno <- phenoData[ , TRAIT]
  names(pheno) <- phenoData$Genotype
  
  # load haplotype info of BL
  load(paste("comp_BL_DHs/focusHaps/", TRAIT, "/info_BL_", direction, "_", TRAIT, ".RData", sep = ""))
  if(direction == "unfav")
  geno_BL_hap <- geno_BL_unfav
  if(direction == "fav")
  geno_BL_hap <- geno_BL_fav
  rownames(geno_BL_hap) <- substr(rownames(geno_BL_hap), 7, nchar(rownames(geno_BL_hap)))
  phen_checks <- pheno[which((substr(names(pheno), 1, 3) != "DH_") & (substr(names(pheno), 1, 3) != "LR_") & (substr(names(pheno), 1, 4) != "F353"))]
  phen_checks <- phen_checks[which(is.na(phen_checks) == FALSE)]
  phen_checks_all <- phen_checks
  hap_temp <- geno_BL_hap[names(phen_checks), hap_focus]
  phen_checks <- phen_checks[which(hap_temp == 0)]
  print(length(phen_checks))
  phen_checks_all
  phen_checks

# test difference phen_checks_without vs phen_checks_with (if at least 2 individuals are carriers)
phen_checks_without <- phen_checks
phen_checks_with <- phen_checks_all[which(phen_checks_all %in% phen_checks == FALSE)]
if(length(phen_checks_with) > 1){
permtest_BL <- permut.test_f(phen_checks_without, phen_checks_with)
print(str(permtest_BL))
}

LR_case <- NULL
xlim_all <- matrix(NA, nrow = 3, ncol = 2)
rownames(xlim_all) <- c("KE", "LL", "PE")
ylim_all <- matrix(NA, nrow = 3, ncol = 2)
rownames(ylim_all) <- c("KE", "LL", "PE")

if(substr(TRAIT,1,3) == "PH_"){
unit_trait <- "cm"
} else {
	if(TRAIT %in% c("MF", "FF")){
	unit_trait <- "days"
	} else {
	unit_trait <- "score"	
	}
}

perm_test <- matrix(NA, nrow = 3, ncol = 2)
colnames(perm_test) <- c("all", "focus")
rownames(perm_test) <- c("KE", "LL", "PE")
n_PerGroup <- matrix(NA, nrow = 4, ncol = 2)
colnames(n_PerGroup) <- c("all", "focus")
rownames(n_PerGroup) <- c("KE", "LL", "PE", "BL")

for(LR in c("KE", "LL", "PE")){
	  # LR
	  phen_LRall <- pheno[which(substr(names(pheno), 1, 5) == paste("DH_", LR, sep = ""))]
	  phen_LRall <- phen_LRall[which(is.na(phen_LRall) == FALSE)]
	  hap_temp <- geno[names(phen_LRall), hap_focus]
	  phen_LR <- phen_LRall[which(hap_temp == 2)]
	  
	  print(LR)
	  print(length(phen_LR))
	  
	  perm_test[LR, "all"] <- permut.test_f(phen_LRall, phen_checks_all)$p_value
	  n_PerGroup[LR, "all"] <- length(phen_LRall)
	  n_PerGroup["BL", "all"] <- length(phen_checks_all)
	  
	  if(length(phen_LR) > 1){
		  LR_case <- c(LR_case, LR)
		  # graphical settings
		  cex <- 0.6
		  cex.axis <- 1.4
		  cex.lab <- 1.4
		  xmin <- min(c(phen_LR, phen_LRall, phen_checks, phen_checks_all))
		  xmax <- max(c(phen_LR, phen_LRall, phen_checks, phen_checks_all))
		  xlim1 <- c(xmin - (xmax - xmin) * 0.11,  xmax + (xmax - xmin) * 0.11)
		  xlim_all[LR, ] <- xlim1
		  perm_test[LR, "focus"] <- permut.test_f(phen_LR, phen_checks)$p_value
		  n_PerGroup[LR, "focus"] <- length(phen_LR)
		  n_PerGroup["BL", "focus"] <- length(phen_checks)

		  png(paste(outfolder, "/PhenoDensityAcross2017_", region, "_", LR, ".png", sep =""), width = 2400, height = 2050, res = 300)
		  par(mar = c(3.5, 3.5, 0.1, 0.4) + 0.1, mgp = c(2.5, 1, 0))
		  
		  d_BL <- density(phen_checks)
		  d_BLall <- density(phen_checks_all)
		  d_LR <- density(phen_LR)
		  d_LRall <- density(phen_LRall)
		  ymax <- max(c(d_BLall$y, d_BL$y, d_LRall$y, d_LR$y)*1.1, na.rm = TRUE)
		  ylim_all[LR, ] <- c(0, ymax)
		  
		  plot(-1000, -1000, xlim = xlim1, lwd = 2.5, ylim = c(0, ymax), main = "", ylab = "Density", xlab = paste(TRAIT, " (", unit_trait, ")", sep =""),
			   cex.lab = cex.lab, cex.axis = cex.axis)
		  polygon(x = c(- 1000, - 1000, + 1000, + 1000), y = c(-1000, + 1000, + 1000, -1000), col = "white", border = "white", lwd = 2)
		  
		  polygon(d_BLall, col = BL.1_fill, border=NA) 
		  polygon(d_LRall, col = LR.1_fill, border=NA)
		  polygon(d_LR, col = LR.2_fill, border=NA)
		  
		  col_leg <- c(BL.1_edge, LR.1_edge, LR.2_edge)
		  col_leg2 <- c(BL.1_fill, LR.1_fill, LR.2_fill)

		  legendlabels <- c(paste0("BL_all"), paste0(LR, "_all"), paste0(LR, "_Focus"))
		  if(length(phen_checks) < 14){
			polygon(d_BL, col = BL.2_fill, border=NA)
		    col_leg <- c(col_leg, BL.2_edge)
		    col_leg2 <- c(col_leg2, BL.2_fill)
			legendlabels <- c(legendlabels, paste("BL_noFocus"))
			legendlabels <- legendlabels[c(1,4,2,3)]
			col_leg <- col_leg[c(1,4,2,3)]
			col_leg2 <- col_leg2[c(1,4,2,3)]
		  }

		  points(d_BLall, col = BL.1_edge, xlim = xlim1, lwd = 2.5, type = "l")
		  points(d_LRall, col = LR.1_edge, xlim = xlim1, lwd = 2.5, type = "l")
		  points(d_LR, col = LR.2_edge, xlim = xlim1, lwd = 2.5, type = "l")
		  if(length(phen_checks) < 14){
			points(d_BL, col = BL.2_edge, xlim = xlim1, lwd = 2.5, type = "l")
		  }
		  abline(v = mean(phen_checks_all), lty = 2, lwd = 2, col = BL.1_edge)
		  abline(v = mean(phen_LRall), lty = 2, lwd = 2, col = LR.1_edge)
		  abline(v = mean(phen_LR), lty = 2, lwd = 2, col = LR.2_edge)		  
		  if(length(phen_checks) < 14){
			points(d_BL, col = BL.2_edge, xlim = xlim1, lwd = 2.5, type = "l")
		    abline(v = mean(phen_checks), lty = 2, lwd = 2, col = BL.2_edge)
		  }

			legend(leg_pos, pch = 22, col = col_leg,
			pt.bg = col_leg2,
			pt.lwd = 2, pt.cex = 2,
			text.col = col_leg, cex = 1.2,
			legend = legendlabels)
		  
		  dev.off()
	}
}
write.table(perm_test, file = paste(outfolder, "/perm_test_DH.BL_Across2017_", region, "_", LR, ".csv", sep = ""), sep = ";", dec = ".", quote = FALSE, col.names = NA, row.names = TRUE)
write.table(n_PerGroup, file = paste(outfolder, "/n_PerGroup_DH.BL_Across2017_", region, "_", LR, ".csv", sep = ""), sep = ";", dec = ".", quote = FALSE, col.names = NA, row.names = TRUE)



#
# also check phenotypic data for individual environments
#
pheno_perEnv <- read.table("Mayer_et_al_phenotypes_BLUEs_perEnvironment.txt.gz", header = TRUE, stringsAsFactors = FALSE)
str(pheno_perEnv)

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
  
  pheno_all <- GWAS_table
  str(pheno_all)

#
# here only for 2017 data
#
pheno_all <- pheno_all[ , c(which((colnames(pheno_all) == "Genotype")),grep("2017", colnames(pheno_all)))]
str(pheno_all)

LR <- LR_case[1]

Value <- NULL
Group <- NULL
Environment <- NULL

p_values_store <- NULL
env_store <- NULL


for(env_i in colnames(pheno_all)[-1]){
print(env_i)

  pheno <- pheno_all[, env_i]
  names(pheno) <- pheno_all$Genotype
  phen_checks <- pheno[which((substr(names(pheno), 1, 3) != "DH_") & (substr(names(pheno), 1, 3) != "LR_") & (substr(names(pheno), 1, 4) != "F353"))]
  phen_checks <- phen_checks[which(is.na(phen_checks) == FALSE)]
  phen_checks_all <- phen_checks
  hap_temp <- geno_BL_hap[names(phen_checks), hap_focus]
  phen_checks <- phen_checks[which(hap_temp == 0)]
 
  print("BL")
  print(length(phen_checks))

	  # LR
	  phen_LRall <- pheno[which(substr(names(pheno), 1, 5) == paste("DH_", LR, sep = ""))]
	  phen_LRall <- phen_LRall[which(is.na(phen_LRall) == FALSE)]
	  hap_temp <- geno[names(phen_LRall), hap_focus]
	  phen_LR <- phen_LRall[which(hap_temp == 2)]
	  
	  print(LR)
	  print(length(phen_LR))
	  
		  permtest <- permut.test_f(phen_checks, phen_LR)
		  p_values_store <- c(p_values_store, permtest$p_value)
		  env_store <- c(env_i, env_store)
		  
	 Value <- c(Value, phen_LR, phen_checks) 
	 Group <- c(Group, rep(paste0(LR, "_Focus"), length(phen_LR)), rep("BL_noFocus", length(phen_checks))) 
	 Environment <- c(Environment, rep(env_i, length(phen_LR)), rep(env_i, length(phen_checks)))

}
names(p_values_store) <- colnames(pheno_all)[-1]
p_values_store

DF <- data.frame(Value = Value,
				 Group = Group,
				 Environment = Environment)
lev <- c("2017.BBG", "2017.EIN", "2017.GOL", "2017.OLI", "2017.ROG", "2017.TOM",
         "2018.EIN", "2018.GOL", "2018.KLW", "2018.ROG", "2018.TOM")
lev <- lev[which(lev %in% colnames(pheno_all))]
DF$Environment <- factor(DF$Environment, levels = lev)
DF$Group <- factor(DF$Group, levels = c(paste0(LR, "_Focus"), "BL_noFocus"))

if(length(phen_checks_all) == length(phen_checks)){
	col_manual <- c(LR.2_fill2, BL.1_fill2)
} else {
	col_manual <- c(LR.2_fill2, BL.2_fill2)
}

ymax <- max(DF$Value) *1.01

		  png(paste(outfolder, "/PhenoPerEnv_", region, "_", LR, ".png", sep =""), width = 3600 * (length(p_values_store) / 6), height = 1800, res = 300)
		  par(mar = c(3.5, 3.5, 0.1, 0.4) + 0.1, mgp = c(2.5, 1, 0))
			ggplot(DF, aes(Environment, Value, fill = Group)) +
			  geom_boxplot(outlier.shape = NA) +
			  scale_fill_manual(values = col_manual) +
			  geom_point(position = position_jitterdodge(), aes(alpha = 0.2), show.legend = FALSE) +
			  theme_bw() +
			  theme(axis.text.x = element_text(size = 14),
			        axis.text.y = element_text(size = 14),
					axis.title=element_text(size=14),
					legend.title = element_text(color = "white"),
					legend.text = element_text(size = 14)) +
			  xlab("") +
			  ylab(paste(TRAIT, " (", unit_trait, ")", sep ="")) +
			  geom_text(x=1, y=ymax, label = round(p_values_store[1], digits = 4)) + 
			  geom_text(x=2, y=ymax, label = round(p_values_store[2], digits = 4)) + 
			  geom_text(x=3, y=ymax, label = round(p_values_store[3], digits = 4)) + 
			  geom_text(x=4, y=ymax, label = round(p_values_store[4], digits = 4)) + 
			  geom_text(x=5, y=ymax, label = round(p_values_store[5], digits = 4)) + 
			  geom_text(x=6, y=ymax, label = round(p_values_store[6], digits = 4)) + 
			  geom_text(x=7, y=ymax, label = round(p_values_store[7], digits = 4)) + 
			  geom_text(x=8, y=ymax, label = round(p_values_store[8], digits = 4)) + 
			  geom_text(x=9, y=ymax, label = round(p_values_store[9], digits = 4)) + 
			  geom_text(x=10, y=ymax, label = round(p_values_store[10], digits = 4)) + 
			  geom_text(x=11, y=ymax, label = round(p_values_store[11], digits = 4))
		  dev.off()


if(length(LR_case) > 1)
LR <- LR_case[2]

Value <- NULL
Group <- NULL
Environment <- NULL

p_values_store <- NULL
env_store <- NULL


for(env_i in colnames(pheno_all)[-1]){
print(env_i)

  pheno <- pheno_all[, env_i]
  names(pheno) <- pheno_all$Genotype
  phen_checks <- pheno[which((substr(names(pheno), 1, 3) != "DH_") & (substr(names(pheno), 1, 3) != "LR_") & (substr(names(pheno), 1, 4) != "F353"))]
  phen_checks <- phen_checks[which(is.na(phen_checks) == FALSE)]
  phen_checks_all <- phen_checks
  hap_temp <- geno_BL_hap[names(phen_checks), hap_focus]
  phen_checks <- phen_checks[which(hap_temp == 0)]
 
  print("BL")
  print(length(phen_checks))

	  # LR
	  phen_LRall <- pheno[which(substr(names(pheno), 1, 5) == paste("DH_", LR, sep = ""))]
	  phen_LRall <- phen_LRall[which(is.na(phen_LRall) == FALSE)]
	  hap_temp <- geno[names(phen_LRall), hap_focus]
	  phen_LR <- phen_LRall[which(hap_temp == 2)]
	  
	  print(LR)
	  print(length(phen_LR))
	  
		  permtest <- permut.test_f(phen_checks, phen_LR)
		  p_values_store <- c(p_values_store, permtest$p_value)
		  env_store <- c(env_i, env_store)
		  
	 Value <- c(Value, phen_LR, phen_checks) 
	 Group <- c(Group, rep(paste0(LR, "_Focus"), length(phen_LR)), rep("BL_noFocus", length(phen_checks))) 
	 Environment <- c(Environment, rep(env_i, length(phen_LR)), rep(env_i, length(phen_checks)))

}
names(p_values_store) <- colnames(pheno_all)[-1]
p_values_store

DF <- data.frame(Value = Value,
				 Group = Group,
				 Environment = Environment)
lev <- c("2017.BBG", "2017.EIN", "2017.GOL", "2017.OLI", "2017.ROG", "2017.TOM",
         "2018.EIN", "2018.GOL", "2018.KLW", "2018.ROG", "2018.TOM")
lev <- lev[which(lev %in% colnames(pheno_all))]
DF$Environment <- factor(DF$Environment, levels = lev)
DF$Group <- factor(DF$Group, levels = c(paste0(LR, "_Focus"), "BL_noFocus"))

if(length(phen_checks_all) == length(phen_checks)){
	col_manual <- c(LR.2_fill2, BL.1_fill2)
} else {
	col_manual <- c(LR.2_fill2, BL.2_fill2)
}

ymax <- max(DF$Value) *1.01

		  png(paste(outfolder, "/PhenoPerEnv_", region, "_", LR, ".png", sep =""), width = 3600 * (length(p_values_store) / 6), height = 1800, res = 300)
		  par(mar = c(3.5, 3.5, 0.1, 0.4) + 0.1, mgp = c(2.5, 1, 0))
			ggplot(DF, aes(Environment, Value, fill = Group)) +
			  geom_boxplot(outlier.shape = NA) +
			  scale_fill_manual(values = col_manual) +
			  geom_point(position = position_jitterdodge(), aes(alpha = 0.2), show.legend = FALSE) +
			  theme_bw() +
			  theme(axis.text.x = element_text(size = 14),
			        axis.text.y = element_text(size = 14),
					axis.title=element_text(size=14),
					legend.title = element_text(color = "white"),
					legend.text = element_text(size = 14)) +
			  xlab("") +
			  ylab(paste(TRAIT, " (", unit_trait, ")", sep ="")) +
			  geom_text(x=1, y=ymax, label = round(p_values_store[1], digits = 4)) + 
			  geom_text(x=2, y=ymax, label = round(p_values_store[2], digits = 4)) + 
			  geom_text(x=3, y=ymax, label = round(p_values_store[3], digits = 4)) + 
			  geom_text(x=4, y=ymax, label = round(p_values_store[4], digits = 4)) + 
			  geom_text(x=5, y=ymax, label = round(p_values_store[5], digits = 4)) + 
			  geom_text(x=6, y=ymax, label = round(p_values_store[6], digits = 4)) + 
			  geom_text(x=7, y=ymax, label = round(p_values_store[7], digits = 4)) + 
			  geom_text(x=8, y=ymax, label = round(p_values_store[8], digits = 4)) + 
			  geom_text(x=9, y=ymax, label = round(p_values_store[9], digits = 4)) + 
			  geom_text(x=10, y=ymax, label = round(p_values_store[10], digits = 4)) + 
			  geom_text(x=11, y=ymax, label = round(p_values_store[11], digits = 4))
		  dev.off()


if(length(LR_case) > 2)
LR <- LR_case[3]

Value <- NULL
Group <- NULL
Environment <- NULL

p_values_store <- NULL
env_store <- NULL


for(env_i in colnames(pheno_all)[-1]){
print(env_i)

  pheno <- pheno_all[, env_i]
  names(pheno) <- pheno_all$Genotype
  phen_checks <- pheno[which((substr(names(pheno), 1, 3) != "DH_") & (substr(names(pheno), 1, 3) != "LR_") & (substr(names(pheno), 1, 4) != "F353"))]
  phen_checks <- phen_checks[which(is.na(phen_checks) == FALSE)]
  phen_checks_all <- phen_checks
  hap_temp <- geno_BL_hap[names(phen_checks), hap_focus]
  phen_checks <- phen_checks[which(hap_temp == 0)]
 
  print("BL")
  print(length(phen_checks))

	  # LR
	  phen_LRall <- pheno[which(substr(names(pheno), 1, 5) == paste("DH_", LR, sep = ""))]
	  phen_LRall <- phen_LRall[which(is.na(phen_LRall) == FALSE)]
	  hap_temp <- geno[names(phen_LRall), hap_focus]
	  phen_LR <- phen_LRall[which(hap_temp == 2)]
	  
	  print(LR)
	  print(length(phen_LR))
	  
		  permtest <- permut.test_f(phen_checks, phen_LR)
		  p_values_store <- c(p_values_store, permtest$p_value)
		  env_store <- c(env_i, env_store)
		  
	 Value <- c(Value, phen_LR, phen_checks) 
	 Group <- c(Group, rep(paste0(LR, "_Focus"), length(phen_LR)), rep("BL_noFocus", length(phen_checks))) 
	 Environment <- c(Environment, rep(env_i, length(phen_LR)), rep(env_i, length(phen_checks)))

}
names(p_values_store) <- colnames(pheno_all)[-1]
p_values_store

DF <- data.frame(Value = Value,
				 Group = Group,
				 Environment = Environment)
lev <- c("2017.BBG", "2017.EIN", "2017.GOL", "2017.OLI", "2017.ROG", "2017.TOM",
         "2018.EIN", "2018.GOL", "2018.KLW", "2018.ROG", "2018.TOM")
lev <- lev[which(lev %in% colnames(pheno_all))]
DF$Environment <- factor(DF$Environment, levels = lev)
DF$Group <- factor(DF$Group, levels = c(paste0(LR, "_Focus"), "BL_noFocus"))

if(length(phen_checks_all) == length(phen_checks)){
	col_manual <- c(LR.2_fill2, BL.1_fill2)
} else {
	col_manual <- c(LR.2_fill2, BL.2_fill2)
}

ymax <- max(DF$Value) *1.01

		  png(paste(outfolder, "/PhenoPerEnv_", region, "_", LR, ".png", sep =""), width = 3600 * (length(p_values_store) / 6), height = 1800, res = 300)
		  par(mar = c(3.5, 3.5, 0.1, 0.4) + 0.1, mgp = c(2.5, 1, 0))
			ggplot(DF, aes(Environment, Value, fill = Group)) +
			  geom_boxplot(outlier.shape = NA) +
			  scale_fill_manual(values = col_manual) +
			  geom_point(position = position_jitterdodge(), aes(alpha = 0.2), show.legend = FALSE) +
			  theme_bw() +
			  theme(axis.text.x = element_text(size = 14),
			        axis.text.y = element_text(size = 14),
					axis.title=element_text(size=14),
					legend.title = element_text(color = "white"),
					legend.text = element_text(size = 14)) +
			  xlab("") +
			  ylab(paste(TRAIT, " (", unit_trait, ")", sep ="")) +
			  geom_text(x=1, y=ymax, label = round(p_values_store[1], digits = 4)) + 
			  geom_text(x=2, y=ymax, label = round(p_values_store[2], digits = 4)) + 
			  geom_text(x=3, y=ymax, label = round(p_values_store[3], digits = 4)) + 
			  geom_text(x=4, y=ymax, label = round(p_values_store[4], digits = 4)) + 
			  geom_text(x=5, y=ymax, label = round(p_values_store[5], digits = 4)) + 
			  geom_text(x=6, y=ymax, label = round(p_values_store[6], digits = 4)) + 
			  geom_text(x=7, y=ymax, label = round(p_values_store[7], digits = 4)) + 
			  geom_text(x=8, y=ymax, label = round(p_values_store[8], digits = 4)) + 
			  geom_text(x=9, y=ymax, label = round(p_values_store[9], digits = 4)) + 
			  geom_text(x=10, y=ymax, label = round(p_values_store[10], digits = 4)) + 
			  geom_text(x=11, y=ymax, label = round(p_values_store[11], digits = 4))
		  dev.off()



##############
##############
##############
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################




