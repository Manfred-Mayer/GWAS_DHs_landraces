###################################################################
###################################################################
####
#### generate barplots for haplotypes with significant effects on multiple traits
#### showing the proportion of genetic variance explained for each trait
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
p_thresh <- 0.01
FDR <- "15"

# graphical parameters
cex <- 0.6
cex.axis <- 1.6
cex.lab <- 1.6

###################################################################################################
###################################################################################################
###################################################################################################

dir.create(paste("TraitCorBivar/barplots_explVar", sep = ""))
outfolder <- paste("TraitCorBivar/barplots_explVar", sep = "")
outfolder

traits <- c("EV_V4", "EV_V6", "PH_V4", "PH_V6",
			"PH_final",
			"FF", "MF",
			"LO", "TILL")

traits1 <- c("PH_V4", "PH_V6", "EV_V4", "EV_V6")
traits1

traits2 <- c("PH_final",
			"FF", "MF",
			"LO", "TILL")
traits2

TRAIT1 <- "PH_V6"

all_fav_withMultEff <- NULL
all_unfav_withMultEff <- NULL

for(TRAIT1 in traits1){
print(TRAIT1)

# generate letter vector for naming haplotype variants
LETTERS_MM <- c(LETTERS[1:26], paste(rep(LETTERS[1:26], rep(26,26)), rep(LETTERS[1:26], 26), sep = ""))
if(TRAIT1 == "PH_V6"){
LETTERS_MM <- c(LETTERS[5:26], paste(rep(LETTERS[1:26], rep(26,26)), rep(LETTERS[1:26], 26), sep = ""))
}
LETTERS_MM

direction <- "fav"

explVar <- NULL
sig <- NULL
eff_sign <- NULL
	for(TRAIT2 in traits2){
		print(TRAIT1)
		print(TRAIT2)
		covariation <- ifelse(TRAIT2 %in% c("PH_final", "LO", "TILL"), "pos", "neg")
		print(covariation)

		if(file.exists(paste("TraitCorBivar/", TRAIT1, ".", TRAIT2, sep = ""))){
			infolder <- paste("TraitCorBivar/", TRAIT1, ".", TRAIT2, sep = "")
			info <- read.table(paste(infolder, "/info_cor_", TRAIT1, "_", TRAIT2, "_main.csv", sep =""), sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)			
		} else {			
			infolder <- paste("TraitCorBivar/", TRAIT2, ".", TRAIT1, sep = "")
			info <- read.table(paste(infolder, "/info_cor_", TRAIT2, "_", TRAIT1, "_main.csv", sep =""), sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)
		}

		# filter only for haplotypes identified for TRAIT1
		# for the respective direction (fav or unfav)
		load(paste("FavUnfav_Stability/", TRAIT1, "/haps_all_", TRAIT1, ".RData", sep = ""))
		if(direction == "fav"){
			dir_haps <- haps_Fav_info$qtl_i
		} else {
			dir_haps <- haps_Unfav_info$qtl_i		
		}
		print(str(dir_haps))
		info <- info[which(info$haplotype %in% dir_haps), ]
		info <- info[which(is.na(info$eff_y1) == FALSE), ]
		explVar_i <- info[ , c(16, 18)]
		rownames(explVar_i) <- info$haplotype
		sig_i <- info[ , c(21, 22)]
		colnames(sig_i) <- gsub("VARexpl_", "sig_", colnames(explVar_i))
		rownames(sig_i) <- rownames(explVar_i)
		eff_sign_i <- apply(info[ , c(4, 7)] ,2, sign)
		colnames(eff_sign_i) <- gsub("VARexpl_", "effSign_", colnames(explVar_i))
		rownames(eff_sign_i) <- rownames(explVar_i)
				
		if(is.null(explVar)){
		  explVar <- explVar_i
		  sig <- sig_i
		  eff_sign <- eff_sign_i
		} else {
		  explVar <- cbind(explVar, explVar_i)
		  sig <- cbind(sig, sig_i)
		  eff_sign <- cbind(eff_sign, eff_sign_i)
		}
}

# extract haplotypes which have a significant effect on any of those other traits
sig_other <- sig[ , -grep(TRAIT1, colnames(sig))]
sig_other <- sig_other[apply(sig_other, 1, function(x) {any(x != 0)}), ]
sig_other
sig_other_haps <- rownames(sig_other)

sig_TRAIT1 <- sig[sig_other_haps, grep(TRAIT1, colnames(sig))]
explVar_TRAIT1 <- explVar[sig_other_haps, grep(TRAIT1, colnames(sig))]
eff_sign_TRAIT1 <- eff_sign[sig_other_haps, grep(TRAIT1, colnames(sig))]

explVar_other <- explVar[sig_other_haps, -grep(TRAIT1, colnames(sig))]
eff_sign_other <- eff_sign[sig_other_haps, -grep(TRAIT1, colnames(sig))]

eff_sign_TRAIT1
apply(explVar_TRAIT1, 1, summary)
VARexpl_TRAIT1 <- apply(explVar_TRAIT1, 1, max)

explVar_other
explVar_other[sig_other == 0] <- 0
explVar_other

explVar <- cbind(VARexpl_TRAIT1, explVar_other)
colnames(explVar) <- gsub("VARexpl_", "", colnames(explVar))

colnames(eff_sign_other) <- gsub("effSign_", "", colnames(eff_sign_other))
# now check if it is an unfavorable effect combination
explVar
deleted_effect <- list()
colnames(sig_other) <- gsub("sig_", "", colnames(sig_other))
for(i in colnames(explVar)[-1]){
	if(direction == "fav"){
		if(i %in% c("FF", "MF")){
		deleted_effect_i <- rownames(explVar)[which((eff_sign_other[,i] == 1) & (sig_other[ , i] == 1))]
		explVar[which(eff_sign_other[,i] == 1),i] <- 0
		} else {
		deleted_effect_i <- rownames(explVar)[which((eff_sign_other[,i] == -1) & (sig_other[ , i] == 1))]
		explVar[which(eff_sign_other[,i] == -1),i] <- 0
		}
	}
	if(direction == "unfav"){
		if(i %in% c("FF", "MF")){
		deleted_effect_i <- rownames(explVar)[which((eff_sign_other[,i] == -1) & (sig_other[ , i] == 1))]
		explVar[which(eff_sign_other[,i] == -1),i] <- 0
		} else {
		deleted_effect_i <- rownames(explVar)[which((eff_sign_other[,i] == 1) & (sig_other[ , i] == 1))]
		explVar[which(eff_sign_other[,i] == 1),i] <- 0
		}
	}
	deleted_effect[[i]] <- deleted_effect_i
}
explVar
colnames(explVar)[which(colnames(explVar) == "TRAIT1")] <- TRAIT1
colnames(explVar)[which(colnames(explVar) == "FF")] <- "FF"
colnames(explVar)[which(colnames(explVar) == "MF")] <- "MF"
colnames(explVar)[which(colnames(explVar) == "LO")] <- "LO"
explVar
explVar <- explVar[apply(explVar[, -1], 1, function(x) {any(x != 0)}), ]
explVar

write.table(explVar, paste(outfolder, "/varExpl_onlySig_", TRAIT1, "vsOthers_", direction, ".csv", sep = ""),
			sep = ";", dec = ".", quote = FALSE, row.names = TRUE, col.names = NA)
save(deleted_effect, file = paste(outfolder, "/deleted_FavEffect_", TRAIT1, "vsOthers_", direction, ".RData", sep = ""))

sig_TRAIT1[rownames(explVar), ]

all_fav_withMultEff <- c(all_fav_withMultEff, rownames(explVar))

# create barplot
cols <- colorRampPalette(c("white",rgb(0,0,0.35, alpha = 1)))(20)
Varexpl <- explVar[order(explVar[, TRAIT1], decreasing = TRUE), ]
Varexpl <- t(Varexpl)
						png(paste(outfolder, "/barplot_varExpl_", TRAIT1, "vsOthers_", direction, ".png", sep = ""), width = 300 + 2400 * (ncol(Varexpl) / 9), height = 2100, res = 300)
						par(mar = c(3.8, 3.9, 0, 0) + 0.1, mgp = c(2.5,0.75,0))
						barplot(height = Varexpl,
								beside = TRUE,
								col = c(cols[20], cols[16], cols[13], cols[10], cols[7], cols[4]),
								names.arg = LETTERS_MM[1:ncol(Varexpl)],
								legend.text = TRUE,
								xlab = paste0("Favorable haplotypes with effect on ", TRAIT1, " and other traits"),
								ylab = "Prop. of genetic variance expl.",
								args.legend = list(x =  "top", cex = 1.5),
								cex.lab = 1.25,
								cex.axis = 1.25,
								cex.names = 1.25
								)
						dev.off()
LETTERS_MM <- LETTERS_MM[-(1:ncol(Varexpl))]


direction <- "unfav"

explVar <- NULL
sig <- NULL
eff_sign <- NULL
	for(TRAIT2 in traits2){
		print(TRAIT1)
		print(TRAIT2)
		covariation <- ifelse(TRAIT2 %in% c("PH_final", "LO", "TILL"), "pos", "neg")
		print(covariation)

		if(file.exists(paste("TraitCorBivar/", TRAIT1, ".", TRAIT2, sep = ""))){
			infolder <- paste("TraitCorBivar/", TRAIT1, ".", TRAIT2, sep = "")
			info <- read.table(paste(infolder, "/info_cor_", TRAIT1, "_", TRAIT2, "_main.csv", sep =""), sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)			
		} else {			
			infolder <- paste("TraitCorBivar/", TRAIT2, ".", TRAIT1, sep = "")
			info <- read.table(paste(infolder, "/info_cor_", TRAIT2, "_", TRAIT1, "_main.csv", sep =""), sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)	
		}

		# filter only for haplotypes identified for TRAIT1
		# for the respective direction (fav or unfav)
		load(paste("FavUnfav_Stability/", TRAIT1, "/haps_all_", TRAIT1, ".RData", sep = ""))
		if(direction == "fav"){
			dir_haps <- haps_Fav_info$qtl_i
		} else {
			dir_haps <- haps_Unfav_info$qtl_i		
		}
		print(str(dir_haps))
		info <- info[which(info$haplotype %in% dir_haps), ]
		info <- info[which(is.na(info$eff_y1) == FALSE), ]
		explVar_i <- info[ , c(16, 18)]
		rownames(explVar_i) <- info$haplotype
		sig_i <- info[ , c(21, 22)]
		colnames(sig_i) <- gsub("VARexpl_", "sig_", colnames(explVar_i))
		rownames(sig_i) <- rownames(explVar_i)
		eff_sign_i <- apply(info[ , c(4, 7)] ,2, sign)
		colnames(eff_sign_i) <- gsub("VARexpl_", "effSign_", colnames(explVar_i))
		rownames(eff_sign_i) <- rownames(explVar_i)
				
		if(is.null(explVar)){
		  explVar <- explVar_i
		  sig <- sig_i
		  eff_sign <- eff_sign_i
		} else {
		  explVar <- cbind(explVar, explVar_i)
		  sig <- cbind(sig, sig_i)
		  eff_sign <- cbind(eff_sign, eff_sign_i)
		}
}

# extract haplotypes which have a significant effect on any of those other traits
sig_other <- sig[ , -grep(TRAIT1, colnames(sig))]
sig_other <- sig_other[apply(sig_other, 1, function(x) {any(x != 0)}), ]
sig_other
sig_other_haps <- rownames(sig_other)

sig_TRAIT1 <- sig[sig_other_haps, grep(TRAIT1, colnames(sig))]
explVar_TRAIT1 <- explVar[sig_other_haps, grep(TRAIT1, colnames(sig))]
eff_sign_TRAIT1 <- eff_sign[sig_other_haps, grep(TRAIT1, colnames(sig))]

explVar_other <- explVar[sig_other_haps, -grep(TRAIT1, colnames(sig))]
eff_sign_other <- eff_sign[sig_other_haps, -grep(TRAIT1, colnames(sig))]

eff_sign_TRAIT1
apply(explVar_TRAIT1, 1, summary)
VARexpl_TRAIT1 <- apply(explVar_TRAIT1, 1, max)

explVar_other
explVar_other[sig_other == 0] <- 0
explVar_other

explVar <- cbind(VARexpl_TRAIT1, explVar_other)
colnames(explVar) <- gsub("VARexpl_", "", colnames(explVar))

colnames(eff_sign_other) <- gsub("effSign_", "", colnames(eff_sign_other))
# now check if it is an unfavorable effect combination
explVar
deleted_effect <- list()
colnames(sig_other) <- gsub("sig_", "", colnames(sig_other))
for(i in colnames(explVar)[-1]){
	if(direction == "fav"){
		if(i %in% c("FF", "MF")){
		deleted_effect_i <- rownames(explVar)[which((eff_sign_other[,i] == 1) & (sig_other[ , i] == 1))]
		explVar[which(eff_sign_other[,i] == 1),i] <- 0
		} else {
		deleted_effect_i <- rownames(explVar)[which((eff_sign_other[,i] == -1) & (sig_other[ , i] == 1))]
		explVar[which(eff_sign_other[,i] == -1),i] <- 0
		}
	}
	if(direction == "unfav"){
		if(i %in% c("FF", "MF")){
		deleted_effect_i <- rownames(explVar)[which((eff_sign_other[,i] == -1) & (sig_other[ , i] == 1))]
		explVar[which(eff_sign_other[,i] == -1),i] <- 0
		} else {
		deleted_effect_i <- rownames(explVar)[which((eff_sign_other[,i] == 1) & (sig_other[ , i] == 1))]
		explVar[which(eff_sign_other[,i] == 1),i] <- 0
		}
	}
	deleted_effect[[i]] <- deleted_effect_i
}
explVar
colnames(explVar)[which(colnames(explVar) == "TRAIT1")] <- TRAIT1
colnames(explVar)[which(colnames(explVar) == "FF")] <- "FF"
colnames(explVar)[which(colnames(explVar) == "MF")] <- "MF"
colnames(explVar)[which(colnames(explVar) == "LO")] <- "LO"
explVar
explVar <- explVar[apply(explVar[, -1], 1, function(x) {any(x != 0)}), ]
explVar

write.table(explVar, paste(outfolder, "/varExpl_onlySig_", TRAIT1, "vsOthers_", direction, ".csv", sep = ""),
			sep = ";", dec = ".", quote = FALSE, row.names = TRUE, col.names = NA)
save(deleted_effect, file = paste(outfolder, "/deleted_FavEffect_", TRAIT1, "vsOthers_", direction, ".RData", sep = ""))

sig_TRAIT1[rownames(explVar), ]

all_unfav_withMultEff <- c(all_unfav_withMultEff, rownames(explVar))

# create barplot
cols <- colorRampPalette(c("white",rgb(0.43,0,0, alpha = 1)))(20)
Varexpl <- explVar[order(explVar[, TRAIT1], decreasing = TRUE), ]
Varexpl <- t(Varexpl)

								#if(TRAIT1 == "PH_V6"){
								#names.arg <- c("B", LETTERS_MM[1:(ncol(Varexpl)-1)])
								#} else {
								names.arg <- LETTERS_MM[1:ncol(Varexpl)]
								#}

						png(paste(outfolder, "/barplot_varExpl_", TRAIT1, "vsOthers_", direction, ".png", sep = ""), width = 300 + 2400 * (ncol(Varexpl) / 9), height = 2100, res = 300)
						par(mar = c(3.8, 3.9, 0, 0) + 0.1, mgp = c(2.5,0.75,0))
						barplot(height = Varexpl,
								beside = TRUE,
								col = c(cols[20], cols[16], cols[13], cols[10], cols[7], cols[4]),
								legend.text = TRUE,
								names.arg = names.arg,
								xlab = paste0("Unfavorable haplotypes with effect on ", TRAIT1, " and other traits"),
								ylab = "Prop. of genetic variance expl.",
								args.legend = list(x =  "top", cex = 1.5),
								cex.lab = 1.25,
								cex.axis = 1.25,
								cex.names = 1.25
								)
						dev.off()


}

#
# check deleted haplotypes due to "favorable" effect combination
#

direction <- "fav"
for(TRAIT1 in traits1){
	load(paste(outfolder, "/deleted_FavEffect_", TRAIT1, "vsOthers_", direction, ".RData", sep = ""))
	print(TRAIT1)
	print(deleted_effect)	
}

direction <- "unfav"
for(TRAIT1 in traits1){
	load(paste(outfolder, "/deleted_FavEffect_", TRAIT1, "vsOthers_", direction, ".RData", sep = ""))
	print(TRAIT1)
	print(deleted_effect)	
}

all_fav_withMultEff <- unique(all_fav_withMultEff)
sort(all_fav_withMultEff)
