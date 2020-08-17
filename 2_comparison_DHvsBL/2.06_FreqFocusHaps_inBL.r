###################################################################
###################################################################
####
#### assess frequencies of focus haplotypes (identified with GWAS in DHs)
#### in breeding lines
#### Plot frequencies, distinguishing favorable, unfavorable, random haplotypes
#### identify favorable haplotypes absent in BL 
#### identify unfavorable haplotypes with high frequencies in BL 
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
library(synbreed)

# arguments
nSNPs <- 10
steps <- 10
p_thresh <- 0.01
FDR <- "15"
minCount <- 3

# graphical parameters
cex <- 0.6
cex.axis <- 1.6
cex.lab <- 1.6

# functions
# function for calculating haplotype similarity (probability of two randomly chosen gametes to show identical haplotypes)
HapSim_f <- function(x, exCount = 0) {
	allele_counts <- table(x)
	allele_freqs <- allele_counts/sum(allele_counts)
	expHet <- 1 - sum(allele_freqs^2)
  HapSim <- 1 - expHet
  return(HapSim)
}

# function for calculating Hudson and Kaplan's minimum number of historical recombination events
HudsonKaplan_Rec_f <- function(geno, pos, n_min){
  fourGamete_pair <- function(x, n_min){
    test <- paste(x[1 : floor(length(x)/2)], x[(floor(length(x)/2) + 1) : length(x)], sep = "_")
    test <- table(test)
	test <- test[which(test >= n_min)]
	if(length(test) == 4){
      return(1)
    } else {
      return(0)
    }
  }
  
  # first, exclude all monomorphic markers
  mono <- which(apply(geno, 2, function(x) {sd(x) == 0}))
  if(length(mono > 0)){
	  geno <- geno[ , -mono]
	  pos <- pos[-mono]
  }
  if(length(pos) > 1){
	  col_pairs <- cbind(colnames(geno)[sequence(1:ncol(geno))],colnames(geno)[rep(1:ncol(geno),1:ncol(geno))])
	  col_pairs <- col_pairs[which(!(col_pairs[,1]==col_pairs[,2])),]
	  
	  if(length(col_pairs) > 2){
		  start_bp <- pos[col_pairs[,1]]
		  end_bp <- pos[col_pairs[,2]]
		  
		  test_mtx <- rbind(geno[, col_pairs[,1]], geno[ , col_pairs[,2]])
		  D <- apply(test_mtx, 2, fourGamete_pair, n_min = n_min)
		  D_info <- cbind(start_bp, end_bp, D)
		  D_info <- D_info[which(D_info[ , "D"] == 1), ]
		  
		  if(length(D_info) == 3){
			n_rec <- 1
		  }
		  if(length(D_info) < 3){
			n_rec <- 0
		  }
		  if(length(D_info) > 3){
			D_info <- D_info[order(D_info[,1], D_info[,2]), ]
			
			# first step of the algorithm
			ex <- NULL
			for(i in 1:nrow(D_info)){
			  temp1 <- D_info[i,1] <= D_info[-i,1]
			  temp2 <- D_info[i,2] >= D_info[-i,2]
			  if(any(temp1 & temp2)){
				ex <- c(ex, i)
			  }
			}
			if(is.null(ex) == FALSE){
			  D_info <- D_info[-ex, ]
			}
		  }
		  
		  if(length(D_info) == 3){
			n_rec <- 1
		  }
		  if(length(D_info) > 3){
			# second step of the algorithm
			i <- 0
			repeat{
			  i <- i + 1
			  ex <- NULL
			  temp1 <- D_info[(i+1):nrow(D_info),1] > D_info[i,1]
			  temp2 <- D_info[(i+1):nrow(D_info),1] < D_info[i,2]
			  if(any(temp1 & temp2)){
				ex <- c(ex, (which(temp1 & temp2) + i))
			  }
			  if(is.null(ex) == FALSE){
				D_info <- D_info[-ex, ]
				i <- 0
			  }
			  if(length(D_info) == 3){
				break
			  }
			  if(i == (nrow(D_info)-1)){
				break
			  }
			}
		  }

		  if(length(D_info) == 3){
			n_rec <- 1
		  } else {
			n_rec <- nrow(D_info)
		  }
		  
	  } else {
		  test_vec <- c(geno[, col_pairs[1]], geno[ , col_pairs[2]])
		  D <- fourGamete_pair(test_vec, n_min = n_min)
		  if(D == 1){
			n_rec <- 1
		  } else {
			n_rec <- 0
		  }
	  }

  } else {
    n_rec <- 0
  }

  return(n_rec)
}

# function for generating Infolist of the corresponding blocks
makeInfolist_genet_f <- function(x, genet_map) {
	genet_map_noNA <- genet_map[which(is.na(genet_map[,3]) == FALSE), ]
	x2 <- x[which(x[,2] >= genet_map_noNA[1,2]), ]
	x2 <- x2[which(x2[,2] <= genet_map_noNA[nrow(genet_map_noNA),2]), ]
	genet_InfoList <- NULL
	for(i in 1:nrow(x2)){	
		# start pos
		dif <- x2[i,2] - genet_map_noNA[,2]
		dif <- dif[order(abs(dif))][1:2]
		phys <- genet_map_noNA[sort(names(dif)),2]
		genet <- genet_map_noNA[sort(names(dif)),3]
		genet_start <- genet[1] + (((x2[i,2] - phys[1]) * diff(genet)) / diff(phys))
		# end pos
		dif <- x2[i,3] - genet_map_noNA[,2]
		dif <- dif[order(abs(dif))][1:2]
		phys <- genet_map_noNA[sort(names(dif)),2]
		genet <- genet_map_noNA[sort(names(dif)),3]
		genet_end <- genet[1] + (((x2[i,3] - phys[1]) * diff(genet)) / diff(phys))
		# size
		genet_size <- genet_end - genet_start
		genet_InfoList <- rbind(genet_InfoList, c(x2[i,1], genet_start, genet_end, genet_size, x2[i,5]))
		}
	colnames(genet_InfoList) <- c("Chr", "Pos_Start_cM", "Pos_End_cM", "Size_cM", "n_Markers")
	rownames(genet_InfoList) <- rownames(x2)
	return(list(genet_InfoList = genet_InfoList, phys_InfoList = x2))
	}

###################################################################

# load data of DHs
load(paste("geno_InfoList_DH_m", nSNPs, "s", steps, ".RData", sep = ""))
str(InfoList)
str(geno)

# load genetic map (genetic and physical position for each marker)
# according to PH207 x EP1 mapping population
# I:\Projekte\MAZE\Science\Genetic_Data_Analyses\GenMapPH207vsEP1_interpolated_600kArray_Positions
load("map_chr_phy_gen_PH207vsEP1_SCAMmpiInterpol_B73v4.RData")
str(map_chr_phy_gen)

# write out each start and end position of each haplotype as genetic as well as physical positions
Info_Lists <- list()
for(CHR in 1:10){
print(paste("CHR", CHR))
InfoList_chr <- InfoList[which(InfoList[,1] == CHR), ]
# only consider first haplotype per window (start and end is the same for each haplotype in a window anyways)
InfoList_chr <- InfoList_chr[which(substr(rownames(InfoList_chr), nchar(rownames(InfoList_chr))-1, nchar(rownames(InfoList_chr))) == "_1"), ]
genet_map_chr <- map_chr_phy_gen[which(map_chr_phy_gen[,1] == CHR), ]
Info_Lists[[CHR]] <- makeInfolist_genet_f(x = InfoList_chr, genet_map = genet_map_chr)
}
str(Info_Lists)
save(Info_Lists, file = "Info_Lists_genet_phys.RData")

###################################################################

dir.create("comp_BL_DHs")
dir.create("comp_BL_DHs/focusHaps")

traits <- c("TILL", "LO",
			"MF", "FF",
			"PH_final",
			"PH_V6", "PH_V4",
			"EV_V6", "EV_V4" 
			)

traits <- c("EV_V4", "EV_V6", "PH_V4", "PH_V6", "TILL", "LO")

for(TRAIT in traits){

infolder_haps <- paste("FavUnfav_Stability/", TRAIT, sep = "")
infolder_haps

dir.create(paste("comp_BL_DHs/focusHaps/", TRAIT, sep = ""))
outFolder <- paste("comp_BL_DHs/focusHaps/", TRAIT, sep = "")
outFolder

# load defined haplotypes
load(paste(infolder_haps, "/haps_all_", TRAIT, ".RData", sep = ""))
haps_Fav_info
haps_Unfav_info
haps_Inter_info

all_markers <- c(haps_Fav_info$qtl_i, haps_Unfav_info$qtl_i, haps_Inter_info$qtl_i)
all_markers
# if in biallelic case, haplotypes were added, we have to do it here also for geno and InfoList
if(any(unique(all_markers) %in% colnames(geno) == FALSE)){
temp_x <- unique(all_markers)
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

###############################################################################################################
###
###
### compare ferquencies
###
###

###
### search for haplotypes in breeding lines
###

# load check genotypic data
# SNP data of BL
load("gpBL.RData")
str(gpBL)

# SNP data of LR
load(paste("gpDH.RData", sep=""))
str(gpDH)

# generate elite line haplotype dataset for "favorable" haplotypes
geno_BL_fav <- NULL
nR_BL_fav <- NULL
HH_BL_fav <- NULL
cM_BL_fav <- NULL
kb_BL_fav <- NULL
for(hap_i in haps_Fav_info$qtl_i){
	chr_qtl <- InfoList[hap_i, "Chr"]
	start_qtl <- InfoList[hap_i, "Pos_Start_bp"]
	end_qtl <- InfoList[hap_i, "Pos_End_bp"]
	IDs_qtl <- rownames(geno)[which(geno[, hap_i] == 2)]

	map_i <- gpBL$map[which(gpBL$map$chr == chr_qtl), ]
	map_i <- map_i[which(map_i$pos >= start_qtl), ]
	map_i <- map_i[which(map_i$pos <= end_qtl), ]

	seq_qtl <- gpDH$geno[ , rownames(map_i)]
	seq_qtl <- apply(seq_qtl, 1, paste0, collapse = "")
	seq_qtl <- unique(seq_qtl[IDs_qtl])
	
	haps_BL <-  gpBL$geno[ , rownames(map_i)]
		
		# nR
		pos_temp <- map_i$pos
		names(pos_temp) <- rownames(map_i)
		nR_BL_fav <- c(nR_BL_fav, HudsonKaplan_Rec_f(geno = haps_BL, pos = pos_temp, n_min = 1))
		
	haps_BL <- apply(haps_BL, 1, paste0, collapse = "")
	
		# HH
		 HH_BL_fav <- c(HH_BL_fav, 1 - HapSim_f(haps_BL, exCount = 0))
	
		# cM
		genet_pos <- Info_Lists[[chr_qtl]][[1]]
		genet_pos <- genet_pos[which(substr(rownames(genet_pos), 1, 13) == substr(hap_i, 1, 13)),]
		cM_BL_fav <- c(cM_BL_fav, genet_pos["Size_cM"])
		
		# kb
		phys_pos <- InfoList
		phys_pos <- phys_pos[which(substr(rownames(phys_pos), 1, 13) == substr(hap_i, 1, 13)),]
		kb_BL_fav <- c(kb_BL_fav, phys_pos[1 , "Size_bp"] / 1000)

	haps_BL <- ifelse(haps_BL == seq_qtl, 2, 0)
	
	geno_BL_fav <- cbind(geno_BL_fav, haps_BL)
}
if(is.null(geno_BL_fav) == FALSE){
colnames(geno_BL_fav) <- haps_Fav_info$qtl_i
str(geno_BL_fav)
names(nR_BL_fav) <- haps_Fav_info$qtl_i
names(HH_BL_fav) <- haps_Fav_info$qtl_i
names(cM_BL_fav) <- haps_Fav_info$qtl_i
names(kb_BL_fav) <- haps_Fav_info$qtl_i
print(summary(nR_BL_fav))
print(summary(HH_BL_fav))
print(summary(cM_BL_fav))
print(summary(kb_BL_fav))
save(geno_BL_fav, nR_BL_fav, HH_BL_fav, cM_BL_fav, kb_BL_fav, file = paste(outFolder, "/info_BL_fav_", TRAIT, ".RData", sep = ""))
}

# generate elite line haplotype dataset for "unfavorable" haplotypes
geno_BL_unfav <- NULL
nR_BL_unfav <- NULL
HH_BL_unfav <- NULL
cM_BL_unfav <- NULL
kb_BL_unfav <- NULL
for(hap_i in haps_Unfav_info$qtl_i){
	chr_qtl <- InfoList[hap_i, "Chr"]
	start_qtl <- InfoList[hap_i, "Pos_Start_bp"]
	end_qtl <- InfoList[hap_i, "Pos_End_bp"]
	IDs_qtl <- rownames(geno)[which(geno[, hap_i] == 2)]

	map_i <- gpBL$map[which(gpBL$map$chr == chr_qtl), ]
	map_i <- map_i[which(map_i$pos >= start_qtl), ]
	map_i <- map_i[which(map_i$pos <= end_qtl), ]

	seq_qtl <- gpDH$geno[ , rownames(map_i)]
	seq_qtl <- apply(seq_qtl, 1, paste0, collapse = "")
	seq_qtl <- unique(seq_qtl[IDs_qtl])
	
	haps_BL <-  gpBL$geno[ , rownames(map_i)]
		
		# nR
		pos_temp <- map_i$pos
		names(pos_temp) <- rownames(map_i)
		nR_BL_unfav <- c(nR_BL_unfav, HudsonKaplan_Rec_f(geno = haps_BL, pos = pos_temp, n_min = 1))
		
	haps_BL <- apply(haps_BL, 1, paste0, collapse = "")
	
		# HH
		 HH_BL_unfav <- c(HH_BL_unfav, 1 - HapSim_f(haps_BL, exCount = 0))
	
		# cM
		genet_pos <- Info_Lists[[chr_qtl]][[1]]
		genet_pos <- genet_pos[which(substr(rownames(genet_pos), 1, 13) == substr(hap_i, 1, 13)),]
		cM_BL_unfav <- c(cM_BL_unfav, genet_pos["Size_cM"])
		
		# kb
		phys_pos <- InfoList
		phys_pos <- phys_pos[which(substr(rownames(phys_pos), 1, 13) == substr(hap_i, 1, 13)),]
		kb_BL_unfav <- c(kb_BL_unfav, phys_pos[1 , "Size_bp"] / 1000)

	haps_BL <- ifelse(haps_BL == seq_qtl, 2, 0)
	
	geno_BL_unfav <- cbind(geno_BL_unfav, haps_BL)
}
if(is.null(geno_BL_unfav) == FALSE){
colnames(geno_BL_unfav) <- haps_Unfav_info$qtl_i
str(geno_BL_unfav)
names(nR_BL_unfav) <- haps_Unfav_info$qtl_i
names(HH_BL_unfav) <- haps_Unfav_info$qtl_i
names(cM_BL_unfav) <- haps_Unfav_info$qtl_i
names(kb_BL_unfav) <- haps_Unfav_info$qtl_i
print(summary(nR_BL_unfav))
print(summary(HH_BL_unfav))
print(summary(cM_BL_unfav))
print(summary(kb_BL_unfav))
save(geno_BL_unfav, nR_BL_unfav, HH_BL_unfav, cM_BL_unfav, kb_BL_unfav, file = paste(outFolder, "/info_BL_unfav_", TRAIT, ".RData", sep = ""))
}


# generate elite line haplotype dataset for "neutral" haplotypes (in the same window as the significant haplotypes but without significant effect)
geno_BL_neut <- NULL
nR_BL_neut <- NULL
HH_BL_neut <- NULL
cM_BL_neut <- NULL
kb_BL_neut <- NULL
names_haps_neut <- NULL
for(hap_ii in c(haps_Fav_info$qtl_i, haps_Unfav_info$qtl_i)){
	# don't take the actual qtl, but all the others which are not significant
	haps_neut_ii <- colnames(geno)[which(substr(colnames(geno), 1, 13) == substr(hap_ii, 1, 13))]
	haps_neut_ii <- haps_neut_ii[which(haps_neut_ii %in% c(haps_Fav_info$qtl_i, haps_Unfav_info$qtl_i, haps_Inter_info$qtl_i) == FALSE)]
	
	for(hap_i in haps_neut_ii){
		names_haps_neut <- c(names_haps_neut, hap_i)
		chr_qtl <- InfoList[hap_i, "Chr"]
		start_qtl <- InfoList[hap_i, "Pos_Start_bp"]
		end_qtl <- InfoList[hap_i, "Pos_End_bp"]
		IDs_qtl <- rownames(geno)[which(geno[, hap_i] == 2)]

		map_i <- gpBL$map[which(gpBL$map$chr == chr_qtl), ]
		map_i <- map_i[which(map_i$pos >= start_qtl), ]
		map_i <- map_i[which(map_i$pos <= end_qtl), ]

		seq_qtl <- gpDH$geno[ , rownames(map_i)]
		seq_qtl <- apply(seq_qtl, 1, paste0, collapse = "")
		seq_qtl <- unique(seq_qtl[IDs_qtl])
		
		haps_BL <-  gpBL$geno[ , rownames(map_i)]
			
			# nR
			pos_temp <- map_i$pos
			names(pos_temp) <- rownames(map_i)
			nR_BL_neut <- c(nR_BL_neut, HudsonKaplan_Rec_f(geno = haps_BL, pos = pos_temp, n_min = 1))
			
		haps_BL <- apply(haps_BL, 1, paste0, collapse = "")
		
			# HH
			 HH_BL_neut <- c(HH_BL_neut, 1 - HapSim_f(haps_BL, exCount = 0))
		
			# cM
			genet_pos <- Info_Lists[[chr_qtl]][[1]]
			genet_pos <- genet_pos[which(substr(rownames(genet_pos), 1, 13) == substr(hap_i, 1, 13)),]
			cM_BL_neut <- c(cM_BL_neut, genet_pos["Size_cM"])
			
			# kb
			phys_pos <- InfoList
			phys_pos <- phys_pos[which(substr(rownames(phys_pos), 1, 13) == substr(hap_i, 1, 13)),]
			kb_BL_neut <- c(kb_BL_neut, phys_pos[1 , "Size_bp"] / 1000)

		haps_BL <- ifelse(haps_BL == seq_qtl, 2, 0)
		
		geno_BL_neut <- cbind(geno_BL_neut, haps_BL)
	}
}
if(is.null(geno_BL_neut) == FALSE){
colnames(geno_BL_neut) <- names_haps_neut
str(geno_BL_neut)
names(nR_BL_neut) <- names_haps_neut
names(HH_BL_neut) <- names_haps_neut
names(cM_BL_neut) <- names_haps_neut
names(kb_BL_neut) <- names_haps_neut
print(summary(nR_BL_neut))
print(summary(HH_BL_neut))
print(summary(cM_BL_neut))
print(summary(kb_BL_neut))
save(geno_BL_neut, nR_BL_neut, HH_BL_neut, cM_BL_neut, kb_BL_neut, file = paste(outFolder, "/info_BL_neut_", TRAIT, ".RData", sep = ""))
}



# generate elite line haplotype dataset for "random" haplotypes
set.seed(232)
rand_haps <- sample(colnames(geno), size = 500, replace = FALSE)
rand_haps <- sort(rand_haps)
rand_haps
geno_BL_rand <- NULL
nR_BL_rand <- NULL
HH_BL_rand <- NULL
cM_BL_rand <- NULL
kb_BL_rand <- NULL
for(hap_i in rand_haps){
	chr_qtl <- InfoList[hap_i, "Chr"]
	start_qtl <- InfoList[hap_i, "Pos_Start_bp"]
	end_qtl <- InfoList[hap_i, "Pos_End_bp"]
	IDs_qtl <- rownames(geno)[which(geno[, hap_i] == 2)]

	map_i <- gpBL$map[which(gpBL$map$chr == chr_qtl), ]
	map_i <- map_i[which(map_i$pos >= start_qtl), ]
	map_i <- map_i[which(map_i$pos <= end_qtl), ]

	seq_qtl <- gpDH$geno[ , rownames(map_i)]
	seq_qtl <- apply(seq_qtl, 1, paste0, collapse = "")
	seq_qtl <- unique(seq_qtl[IDs_qtl])
	
	haps_BL <-  gpBL$geno[ , rownames(map_i)]
		
		# nR
		pos_temp <- map_i$pos
		names(pos_temp) <- rownames(map_i)
		nR_BL_rand <- c(nR_BL_rand, HudsonKaplan_Rec_f(geno = haps_BL, pos = pos_temp, n_min = 1))
		
	haps_BL <- apply(haps_BL, 1, paste0, collapse = "")
	
		# HH
		 HH_BL_rand <- c(HH_BL_rand, 1 - HapSim_f(haps_BL, exCount = 0))
	
		# cM
		genet_pos <- Info_Lists[[chr_qtl]][[1]]
		genet_pos <- genet_pos[which(substr(rownames(genet_pos), 1, 13) == substr(hap_i, 1, 13)),]
		cM_BL_rand <- c(cM_BL_rand, genet_pos["Size_cM"])
		
		# kb
		phys_pos <- InfoList
		phys_pos <- phys_pos[which(substr(rownames(phys_pos), 1, 13) == substr(hap_i, 1, 13)),]
		if(length(phys_pos) > 5){
		kb_BL_rand <- c(kb_BL_rand, phys_pos[1 , "Size_bp"] / 1000)
		} else {
		kb_BL_rand <- c(kb_BL_rand, phys_pos["Size_bp"] / 1000)		
		}

	haps_BL <- ifelse(haps_BL == seq_qtl, 2, 0)
	
	geno_BL_rand <- cbind(geno_BL_rand, haps_BL)
}
colnames(geno_BL_rand) <- rand_haps
str(geno_BL_rand)
names(nR_BL_rand) <- rand_haps
names(HH_BL_rand) <- rand_haps
names(cM_BL_rand) <- rand_haps
names(kb_BL_rand) <- rand_haps
print(summary(nR_BL_rand))
print(summary(HH_BL_rand))
print(summary(cM_BL_rand))
print(summary(kb_BL_rand))
save(geno_BL_rand, nR_BL_rand, HH_BL_rand, cM_BL_rand, kb_BL_rand, file = paste(outFolder, "/info_BL_rand_", TRAIT, ".RData", sep = ""))


}


################################################
###
### summarize all early traits
###
################################################

name_set <- "earlyTraits"
dir.create(paste("comp_BL_DHs/focusHaps/", name_set, sep = ""))
outFolder <- paste("comp_BL_DHs/focusHaps/", name_set, sep = "")
outFolder

traits <- c("EV_V4", "EV_V6", "PH_V4", "PH_V6")
traits

geno_BL_fav_all <- NULL
geno_BL_unfav_all <- NULL
geno_BL_rand_all <- NULL

nR_BL_fav_all <- NULL
nR_BL_unfav_all <- NULL
nR_BL_rand_all <- NULL

HH_BL_fav_all <- NULL
HH_BL_unfav_all <- NULL
HH_BL_rand_all <- NULL

cM_BL_fav_all <- NULL
cM_BL_unfav_all <- NULL
cM_BL_rand_all <- NULL

kb_BL_fav_all <- NULL
kb_BL_unfav_all <- NULL
kb_BL_rand_all <- NULL

for(TRAIT in traits){

infolder <- paste("comp_BL_DHs/focusHaps/", TRAIT, sep = "")

load(paste(infolder, "/info_BL_fav_", TRAIT, ".RData", sep = ""))
load(paste(infolder, "/info_BL_unfav_", TRAIT, ".RData", sep = ""))
load(paste(infolder, "/info_BL_rand_", TRAIT, ".RData", sep = ""))

infolder_haps <- paste("FavUnfav_Stability/", TRAIT, sep = "")
load(paste(infolder_haps, "/haps_all_", TRAIT, ".RData", sep = ""))
all_markers <- c(haps_Fav_info$qtl_i, haps_Unfav_info$qtl_i, haps_Inter_info$qtl_i)
# if in biallelic case, haplotypes were added, we have to do it here also for geno and InfoList
if(any(unique(all_markers) %in% colnames(geno) == FALSE)){
temp_x <- unique(all_markers)
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

geno_BL_fav_all <- cbind(geno_BL_fav_all, geno_BL_fav)
nR_BL_fav_all <- c(nR_BL_fav_all, nR_BL_fav)
HH_BL_fav_all <- c(HH_BL_fav_all, HH_BL_fav)
cM_BL_fav_all <- c(cM_BL_fav_all, cM_BL_fav)
kb_BL_fav_all <- c(kb_BL_fav_all, kb_BL_fav)

geno_BL_unfav_all <- cbind(geno_BL_unfav_all, geno_BL_unfav)
nR_BL_unfav_all <- c(nR_BL_unfav_all, nR_BL_unfav)
HH_BL_unfav_all <- c(HH_BL_unfav_all, HH_BL_unfav)
cM_BL_unfav_all <- c(cM_BL_unfav_all, cM_BL_unfav)
kb_BL_unfav_all <- c(kb_BL_unfav_all, kb_BL_unfav)

}

# random are the same in every run
geno_BL_rand_all <- cbind(geno_BL_rand_all, geno_BL_rand)
nR_BL_rand_all <- c(nR_BL_rand_all, nR_BL_rand)
HH_BL_rand_all <- c(HH_BL_rand_all, HH_BL_rand)
cM_BL_rand_all <- c(cM_BL_rand_all, cM_BL_rand)
kb_BL_rand_all <- c(kb_BL_rand_all, kb_BL_rand)


summary(nR_BL_fav_all)
summary(nR_BL_unfav_all)
summary(nR_BL_rand_all)

summary(HH_BL_fav_all)
summary(HH_BL_unfav_all)
summary(HH_BL_rand_all)

summary(cM_BL_fav_all)
summary(cM_BL_unfav_all)
summary(cM_BL_rand_all)

summary(kb_BL_fav_all)
summary(kb_BL_unfav_all)
summary(kb_BL_rand_all)

haps_fav_all <- colnames(geno_BL_fav_all)
head(sort(table(haps_fav_all), decreasing = TRUE), n = 10)
length(haps_fav_all)
length(unique(haps_fav_all))

haps_unfav_all <- colnames(geno_BL_unfav_all)
head(sort(table(haps_unfav_all), decreasing = TRUE), n = 10)
length(haps_unfav_all)
length(unique(haps_unfav_all))

haps_rand_all <- colnames(geno_BL_rand_all)
head(sort(table(haps_rand_all), decreasing = TRUE), n = 10)
length(haps_rand_all)
length(unique(haps_rand_all))

# filter for non-overlapping haplotypes
geno_BL_fav_all <- geno_BL_fav_all[ , unique(colnames(geno_BL_fav_all))]
geno_BL_unfav_all <- geno_BL_unfav_all[ , unique(colnames(geno_BL_unfav_all))]
geno_BL_rand_all <- geno_BL_rand_all[ , unique(colnames(geno_BL_rand_all))]
str(geno_BL_fav_all)
str(geno_BL_unfav_all)
str(geno_BL_rand_all)

#
# calculate frequencies
#

str(geno_BL_fav_all)
str(geno_BL_unfav_all)
str(geno_BL_rand_all)

#
# calculate frequencies including
#

AF_f <- function(x){
        y <- sum(length(which(x==2))*2+length(which(x==1)))/(2*(length(x)-length(which(is.na(x)))))
        return(y)
}

### elite
freq_fav_elite_all <- apply(geno_BL_fav_all, 2, AF_f)
freq_unfav_elite_all <- apply(geno_BL_unfav_all, 2, AF_f)
freq_rand_elite_all <- apply(geno_BL_rand_all, 2, AF_f)
summary(freq_fav_elite_all)
summary(freq_unfav_elite_all)
summary(freq_rand_elite_all)

### LR
freq_fav_LR_all <- apply(geno[ ,colnames(geno_BL_fav_all)], 2, AF_f)
freq_unfav_LR_all <- apply(geno[ ,colnames(geno_BL_unfav_all)], 2, AF_f)
freq_rand_LR_all <- apply(geno[ ,colnames(geno_BL_rand_all)], 2, AF_f)
summary(freq_fav_LR_all)
summary(freq_unfav_LR_all)
summary(freq_rand_LR_all)


##
## permutation test for significant differences between fav/unfav/rand in elite panel
##

# permut est function (two.sided)
permut.test_f <- function(x, y, n = 10000, alternative = "two.sided"){
# remove NAs
x <- na.omit(x)
y <- na.omit(y)
# count entries per vector
nx <- length(x)
ny <- length(y)
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

set.seed(2121)

summary(freq_fav_elite_all)
summary(freq_unfav_elite_all)
summary(freq_rand_elite_all)

permTest_favRand <- permut.test_f(x = freq_fav_elite_all, y = freq_rand_elite_all)
permTest_favUnfav <- permut.test_f(x = freq_fav_elite_all, y = freq_unfav_elite_all)
permTest_unfavRand <- permut.test_f(x = freq_unfav_elite_all, y = freq_rand_elite_all)

str(permTest_favRand)
str(permTest_favUnfav)
str(permTest_unfavRand)

save(freq_fav_elite_all, freq_unfav_elite_all, freq_rand_elite_all, file = paste(outFolder, "/info_freqs_", name_set, "_nonDup.RData", sep = ""))







################################################
###
### plotting
###
################################################

##
## generate vectors of "independent" haplotypes (dist > 1Mb and/or r2<0.8)
##

#
# check overlapping regions between traits
#

traits <- c("PH_V4", "PH_V6", "EV_V4", "EV_V6")
traits

HAP_list <- list()

for (TRAIT in traits){
  
	infolder_haps <- paste("FavUnfav_Stability/", TRAIT, sep = "")
	load(paste(infolder_haps, "/haps_all_", TRAIT, ".RData", sep = ""))
  haps_fav_i <- haps_Fav_info$qtl_i
  haps_unfav_i <- haps_Unfav_info$qtl_i
  haps_inter_i <- haps_Inter_info$qtl_i
  haps_all <- c(haps_fav_i, haps_unfav_i, haps_inter_i)
  haps <- list()
  haps[["fav"]] <- haps_fav_i
  haps[["unfav"]] <- haps_unfav_i
  haps[["inter"]] <- haps_inter_i
  haps[["all"]] <- haps_all
  HAP_list[[TRAIT]] <- haps
  
}
str(HAP_list)

HAP_list_fav <- unlist(lapply(HAP_list, function(x) {x[[1]]}))
str(HAP_list_fav)
HAP_list_unfav <- unlist(lapply(HAP_list, function(x) {x[[2]]}))
str(HAP_list_unfav)

HAP_list_fav[which(HAP_list_fav %in% names(freq_fav_elite_all) == FALSE)]
HAP_list_unfav[which(HAP_list_unfav %in% names(freq_unfav_elite_all) == FALSE)]

names(freq_fav_elite_all)[which(names(freq_fav_elite_all) %in% HAP_list_fav == FALSE)]
names(freq_unfav_elite_all)[which(names(freq_unfav_elite_all) %in% HAP_list_unfav == FALSE)]

new_favorable_haps <- names(freq_fav_elite_all)[which(freq_fav_elite_all == 0)]
new_favorable_haps
for(i in names(HAP_list)){
  print(i)
  print(HAP_list[[i]]$fav[which(HAP_list[[i]]$fav %in% new_favorable_haps)])
}


common_unfavorable_haps <- names(freq_unfav_elite_all)[which(freq_unfav_elite_all > 0.25)]
common_unfavorable_haps
for(i in names(HAP_list)){
  print(i)
  print(HAP_list[[i]]$unfav[which(HAP_list[[i]]$unfav %in% common_unfavorable_haps)])
}


common_favorable_haps <- names(freq_fav_elite_all)[which(freq_fav_elite_all > 0.25)]
common_favorable_haps
for(i in names(HAP_list)){
  print(i)
  print(HAP_list[[i]]$fav[which(HAP_list[[i]]$fav %in% common_favorable_haps)])
}


#
# now load the according associated genomic regions
#

REG_list <- list()

for (TRAIT in traits){
  
  regs <- read.table(paste(TRAIT, "/QTLregs/finalQTL/finalQTLregs_", TRAIT, ".csv", sep = ""),
                     sep = ";",
                     dec = ".",
                     header = TRUE)
  rownames(regs) <- regs$qtl_i
  REG_list[[TRAIT]] <- regs
 
}
str(REG_list)

#
# calculate for traits overlapping associations, defined as:
# focus haplotypes within 1Mb and with r2 > 0.8
#
#

# filter for the 899 gphenotyped ones
load("pheno_perEnv_list_DHs.RData")
phenotyped <- pheno_perEnv_list[[1]]$Geno
ID_set <- intersect(rownames(geno), phenotyped)

# generate pairs
traits_2 <- combn(traits, 2)
traits_2

#
# for favorables
#
direction <- "fav"

# pairs
pairs_overlap_list_fav <- list()

for(trait_i in 1:ncol(traits_2)){
  TRAIT1 <- traits_2[1,trait_i]
  TRAIT2 <- traits_2[2,trait_i]
  print(paste(TRAIT1, "vs", TRAIT2))
  info_r2 <- NULL
  
  # get positions of focus haplotypes
  regs1 <- InfoList[HAP_list[[TRAIT1]][[direction]], 1:3]
  regs2 <- InfoList[HAP_list[[TRAIT2]][[direction]], 1:3]
  
  # extend for 0.5 Mb up and downstream
  regs1_ext1Mb <- regs1
  regs1_ext1Mb[, 2] <- regs1_ext1Mb[,2] - 2500000
  regs1_ext1Mb[, 3] <- regs1_ext1Mb[,3] + 2500000
  regs2_ext1Mb <- regs2
  regs2_ext1Mb[, 2] <- regs2_ext1Mb[,2] - 2500000
  regs2_ext1Mb[, 3] <- regs2_ext1Mb[,3] + 2500000
  
  # check for pairs within 1Mb distance and in case, calculate r2
  for(regs1_i in 1:nrow(regs1_ext1Mb)){
    regs2_chr <- rownames(regs2_ext1Mb)[which(regs2_ext1Mb[,1] == regs1_ext1Mb[regs1_i,1])]
    if(length(regs2_chr) > 0){
      haps2 <- regs2_chr[which((regs2_ext1Mb[regs2_chr, 2] < regs1_ext1Mb[regs1_i,3]) & (regs2_ext1Mb[regs2_chr, 3] > regs1_ext1Mb[regs1_i,2]))]
      if(length(haps2) > 0){
        haps1 <- rownames(regs1_ext1Mb)[regs1_i]
        r2 <- (cor(geno[ID_set, haps1], geno[ID_set, haps2]))^2
        info_r2_i <- data.frame(TRAIT1 = rep(haps1, length(haps2)),
                                TRAIT2 = haps2,
                                r2 = t(r2),
                                chr = regs1_ext1Mb[regs1_i,1],
                                start = min(c(regs1[haps1,2], regs2[haps2, 2])),
                                end = max(c(regs1[haps1,3], regs2[haps2, 3])),
                                stringsAsFactors = FALSE)
        info_r2 <- rbind(info_r2, info_r2_i)
      }
    }
  }
  if(is.null(info_r2) == FALSE){
    colnames(info_r2) <- c(TRAIT1, TRAIT2, "r2", "chr", "start_bp", "end_bp")
    rownames(info_r2) <- 1:nrow(info_r2)
    print(info_r2[which(info_r2$r2 > 0.4),])
  }
  pairs_overlap_list_fav[[paste(TRAIT1, TRAIT2, sep = ".")]] <- info_r2
}



#
# for unfavorables
#

direction <- "unfav"

# pairs
pairs_overlap_list_unfav <- list()

for(trait_i in 1:ncol(traits_2)){
  TRAIT1 <- traits_2[1,trait_i]
  TRAIT2 <- traits_2[2,trait_i]
  print(paste(TRAIT1, "vs", TRAIT2))
  info_r2 <- NULL
  
  # get positions of focus haplotypes
  regs1 <- InfoList[HAP_list[[TRAIT1]][[direction]], 1:3]
  regs2 <- InfoList[HAP_list[[TRAIT2]][[direction]], 1:3]
  
  # extend for 0.5 Mb up and downstream
  regs1_ext1Mb <- regs1
  regs1_ext1Mb[, 2] <- regs1_ext1Mb[,2] - 2500000
  regs1_ext1Mb[, 3] <- regs1_ext1Mb[,3] + 2500000
  regs2_ext1Mb <- regs2
  regs2_ext1Mb[, 2] <- regs2_ext1Mb[,2] - 2500000
  regs2_ext1Mb[, 3] <- regs2_ext1Mb[,3] + 2500000
  
  # check for pairs within 1Mb distance and in case, calculate r2
  for(regs1_i in 1:nrow(regs1_ext1Mb)){
    regs2_chr <- rownames(regs2_ext1Mb)[which(regs2_ext1Mb[,1] == regs1_ext1Mb[regs1_i,1])]
    if(length(regs2_chr) > 0){
      haps2 <- regs2_chr[which((regs2_ext1Mb[regs2_chr, 2] < regs1_ext1Mb[regs1_i,3]) & (regs2_ext1Mb[regs2_chr, 3] > regs1_ext1Mb[regs1_i,2]))]
      if(length(haps2) > 0){
        haps1 <- rownames(regs1_ext1Mb)[regs1_i]
        r2 <- (cor(geno[ID_set, haps1], geno[ID_set, haps2]))^2
        info_r2_i <- data.frame(TRAIT1 = rep(haps1, length(haps2)),
                                TRAIT2 = haps2,
                                r2 = t(r2),
                                chr = regs1_ext1Mb[regs1_i,1],
                                start = min(c(regs1[haps1,2], regs2[haps2, 2])),
                                end = max(c(regs1[haps1,3], regs2[haps2, 3])),
                                stringsAsFactors = FALSE)
        info_r2 <- rbind(info_r2, info_r2_i)
      }
    }
  }
  if(is.null(info_r2) == FALSE){
    colnames(info_r2) <- c(TRAIT1, TRAIT2, "r2", "chr", "start_bp", "end_bp")
    rownames(info_r2) <- 1:nrow(info_r2)
    print(info_r2[which(info_r2$r2 > 0.4),])
  }
  pairs_overlap_list_unfav[[paste(TRAIT1, TRAIT2, sep = ".")]] <- info_r2
}


#
# for changing sign
#
direction <- "inter"

# pairs
pairs_overlap_list_inter <- list()

for(trait_i in 1:ncol(traits_2)){
  TRAIT1 <- traits_2[1,trait_i]
  TRAIT2 <- traits_2[2,trait_i]
  print(paste(TRAIT1, "vs", TRAIT2))
  info_r2 <- NULL
  
  # get positions of focus haplotypes
  regs1 <- InfoList[HAP_list[[TRAIT1]][[direction]], 1:3]
  regs2 <- InfoList[HAP_list[[TRAIT2]][[direction]], 1:3]
  
  # extend for 0.5 Mb up and downstream
  regs1_ext1Mb <- regs1
  regs1_ext1Mb[, 2] <- regs1_ext1Mb[,2] - 2500000
  regs1_ext1Mb[, 3] <- regs1_ext1Mb[,3] + 2500000
  regs2_ext1Mb <- regs2
  regs2_ext1Mb[, 2] <- regs2_ext1Mb[,2] - 2500000
  regs2_ext1Mb[, 3] <- regs2_ext1Mb[,3] + 2500000
  
  # check for pairs within 1Mb distance and in case, calculate r2
  for(regs1_i in 1:nrow(regs1_ext1Mb)){
    regs2_chr <- rownames(regs2_ext1Mb)[which(regs2_ext1Mb[,1] == regs1_ext1Mb[regs1_i,1])]
    if(length(regs2_chr) > 0){
      haps2 <- regs2_chr[which((regs2_ext1Mb[regs2_chr, 2] < regs1_ext1Mb[regs1_i,3]) & (regs2_ext1Mb[regs2_chr, 3] > regs1_ext1Mb[regs1_i,2]))]
      if(length(haps2) > 0){
        haps1 <- rownames(regs1_ext1Mb)[regs1_i]
        r2 <- (cor(geno[ID_set, haps1], geno[ID_set, haps2]))^2
        info_r2_i <- data.frame(TRAIT1 = rep(haps1, length(haps2)),
                                TRAIT2 = haps2,
                                r2 = t(r2),
                                chr = regs1_ext1Mb[regs1_i,1],
                                start = min(c(regs1[haps1,2], regs2[haps2, 2])),
                                end = max(c(regs1[haps1,3], regs2[haps2, 3])),
                                stringsAsFactors = FALSE)
        info_r2 <- rbind(info_r2, info_r2_i)
      }
    }
  }
  if(is.null(info_r2) == FALSE){
    colnames(info_r2) <- c(TRAIT1, TRAIT2, "r2", "chr", "start_bp", "end_bp")
    rownames(info_r2) <- 1:nrow(info_r2)
    print(info_r2[which(info_r2$r2 > 0.4),])
  }
  pairs_overlap_list_inter[[paste(TRAIT1, TRAIT2, sep = ".")]] <- info_r2
}



#
# for all
#
direction <- "all"

# pairs
pairs_overlap_list_all <- list()

for(trait_i in 1:ncol(traits_2)){
  TRAIT1 <- traits_2[1,trait_i]
  TRAIT2 <- traits_2[2,trait_i]
  print(paste(TRAIT1, "vs", TRAIT2))
  info_r2 <- NULL
  
  # get positions of focus haplotypes
  regs1 <- InfoList[HAP_list[[TRAIT1]][[direction]], 1:3]
  regs2 <- InfoList[HAP_list[[TRAIT2]][[direction]], 1:3]
  
  # extend for 0.5 Mb up and downstream
  regs1_ext1Mb <- regs1
  regs1_ext1Mb[, 2] <- regs1_ext1Mb[,2] - 2500000
  regs1_ext1Mb[, 3] <- regs1_ext1Mb[,3] + 2500000
  regs2_ext1Mb <- regs2
  regs2_ext1Mb[, 2] <- regs2_ext1Mb[,2] - 2500000
  regs2_ext1Mb[, 3] <- regs2_ext1Mb[,3] + 2500000
  
  # check for pairs within 1Mb distance and in case, calculate r2
  for(regs1_i in 1:nrow(regs1_ext1Mb)){
    regs2_chr <- rownames(regs2_ext1Mb)[which(regs2_ext1Mb[,1] == regs1_ext1Mb[regs1_i,1])]
    if(length(regs2_chr) > 0){
      haps2 <- regs2_chr[which((regs2_ext1Mb[regs2_chr, 2] < regs1_ext1Mb[regs1_i,3]) & (regs2_ext1Mb[regs2_chr, 3] > regs1_ext1Mb[regs1_i,2]))]
      if(length(haps2) > 0){
        haps1 <- rownames(regs1_ext1Mb)[regs1_i]
        r2 <- (cor(geno[ID_set, haps1], geno[ID_set, haps2]))^2
        info_r2_i <- data.frame(TRAIT1 = rep(haps1, length(haps2)),
                                TRAIT2 = haps2,
                                r2 = t(r2),
                                chr = regs1_ext1Mb[regs1_i,1],
                                start = min(c(regs1[haps1,2], regs2[haps2, 2])),
                                end = max(c(regs1[haps1,3], regs2[haps2, 3])),
                                stringsAsFactors = FALSE)
        info_r2 <- rbind(info_r2, info_r2_i)
      }
    }
  }
  if(is.null(info_r2) == FALSE){
    colnames(info_r2) <- c(TRAIT1, TRAIT2, "r2", "chr", "start_bp", "end_bp")
    rownames(info_r2) <- 1:nrow(info_r2)
    print(info_r2[which(info_r2$r2 > 0.4),])
  }
  pairs_overlap_list_all[[paste(TRAIT1, TRAIT2, sep = ".")]] <- info_r2
}


#
# for all detected pairs, delete randomly one haplotype
# how many "unrelated" haplotypes are left
#

fav_haps_all <- names(freq_fav_elite_all)
unfav_haps_all <- names(freq_unfav_elite_all)
str(fav_haps_all)
str(unfav_haps_all)

set.seed(232)
r2_thresh_i <- "0.8"

  haps_pairs <- NULL
  r2_thresh <- as.numeric(r2_thresh_i)
  ex <- NULL
    for(i in 1:ncol(traits_2)){
      haps_i <- pairs_overlap_list_fav[[i]]
      haps_i <- haps_i[which(haps_i[,1] != haps_i[,2]), ]
      haps_i <- cbind(haps_i[which(haps_i$r2 >= r2_thresh), 1], haps_i[which(haps_i$r2 >= r2_thresh), 2])
      if(nrow(haps_i) > 0){
        for(ii in 1:nrow(haps_i)){
          if(any(haps_i[ii,] %in% ex) == FALSE){
            ex <- c(ex, sample(haps_i[ii,], size = 1))
            haps_pairs <- rbind(haps_pairs, haps_i[ii,])
          }
        }
      }
    }
    freq_fav_elite_all[haps_pairs]
    ex
    fav_haps_indep <- fav_haps_all[which(fav_haps_all %in% ex == FALSE)]
    

    haps_pairs <- NULL
    r2_thresh <- as.numeric(r2_thresh_i)
    ex <- NULL
    for(i in 1:ncol(traits_2)){
      haps_i <- pairs_overlap_list_unfav[[i]]
      haps_i <- haps_i[which(haps_i[,1] != haps_i[,2]), ]
      haps_i <- cbind(haps_i[which(haps_i$r2 >= r2_thresh), 1], haps_i[which(haps_i$r2 >= r2_thresh), 2])
      if(nrow(haps_i) > 0){
        for(ii in 1:nrow(haps_i)){
          if(any(haps_i[ii,] %in% ex) == FALSE){
            ex <- c(ex, sample(haps_i[ii,], size = 1))
            haps_pairs <- rbind(haps_pairs, haps_i[ii,])
          }
        }
      }
    }
    freq_unfav_elite_all[haps_pairs]
    ex
    unfav_haps_indep <- unfav_haps_all[which(unfav_haps_all %in% ex == FALSE)]
    
str(fav_haps_indep)
str(unfav_haps_indep)

summary(freq_fav_elite_all)
summary(freq_unfav_elite_all)

summary(freq_fav_elite_all[fav_haps_indep])
summary(freq_unfav_elite_all[unfav_haps_indep])
summary(freq_rand_elite_all)

length(which(freq_unfav_elite_all[unfav_haps_indep] > quantile(freq_rand_elite_all, 0.75)))
length(which(freq_unfav_elite_all[unfav_haps_indep] > quantile(freq_rand_elite_all, 0.75))) / length(freq_unfav_elite_all[unfav_haps_indep])

length(which(freq_fav_elite_all[fav_haps_indep] == 0))
length(which(freq_fav_elite_all[fav_haps_indep] == 0)) / length(freq_fav_elite_all[fav_haps_indep])


##
## test difference between distributions
##
# Mann-Whitney
wilcox.test(x = freq_fav_elite_all[fav_haps_indep], y = freq_rand_elite_all, alternative = "two.sided")
wilcox.test(x = freq_unfav_elite_all[unfav_haps_indep], y = freq_rand_elite_all, alternative = "two.sided")
wilcox.test(x = freq_fav_elite_all[fav_haps_indep], y = freq_unfav_elite_all[unfav_haps_indep], alternative = "two.sided")

##
## draw frequency plot
##

name_set <- "earlyTraits"

freq_fav_elite_all_plot <- freq_fav_elite_all[fav_haps_indep]
freq_unfav_elite_all_plot <- freq_unfav_elite_all[unfav_haps_indep]

pos_vec <- freq_fav_elite_all_plot
neg_vec <- freq_unfav_elite_all_plot
rand_vec <- freq_rand_elite_all

str(rand_vec)
str(pos_vec)
str(neg_vec)

xlim1 <- c(0,1)
xlim2 <- c(-0.5,1)

d_rand <- density(rand_vec)
str(d_rand)
d_rand$y <- d_rand$y[which(d_rand$x > 0)]
d_rand$x <- d_rand$x[which(d_rand$x > 0)]
d_rand$y <- d_rand$y[which(d_rand$x < 1)]
d_rand$x <- d_rand$x[which(d_rand$x < 1)]
d_rand$x <- c(0, d_rand$x, 1)
d_rand$y <- c(0, d_rand$y, 0)

d_pos <- density(pos_vec)
str(d_pos)
d_pos$y <- d_pos$y[which(d_pos$x > 0)]
d_pos$x <- d_pos$x[which(d_pos$x > 0)]
d_pos$x <- c(0, d_pos$x)
d_pos$y <- c(0, d_pos$y)

d_neg <- density(neg_vec)
str(d_neg)
d_neg$y <- d_neg$y[which(d_neg$x > 0)]
d_neg$x <- d_neg$x[which(d_neg$x > 0)]
d_neg$x <- c(0, d_neg$x)
d_neg$y <- c(0, d_neg$y)

png(paste(outFolder, "/density_RandFavUnfav_", name_set, "_indep.png", sep =""), width = 2100, height = 1800, res = 300)
par(mar = c(4, 4, 0, 0) + 0.1, mgp = c(2.9, 1, 0))

plot(d_rand, col = rgb(0.2, 0.2, 0.2, alpha = 1), xlim = xlim1, lwd = 2, ylim = c(0, max(c(d_rand$y, d_pos$y, d_neg$y), na.rm = TRUE)),
     main = "", xlab = paste("Haplotype frequency in breeding lines", sep =""), bty = "n",
     cex.axis = 1.2, cex.lab = 1.5)
polygon(d_rand, col = rgb(0.2, 0.2, 0.2, alpha = 0.2), border=NA)

points(d_neg, col = "red", xlim = xlim2, lwd = 2, type = "l")
polygon(d_neg, col = rgb(1, 0, 0, alpha = 0.15), border=NA)

points(d_pos, col = "blue", xlim = xlim2, lwd = 2, type = "l")
polygon(d_pos, col = rgb(0, 0, 1, alpha = 0.15), border=NA)

legend("center", pch = 22, col = c("blue", "red", rgb(0.2, 0.2, 0.2, alpha = 1)),
       pt.bg = c(rgb(0, 0, 1, alpha = 0.15), rgb(1, 0, 0, alpha = 0.15), rgb(0.2, 0.2, 0.2, alpha = 0.2)),
       pt.lwd = 2, pt.cex = 2,
       text.col = c("blue", "red", rgb(0.2, 0.2, 0.2, alpha = 1)), cex = 1.2,
       legend = c("Favorable", "Unfavorable", "Random"),
       bty = "n")

dev.off()

###
###################################################################
