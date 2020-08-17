###################################################################
###################################################################
####
#### calculate minimum number of historical recombination events according to Hudson and Kaplan 1985
#### for DHs of KE, LL, PE separately as well as overall
#### for BLs
####
#### use same windows as for haplotypes
####
#### Manfred Mayer (Technical University of Munich, Plant Breeding)
#### manfred.mayer@tum.de
####
#### date: 30.06.2020
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

# load packages
library(synbreed)

# function for calculating Hudson and Kaplan's minimum number of historical recombination events
# geno = genotype matrix (rows = individuals, columns = markers)
# pos = position of markers (in bp)
# n_min = number of times a allele combination has to be present to be counted
HudsonKaplan_Rec_f <- function(geno, pos, n_min = 1){
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

###################################################################
###

# output folder
dir.create(paste("diversityParams", sep = ""))
dir.create(paste("diversityParams/nR", sep = ""))
outfolder <- paste("diversityParams/nR", sep = "")
outfolder

# load genotypic data of DH and BL lines
load("gpDH.BL.RData")
str(gpDH.BL)

# load haplotype window information
load(paste("geno_InfoList_DH.BL_m", nSNPs, "s", steps, ".RData", sep = ""))
str(geno)
rm(geno)
str(InfoList)
# do calculations for every window (that means in the haplotype marker set use only those with "_1" at the end)
InfoList <- InfoList[which(substr(rownames(InfoList), 14, 99) == "_1"), ]
rownames(InfoList) <- substr(rownames(InfoList), 1, 13)
str(InfoList)

info_nR <- matrix(NA, nrow = 5, ncol = nrow(InfoList))
rownames(info_nR) <- c("KE", "LL", "PE", "DHall", "BL")
colnames(info_nR) <- rownames(InfoList)
str(info_nR)

for(hap_i in rownames(InfoList)){
print(paste0("hap ", which(rownames(InfoList) == hap_i), " of ", nrow(InfoList)))

	chr_qtl <- InfoList[hap_i, "Chr"]
	start_qtl <- InfoList[hap_i, "Pos_Start_bp"]
	end_qtl <- InfoList[hap_i, "Pos_End_bp"]

	map_i <- gpDH.BL$map[which(gpDH.BL$map$chr == chr_qtl), ]
	map_i <- map_i[which(map_i$pos >= start_qtl), ]
	map_i <- map_i[which(map_i$pos <= end_qtl), ]
	
	geno_sub <-  gpDH.BL$geno[ , rownames(map_i)]
	pos_sub <- map_i$pos
	names(pos_sub) <- rownames(map_i)

	# calculate nR for each landrace separately
	for(LR in c("KE", "LL", "PE")){
	geno <- geno_sub[which(substr(rownames(geno_sub), 1, 5) == paste0("DH_", LR)), ]
	info_nR[LR, hap_i] <- HudsonKaplan_Rec_f(geno = geno, pos = pos_sub)
	rm(geno)
	}

	# calculate nR for all landraces together
	geno <- geno_sub[which(substr(rownames(geno_sub), 1, 3) == paste0("DH_")), ]
	info_nR["DHall", hap_i]  <- HudsonKaplan_Rec_f(geno = geno, pos = pos_sub)
	rm(geno)

	# calculate nR for BLs
	geno <- geno_sub[which(substr(rownames(geno_sub), 1, 3) == paste0("EL_")), ]
	info_nR["BL", hap_i]  <- HudsonKaplan_Rec_f(geno = geno, pos = pos_sub)
	rm(geno)

}

save(info_nR, file = paste(outfolder, "/info_nR.RData", sep = ""))

# average over all windows
info_nR_avg <- matrix(NA, nrow = 5, ncol = 1)
rownames(info_nR_avg) <- c("KE", "LL", "PE", "DHall", "BL")
colnames(info_nR_avg) <- "nR"
info_nR_avg

for(group_i in rownames(info_nR)){
	info_nR_avg[group_i, "nR"] <- mean(info_nR[group_i, ])
}
info_nR_avg

write.table(info_nR_avg, file = paste0(outfolder, "/info_nR_avg.csv"), sep = ";", dec = ".", quote = FALSE, col.names = FALSE, row.names = TRUE)

apply(info_nR, 1, summary)

###
###################################################################
