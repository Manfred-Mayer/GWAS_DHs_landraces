###################################################################
###################################################################
####
#### construction of haplotypes
#### for the combined set of DHs and BLs
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


# load packages
library(synbreed)
library(zoo)

# arguments
nSNPs <- 10
steps <- 10


# functions
# function for recoding geno into "haplotype alleles"
recode_hapAlleles_f <- function(markers, geno){
	geno_wind <- geno[, markers]
	geno_wind_unique <- unique(geno_wind)
	haps <- apply(geno_wind, 1, paste0, collapse ="")
	return(haps)
	}

# function for writing a Markerlist containing the names of markers within sliding windows along each chromosome
FixedLengthMarkerlist_f <- function(map, nSNPs, steps) {
	Markerlist <- NULL
	for (Chr in unique(map$chr)){	
	SNPs_chr <- rownames(map)[map$chr == Chr]
	pos_chr <- map$pos[map$chr == Chr]
	names(pos_chr) <- SNPs_chr
	Markerlist <- rbind(Markerlist, rollapply(SNPs_chr, width = nSNPs, FUN = function(x){x}, by = steps))
	}
	return(Markerlist)
	}

# function for generating Infolist of the corresponding blocks
makeInfolist_f <- function(x, map) {
	map_temp <- map[which(rownames(map) %in% x),]
	Chr <- unique(map_temp$chr)
	if(length(Chr) != 1){
		print(x)
		print("block goes across chromosomes!")
		}
	pos <- map_temp$pos
	names(pos) <- rownames(map_temp)
	info <- c(Chr, pos[x[1]], pos[x[length(x)]], pos[x[length(x)]] - pos[x[1]] + 1, length(x))
	names(info) <- c("Chr", "Pos_Start_bp", "Pos_End_bp", "Size_bp", "n_Markers")
	return(info)
	}

# function for recoding vector of haplotype alleles (1,2,3,4,5,...) into binary format with one column per haplotype allele
# (if you compare it with the function used in 1.02, you notice, that here no MAF filter is used, i.e. every haplotype is counted)
binary_coding_f <- function(x) {
haps <- unique(x)
binary <- matrix(x, nrow = length(x), ncol = length(haps))
binary <- t(apply(binary, 1, helper_binary_coding_f, haps = haps))
if(nrow(binary) == 1){
binary <- t(binary)
}
colnames(binary) <- haps
return(binary)
}
helper_binary_coding_f <- function(y, haps){
	y[y != haps] <- 0
	y[y == haps] <- 2	
	return(y)
}


# load genotypic data
load("gpDH.BL.RData")

# do this for all 10 chromosomes separately
for(CHR in 1:10){
print(paste("CHR", CHR))

# filter map and geno for CHR
map <- gpDH.BL$map[which(gpDH.BL$map$chr == CHR), ]
str(map)
geno <- rbind(gpDH.BL$geno[ , rownames(map)])
str(geno)

# define windows
# write Markerlist according to the specified windows
	Marker_List <- FixedLengthMarkerlist_f(map = map, nSNPs = nSNPs, steps = steps)
	rownames(Marker_List) <- paste("wind", seq(1, nrow(Marker_List), 1), sep="")
	rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 5)] <- paste(substr(rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 5)], 1, 4), "00000", substr(rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 5)], 5, 5), sep = "")
	rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 6)] <- paste(substr(rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 6)], 1, 4), "0000", substr(rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 6)], 5, 6), sep = "")
	rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 7)] <- paste(substr(rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 7)], 1, 4), "000", substr(rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 7)], 5, 7), sep = "")
	rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 8)] <- paste(substr(rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 8)], 1, 4), "00", substr(rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 8)], 5, 8), sep = "")
	rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 9)] <- paste(substr(rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 9)], 1, 4), "0", substr(rownames(Marker_List)[which(nchar(rownames(Marker_List)) == 9)], 5, 9), sep = "")
	str(Marker_List)

# define Info_List for these windows (position, size, ...)
	Info_List <- t(apply(Marker_List, 1, makeInfolist_f, map = map))
	str(Info_List)

	#################################################################################################################
	# define haplotypes per window
	haplo <- apply(Marker_List, 1, recode_hapAlleles_f, geno = geno)
	colnames(haplo) <- rownames(Marker_List)
	rownames(haplo) <- rownames(geno)
	str(haplo)

	#################################################################################################################
	# code binary
	haplo_binary_list <- apply(haplo, 2, binary_coding_f)

		# convert to matrix and assign successive numbers to haplotypes within a window
		i_chr <- which(Info_List[,1] == CHR)
		haplo_binary_list_chr <- haplo_binary_list[i_chr]
		haplo_binary <- matrix(unlist(haplo_binary_list_chr), nrow = nrow(haplo))
		ncol_f <- function(x){
			if(length(dim(x)) == 0){
				if(is.null(x)){
				return(0)
				} else {
				return(1)
				}
			} else {
				return(ncol(x))
			}
		}
		nhap_per_block <- unlist(lapply(haplo_binary_list_chr, ncol_f))
		numbers_within <- NULL
		for(i in 1:length(nhap_per_block)){
		if(nhap_per_block[i] != 0){
		numbers_within <- c(numbers_within, seq(1, nhap_per_block[i], 1))
		}
		}
		colnames(haplo_binary) <- paste(rep(names(haplo_binary_list_chr), nhap_per_block), "_", numbers_within, sep = "")
		rownames(haplo_binary) <- rownames(haplo)
		rm(haplo_binary_list_chr, nhap_per_block, numbers_within)
		gc()

		all.equal(rownames(haplo_binary), rownames(haplo))
		str(haplo_binary)

		# adjust names of haplotypes according to "wind_CHROMOSOENUMBER_WINDOWNUMBER_HAPLOTYPENUMBERWITHINWINDOW"
		tempExpr <- unlist(gregexpr("_",colnames(haplo_binary)))
		hap_names_chr <- substr(colnames(haplo_binary), 1, tempExpr - 1)
		n_hap_per_wind <- table(hap_names_chr)
		InfoList_chr <- Info_List[hap_names_chr, ]

	if(CHR %in% 1:9){
		if(length(n_hap_per_wind) > 9){
		hap_names_chr_new <- paste("wind_0", CHR, "_0000", 1:9, sep = "")
			if(length(n_hap_per_wind) > 99){
			hap_names_chr_new <- c(hap_names_chr_new, paste("wind_0", CHR, "_000", 10:99, sep = ""))
				if(length(n_hap_per_wind) > 999){
				hap_names_chr_new <- c(hap_names_chr_new, paste("wind_0", CHR, "_00", 100:999, sep = ""))
					if(length(n_hap_per_wind) > 9999){
					hap_names_chr_new <- c(hap_names_chr_new, paste("wind_0", CHR, "_0", 1000:9999, sep = ""))
					hap_names_chr_new <- c(hap_names_chr_new, paste("wind_0", CHR, "_", 10000:length(n_hap_per_wind), sep = ""))
					} else {
					hap_names_chr_new <- c(hap_names_chr_new, paste("wind_0", CHR, "_0", 1000:length(n_hap_per_wind), sep = ""))
					}			
				} else {
				hap_names_chr_new <- c(hap_names_chr_new, paste("wind_0", CHR, "_00", 100:length(n_hap_per_wind), sep = ""))
				}		
			} else {
			hap_names_chr_new <- c(hap_names_chr_new, paste("wind_0", CHR, "_000", 10:length(n_hap_per_wind), sep = ""))
			}
		} else {
		hap_names_chr_new <- paste("wind_0", CHR, "_0000", 1:length(n_hap_per_wind), sep = "")	
		}
	} else{
		if(length(n_hap_per_wind) > 9){
		hap_names_chr_new <- paste("wind_", CHR, "_0000", 1:9, sep = "")
			if(length(n_hap_per_wind) > 99){
			hap_names_chr_new <- c(hap_names_chr_new, paste("wind_", CHR, "_000", 10:99, sep = ""))
				if(length(n_hap_per_wind) > 999){
				hap_names_chr_new <- c(hap_names_chr_new, paste("wind_", CHR, "_00", 100:999, sep = ""))
					if(length(n_hap_per_wind) > 9999){
					hap_names_chr_new <- c(hap_names_chr_new, paste("wind_", CHR, "_0", 1000:9999, sep = ""))
					hap_names_chr_new <- c(hap_names_chr_new, paste("wind_", CHR, "_", 10000:length(n_hap_per_wind), sep = ""))
					} else {
					hap_names_chr_new <- c(hap_names_chr_new, paste("wind_", CHR, "_0", 1000:length(n_hap_per_wind), sep = ""))
					}			
				} else {
				hap_names_chr_new <- c(hap_names_chr_new, paste("wind_", CHR, "_00", 100:length(n_hap_per_wind), sep = ""))
				}		
			} else {
			hap_names_chr_new <- c(hap_names_chr_new, paste("wind_", CHR, "_000", 10:length(n_hap_per_wind), sep = ""))
			}
		} else {
		hap_names_chr_new <- paste("wind_", CHR, "_0000", 1:length(n_hap_per_wind), sep = "")	
		}
	}

	numbers_within <- NULL
	for(i in 1:length(n_hap_per_wind)){
	if(n_hap_per_wind[i] != 0){
	numbers_within <- c(numbers_within, seq(1, n_hap_per_wind[i], 1))
	}
	}
	colnames(haplo_binary) <- paste(rep(hap_names_chr_new, n_hap_per_wind), "_", numbers_within, sep = "")
	rownames(InfoList_chr) <- colnames(haplo_binary)

	print(str(haplo_binary))
	print(str(InfoList_chr))
	InfoList <- InfoList_chr

	# save geno and InfoList
	save(InfoList, file = paste("InfoList_DH.BL_chr", CHR, "_m", nSNPs, "s", steps, ".RData", sep = ""))
	save(haplo_binary, file = paste("haplo_binary_DH.BL_chr", CHR, "_m", nSNPs, "s", steps, ".RData", sep = ""))
	
	# also save the Marker_List and haplo objects (updated with the respective window names)
	rownames(Marker_List) <- unique(substr(rownames(InfoList), 1, 13))
	save(Marker_List, file = paste("Marker_List_DH.BL_chr", CHR, "_m", nSNPs, "s", steps, ".RData", sep=""))
	colnames(haplo) <- rownames(Marker_List)
	save(haplo, file = paste("haplo_DH.BL_chr", CHR, "_m", nSNPs, "s", steps, ".RData", sep=""))

}

###
#############################################################################################
