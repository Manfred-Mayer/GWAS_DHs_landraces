###################################################################
###################################################################
####
#### plot the trait associated genomic regions as defined in script 1.08
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

#
traits <- c("TILL", "LO",
			"MF", "FF",
			"PH_final",
			"PH_V6", "PH_V4",
			"EV_V6", "EV_V4" 
			)
white_space_size <- 25000000

# define centromere positions (B73v4)
Start_bp <- c(136.77, 95.51, 85.78, 109.07, 104.54, 52.3, 56.38, 50.53, 53.75, 57.36, 51.39)
End_bp <- c(137.12, 97.49, 86.93, 110.5, 106.82, 53.11, 56.68, 52.07, 55.39, 57.76, 52.78)
centromers <- cbind(Start_bp, End_bp)
rownames(centromers) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr9_2", "chr10")
centromers <- centromers[-grep("_", rownames(centromers)), ]
centromers_mean <- apply(centromers, 1, mean) * 1000000
chr_size <- c(307.041717, 244.442276, 235.667834, 246.994605, 223.902240, 174.033170, 182.381542, 181.122637, 159.769782, 150.982314)
chr_size_temp <- c(0, chr_size)*1000000
chr_ends <- 0
chr_add <- NULL
white_space <- list()
white_space[[1]] <- c(-999999999999, 0)
chr_size_white_space <- c(0, chr_size)*1000000
for(i in 2:11){
white_space[[i]] <- c(sum(chr_size_temp[1:i]), sum(chr_size_temp[1:i]) + white_space_size)
chr_size_temp[i] <- chr_size_temp[i] + white_space_size
}
white_space[[11]][2] <- 999999999999
for(CHR in 1:10){
chr_add <- c(chr_add, sum(chr_size_temp[1:CHR]))
chr_ends <- c(chr_ends, sum(chr_size_temp[1:CHR]) - white_space_size, sum(chr_size_temp[1:CHR]))
}
chr_ends <- c(chr_ends, sum(chr_size_temp[1:11]) - white_space_size)
names(chr_add) <- 1:10
centromers_mean <- centromers_mean + chr_add

	png(paste("GWASregs_Summary_FDR", FDR, "_p", p_thresh, ".png", sep =""), width = 6000, height = 3000, res = 300)
	par(mar = c(5, 7, 4, 2) + 0.1)
	
	plot(x = -99999999, y = -99999999, xlim = c(50000000, sum(chr_size_temp)), ylim = c(0.3, (length(traits))+0.5), main = "", ylab = "", xlab = "Chromosome", xaxt = "n", yaxt = "n", cex.axis = cex.axis, cex.lab = cex.lab, cex = cex)
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(0.95,0.95,0.95, alpha = 1))
	axis(side = 2, at = 1 : length(traits), labels = traits, cex.axis = cex.axis, cex.lab = cex.lab, cex = cex, las = 1)
	abline(h = seq(0.5, length(traits)+0.5, 1), lwd = 0.5)
	abline(v = chr_ends)
			
		for(TRAIT in traits){
			start_pos_all <- NULL
			end_pos_all <- NULL
				# input folder
					infolder <- paste(TRAIT, "/QTLregs/finalQTL", sep = "")
				if(file.exists(paste(infolder, "/finalQTLregs_", TRAIT, ".csv", sep = ""))){
				regs_temp <- read.table(paste(infolder, "/finalQTLregs_", TRAIT, ".csv", sep = ""), sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)
				start_pos_temp <- regs_temp$start_qtl
				start_pos_temp <- start_pos_temp + chr_add[regs_temp$chr_qtl]
				end_pos_temp <- regs_temp$end_qtl
				end_pos_temp <- end_pos_temp + chr_add[regs_temp$chr_qtl]
				start_pos_all <- c(start_pos_all, start_pos_temp)
				end_pos_all <- c(end_pos_all, end_pos_temp)
				}
			if(length(start_pos_all) > 0){
			y_pos_lower <- rep(which(traits == TRAIT) - 0.45, length(start_pos_all))
			y_pos_upper <- rep(which(traits == TRAIT) + 0.45, length(start_pos_all))
			for(i in 1:length(start_pos_all)){
				polygon(x = c(start_pos_all[i], start_pos_all[i], end_pos_all[i], end_pos_all[i]), y = c(y_pos_lower[i], y_pos_upper[i], y_pos_upper[i], y_pos_lower[i]), col = "black", border = "black", lwd = 2)
			}			
			}
			}

	polygon(x = c(- 99999999, - 99999999, 9999999999999, 9999999999999), y = c(-100, 0.5, 0.5, -100), col = "white", border = NA)
	polygon(x = c(- 99999999, - 99999999, 9999999999999, 9999999999999), y = c(100, length(traits)+0.5, length(traits)+0.5, 100), col = "white", border = NA)
	points(x = centromers_mean, y = rep(0.33, length(centromers_mean)), pch = 17, col = rgb(0, 0, 0, alpha = 0.5), cex = 2)
	abline(h = 0.1, lwd = 1.2)
	for(i in 1:length(white_space)){
		polygon(x = c(white_space[[i]][1], white_space[[i]][1], white_space[[i]][2], white_space[[i]][2]), y = c(-100, 100, 100, -100), col = "white", border = NA)	
	}	
	
	dev.off()
