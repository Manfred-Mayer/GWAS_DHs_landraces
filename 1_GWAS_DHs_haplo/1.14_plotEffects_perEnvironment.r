###################################################################
###################################################################
####
#### plot the for each trait the environment-specific effects
#### for the final set of focus haplotypes
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
TRAIT <- "PH_V6"
hap1 <- NULL
hap2 <- NULL
hap3 <- NULL
p_thresh <- 0.01
FDR <- "15"

# graphical parameters
cex <- 0.6
cex.axis <- 1.6
cex.lab <- 1.6

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
TRAIT
hap1
hap2
hap3

# input folder
infolder <- paste(TRAIT, "/QTLregs/finalQTL", sep = "")
infolder

dir.create(paste("Stability_GWASplots", sep = ""))
dir.create(paste("Stability_GWASplots/", TRAIT, sep = ""))
outfolder <- paste("Stability_GWASplots/", TRAIT, sep = "")
outfolder

Envs <- read.table(paste0(TRAIT, "/", list.files(path = TRAIT, pattern = ".phenotypes_order.txt")), stringsAsFactors = FALSE)
Envs <- Envs[,1]
Envs

# load hap set final
load(paste(infolder, "/final_HapSet_", TRAIT, ".RData", sep =""))
str(info_qtls_sig_95_final)

info_qtls_sig_95_all <- info_qtls_sig_95_final
HapSet_info_all <- info_qtls_all_final

n_env_all <- length(Envs[which(Envs != "Across")])
n_env_all

#####
#####
#####	generate overview plot
#####
#####

str(info_qtls_sig_95_all)

white_space_size <- 25000000
max_cex <- 5.5

# define centromere positions (B73v4)
Start_bp <- c(136.77, 95.51, 85.78, 109.07, 104.54, 52.3, 56.38, 50.53, 53.75, 57.36, 51.39)
End_bp <- c(137.12, 97.49, 86.93, 110.5, 106.82, 53.11, 56.68, 52.07, 55.39, 57.76, 52.78)
centromers <- cbind(Start_bp, End_bp)
chr_size <- c(307.041717, 244.442276, 235.667834, 246.994605, 223.902240, 174.033170, 182.381542, 181.122637, 159.769782, 150.982314)
chr_size_temp <- c(0, chr_size)*1000000
chr_axis_at <- NULL
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
chr_axis_at <- c(chr_axis_at, mean(c(sum(chr_size_temp[1:CHR]) - white_space_size, sum(chr_size_temp[1:(CHR+1)]))))
chr_ends <- c(chr_ends, sum(chr_size_temp[1:CHR]) - white_space_size, sum(chr_size_temp[1:CHR]))
}
chr_ends <- c(chr_ends, sum(chr_size_temp[1:11]) - white_space_size)
names(chr_add) <- 1:10

# favorable color "blue"
# unfavorable color "red"
# adjust for each trait

	# highlight all markers with a positive effect (agronomical sense, not mathematically)
	if(gsub("BLUE.", "", TRAIT) %in% gsub("BLUE.", "", c(	'EV_V3', 'EV_V4', 'EV_V6',
															'PH_V3', 'PH_V4', 'PH_V6',
															'PH_final',
															'EarTasDist',
															'Purp_Stm_V3', 'Purp_Lvs_V3',
															'NDL_V4', 'NDL_V6',
															'TasLength', 'SpikeLength', 'NumBranch', 'BranchSpace', 'SpikeTasRatio', 'TasAngle',
															'CT',
															'FMY_V4', 'DMY_V4', 'FMY_V6', 'DMY_V6',
															'SPAD_V3', 'SPAD_V4', 'SPAD_V6',
															'Fv_Fm_V4', 'Fv_Fm_V6', 'Rust',
															'DMC', 'DMY',
															'WSC',
															'Purp_Coloring',
															'Rust'))){
	
		col_pos <- rgb(0,0,1, alpha = 0.55)
		col_neg <- rgb(1,0,0, alpha = 0.55)
	
	} else {
	#EH, DtSILK, DtTAS, ASI, RL_R6, TILL, Drought.heat_stress, RL_1, SL,

		col_pos <- rgb(1,0,0, alpha = 0.55)
		col_neg <- rgb(0,0,1, alpha = 0.55)
	
	}

	png(paste(outfolder, "/finalQTLset_all_", TRAIT , "_EffSD_FDR", FDR, "_p", p_thresh, "_CI95_forPaper.png", sep =""), width = 4000 * 3/2, height = 500 + (1600/11 * n_env_all)  * 3/2, res = 200  * 3/2)
	par(mar = c(5, 7, 4, 2) + 0.1)
	
	plot(x = -99999999, y = -99999999, xlim = c(50000000, sum(chr_size_temp) + 150000000), ylim = c(0.5, (length(which(Envs != "Across")))+0.5), main = paste(TRAIT, " QTLs all (FDR", FDR, "%)", sep = ""), ylab = "", xlab = "Chromosome", xaxt = "n", yaxt = "n", cex.axis = cex.axis, cex.lab = cex.lab, cex = cex)
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(0.95,0.95,0.95, alpha = 1))
	axis(side = 1, at = chr_axis_at, labels = 1:10, cex.axis = cex.axis, cex.lab = cex.lab, cex = cex)
	axis(side = 2, at = 1 : length(which(Envs != "Across")), labels = Envs[which(Envs != "Across")], cex.axis = cex.axis, cex.lab = cex.lab, cex = cex, las = 1)
	abline(h = seq(0.5, length(which(Envs != "Across"))+0.5, 1))
	abline(v = chr_ends)
	for(i in 1:length(white_space)){
		polygon(x = c(white_space[[i]][1], white_space[[i]][1], white_space[[i]][2], white_space[[i]][2]), y = c(-100, 100, 100, -100), col = "white", border = NA)	
	}	
	polygon(x = c(- 99999999, - 99999999, 9999999999999, 9999999999999), y = c(-100, 0.5, 0.5, -100), col = "white", border = NA)
	polygon(x = c(- 99999999, - 99999999, 9999999999999, 9999999999999), y = c(100, length(which(Envs != "Across"))+0.5, length(which(Envs != "Across"))+0.5, 100), col = "white", border = NA)
			
		x_pos_all <- NULL
		for(env_i in Envs[which(Envs != "Across")]){
			Y_temp <- read.table(paste(TRAIT, "/", list.files(path = TRAIT, pattern = ".pheno.txt"), sep = ""), stringsAsFactors = FALSE)
			colnames(Y_temp) <- Envs
			Y_temp <- Y_temp[ , env_i]
			res <- info_qtls_sig_95_all[which(info_qtls_sig_95_all$Environment == env_i), ]
			x_pos <- res$Pos_bp
			x_pos <- x_pos + chr_add[res$Chr]
			#points(x = x_pos, y = y_pos, col = col_eff, cex = cex_scale, pch = 19)
			x_pos_all <- c(x_pos_all, x_pos)
		}
		x_pos_all <- unique(x_pos_all)
		abline(v = x_pos_all, lty = 1, lwd = 0.45, col = "gray20")

	polygon(x = c(- 99999999, - 99999999, 9999999999999, 9999999999999), y = c(-100, 0.5, 0.5, -100), col = "white", border = NA)
	polygon(x = c(- 99999999, - 99999999, 9999999999999, 9999999999999), y = c(100, length(which(Envs != "Across"))+0.5, length(which(Envs != "Across"))+0.5, 100), col = "white", border = NA)
		
		sd_env <- NULL
		for(env_i in Envs[which(Envs != "Across")]){
			Y_temp <- read.table(paste(TRAIT, "/", list.files(path = TRAIT, pattern = ".pheno.txt"), sep = ""), stringsAsFactors = FALSE)
			colnames(Y_temp) <- Envs
			Y_temp <- Y_temp[ , env_i]
			res <- info_qtls_sig_95_all[which(info_qtls_sig_95_all$Environment == env_i), ]
			x_pos <- res$Pos_bp
			x_pos <- x_pos + chr_add[res$Chr]
			eff <- res$Qeff / sd(Y_temp, na.rm = TRUE)
			sd_env <- c(sd_env, sd(Y_temp, na.rm = TRUE))
			names(sd_env)[length(sd_env)] <- env_i
			col_eff <- ifelse(eff > 0, col_pos, col_neg)
			# for making the plot look nicer limit cex up to 1 sd effect size
			eff[which(eff < (-1))] <- -1
			eff[which(eff > ( 1))] <-  1
			cex_scale <- abs(eff) * max_cex
			y_pos <- rep(which(Envs[which(Envs != "Across")] == env_i), length(x_pos))
			y_pos <- y_pos + ifelse(col_eff == rgb(0,0,1, alpha = 0.55), +0.25, -0.25)
			points(x = x_pos, y = y_pos, col = col_eff, cex = cex_scale, pch = 19)
		}
	
	# legend
	leg_lab <- c(">1.000", " 0.750", " 0.500", " 0.250", " 0.125")
	leg <- c(1, 0.75, 0.50, 0.25, 0.125)
	cex_scale <- leg * max_cex
	y_leg <- seq((length(Envs[which(Envs != "Across")])/(1 + 1/11*length(Envs[which(Envs != "Across")])))/length(leg),length(Envs[which(Envs != "Across")])/(1 + 1/11*length(Envs[which(Envs != "Across")])), (length(Envs[which(Envs != "Across")])/(1 + 1/11*length(Envs[which(Envs != "Across")])))/length(leg)) + 0.5
	points(x = rep(sum(chr_size_temp) + 70000000, length(leg)), y = y_leg, cex = cex_scale[length(leg):1], col = "black", pch = 19)
	text(x = rep(sum(chr_size_temp) + 170000000, length(leg)), y = y_leg[length(leg):1], cex = 1.3, col = "black", labels = leg_lab)
	text(x = sum(chr_size_temp) + 120000000, y = y_leg[length(leg)] + (length(Envs[which(Envs != "Across")])/(1 + 1/11*length(Envs[which(Envs != "Across")])))/length(leg), cex = 1.4, labels = "Effect size (sd)")

	dev.off()

sd_env
#####
##### plot CI (95%) for specific markers
#####
	load(paste("FavUnfav_Stability/", TRAIT, "/haps_all_", TRAIT, ".RData", sep = ""))
	haps_Fav_info
	
	# ylim has to be adjusted accordingly
	# ylim_all <- c(-0.914, +0.914)

	if(nrow(haps_Fav_info)>0){
	  if(is.null(hap1)){
	  # if no haplotype is specified, plot the most "favorable" marker (n_environments)
	  temp <- haps_Fav_info[order(haps_Fav_info$n_env_fav, haps_Fav_info$max_eff_sd, decreasing = TRUE), ]
	  mk <- temp$qtl_i[1]
	  } else {
	  mk <- hap1
	  }
	Qeff <- info_qtls_all_final$Qeff[which(info_qtls_all_final$Marker == mk)]
	seQeff <- info_qtls_all_final$seQeff[which(info_qtls_all_final$Marker == mk)]
	lower <- Qeff - qnorm(0.975, mean = 0, sd = 1) * seQeff
	upper <- Qeff + qnorm(0.975, mean = 0, sd = 1) * seQeff	
	labE <- unique(info_qtls_all_final$Environment[which(info_qtls_all_final$Marker == mk)])
	sd_env <- sd_env[labE]
	Qeff <- Qeff / sd_env
	lower <- lower / sd_env
	upper <- upper / sd_env
	CI <- data.frame(labE, Qeff, lower, upper)

	color <- rep("gray40", length(Qeff))
	color[which(Qeff > 0 & lower > 0)] <- "blue"
	color[which(Qeff < 0 & upper < 0)] <- "red"

		png(paste(outfolder, "/", mk, "_", TRAIT , "_CI95.png", sep =""), width = 2100, height = 2400, res = 300)
		par(mar = c(5, 7, 4, 2) + 0.1)
		print(
		ggplot(CI, aes(x = c(1:length(labE)), y = Qeff)) + 
		  #ylim(ylim_all) +
		  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, size = 1.3) +
		  geom_point(colour=color, size = 5) +
		  theme(text = element_text(size = 15),
			axis.title.x = element_text(size = 20),
			axis.text.x = element_text(size = 16, angle = 90),
			axis.text.y = element_text(size = 16),
			axis.title.y = element_text(size = 20)) + 
		  xlab("Environment") +
		  scale_x_continuous(breaks = c(1:length(labE)), labels = labE) +
		  geom_hline(yintercept = 0, col="orange3", linetype = "dashed", size = 1.8) +
		  ylab(paste(TRAIT, " (SD)", sep = ""))
		  )
		dev.off()
	}

#####

		if(nrow(haps_Unfav_info)>0){
		if(is.null(hap2)){
		  # if no haplotype is specified, plot the most "unfavorable" marker (n_environments)
		  temp <- haps_Unfav_info[order(haps_Unfav_info$n_env_unfav, abs(haps_Unfav_info$min_eff_sd), decreasing = TRUE), ]
		  mk <- temp$qtl_i[1]
		} else {
		  #####	
		  mk <- hap2
		}
	Qeff <- info_qtls_all_final$Qeff[which(info_qtls_all_final$Marker == mk)]
	seQeff <- info_qtls_all_final$seQeff[which(info_qtls_all_final$Marker == mk)]
	lower <- Qeff - qnorm(0.975, mean = 0, sd = 1) * seQeff
	upper <- Qeff + qnorm(0.975, mean = 0, sd = 1) * seQeff	
	labE <- unique(info_qtls_all_final$Environment[which(info_qtls_all_final$Marker == mk)])
	sd_env <- sd_env[labE]
	Qeff <- Qeff / sd_env
	lower <- lower / sd_env
	upper <- upper / sd_env
	CI <- data.frame(labE, Qeff, lower, upper)

	color <- rep("gray40", length(Qeff))
	color[which(Qeff > 0 & lower > 0)] <- "blue"
	color[which(Qeff < 0 & upper < 0)] <- "red"

		png(paste(outfolder, "/", mk, "_", TRAIT , "_CI95.png", sep =""), width = 2100, height = 2400, res = 300)
		par(mar = c(5, 7, 4, 2) + 0.1)
		print(
		ggplot(CI, aes(x = c(1:length(labE)), y = Qeff)) + 
		  #ylim(ylim_all) +
		  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, size = 1.3) +
		  geom_point(colour=color, size = 5) +
		  theme(text = element_text(size = 15),
			axis.title.x = element_text(size = 20),
			axis.text.x = element_text(size = 16, angle = 90),
			axis.text.y = element_text(size = 16),
			axis.title.y = element_text(size = 20)) + 
		  xlab("Environment") +
		  scale_x_continuous(breaks = c(1:length(labE)), labels = labE) +
		  geom_hline(yintercept = 0, col="orange3", linetype = "dashed", size = 1.8) +
		  ylab(paste(TRAIT, " (SD)", sep = ""))
		  )
		dev.off()
}
	
#####

	if(nrow(haps_Inter_info)>0){
		if(is.null(hap3)){
		  # if no haplotype is specified, plot the most "GxE" marker (n_environments)
		  temp <- haps_Inter_info[order(haps_Inter_info$n_env_neut, decreasing = FALSE), ]
		  mk <- temp$qtl_i[1]
		} else {
		  #####	
		  mk <- hap3
		}
	Qeff <- info_qtls_all_final$Qeff[which(info_qtls_all_final$Marker == mk)]
	seQeff <- info_qtls_all_final$seQeff[which(info_qtls_all_final$Marker == mk)]
	lower <- Qeff - qnorm(0.975, mean = 0, sd = 1) * seQeff
	upper <- Qeff + qnorm(0.975, mean = 0, sd = 1) * seQeff	
	labE <- unique(info_qtls_all_final$Environment[which(info_qtls_all_final$Marker == mk)])
	sd_env <- sd_env[labE]
	Qeff <- Qeff / sd_env
	lower <- lower / sd_env
	upper <- upper / sd_env
	CI <- data.frame(labE, Qeff, lower, upper)
	summary(CI)
	color <- rep("gray40", length(Qeff))
	color[which(Qeff > 0 & lower > 0)] <- "blue"
	color[which(Qeff < 0 & upper < 0)] <- "red"

		png(paste(outfolder, "/", mk, "_", TRAIT , "_CI95.png", sep =""), width = 2100, height = 2400, res = 300)
		par(mar = c(5, 7, 4, 2) + 0.1)
		print(
		ggplot(CI, aes(x = c(1:length(labE)), y = Qeff)) + 
		  # ylim(ylim_all) +
		  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, size = 1.3) +
		  geom_point(colour=color, size = 5) +
		  theme(text = element_text(size = 15),
			axis.title.x = element_text(size = 20),
			axis.text.x = element_text(size = 16, angle = 90),
			axis.text.y = element_text(size = 16),
			axis.title.y = element_text(size = 20)) + 
		  xlab("Environment") +
		  scale_x_continuous(breaks = c(1:length(labE)), labels = labE) +
		  geom_hline(yintercept = 0, col="orange3", linetype = "dashed", size = 1.8) +
		  ylab(paste(TRAIT, " (SD)", sep = ""))
		  )
		dev.off()
		}
		
