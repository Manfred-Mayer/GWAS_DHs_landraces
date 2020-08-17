###################################################################
###################################################################
####
#### for a given window
#### plot the frequencies of all haplotypes within DHs and BLs
#### plot the effects of all alternative haplotypes relative to the focus haplotype (per environment)
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
library(plotrix)

# arguments
nSNPs <- 10
steps <- 10
p_thresh <- 0.01
FDR <- "15"
minCount <- 3

# choose according to the focus haplotype you want to plot
TRAIT <- "PH_V4"
hap_focus <- "wind_03_01592_6"

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
TRAIT
hap_focus

# graphical parameters
cex <- 0.6
cex.axis <- 1.6
cex.lab <- 1.6

# generate letter vector for naming haplotype variants
letters_MM <- c(letters[1:26], paste(rep(letters[1:26], rep(26,26)), rep(letters[1:26], 26), sep = ""))
letters_MM <- letters_MM[-1]
letters_MM

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##############
##############

wind <- substr(hap_focus, 1, 13)
wind
region <- paste0(TRAIT, ".", hap_focus)
region

# input folder
infolder <- paste("wholeWindModel/", TRAIT, "/", region, sep = "")
infolder

# output folder
outfolder <- paste("wholeWindModel/", TRAIT, "/", region, sep = "")
outfolder

#################################################################################################################
# load DH data
load(paste("geno_InfoList_DH_m", nSNPs, "s", steps, ".RData", sep = ""))
str(InfoList)
str(geno)

# generate haplotype vector for the respective window (coded as factor)
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

# environments for that trait
Envs <- read.table(paste0(TRAIT, "/", list.files(path = TRAIT, pattern = ".phenotypes_order.txt")), stringsAsFactors = FALSE)
Envs <- Envs[,1]
# in model we only need the within environment blues, so exclude the across environment blues
Envs <- Envs[-which(Envs == "Across")]
Envs

# load results from 1.11_modelWholeWind
load(paste(infolder, "/info_", region, "_v3.RData", sep = ""))
str(info)
info <- info[grep(wind, rownames(info)),]
str(info)
summary(info$HAPeff_sd)

# load check genotypic data
# SNP data of BL
load("gpBL.RData")
str(gpBL)

# SNP data of LR
load(paste("gpDH.RData", sep=""))
str(gpDH)

	chr_qtl <- InfoList[hap_focus, "Chr"]
	start_qtl <- InfoList[hap_focus, "Pos_Start_bp"]
	end_qtl <- InfoList[hap_focus, "Pos_End_bp"]
	IDs_qtl <- rownames(geno)[which(geno[, hap_focus] == 2)]

	map_i <- gpBL$map[which(gpBL$map$chr == chr_qtl), ]
	map_i <- map_i[which(map_i$pos >= start_qtl), ]
	map_i <- map_i[which(map_i$pos <= end_qtl), ]

	geno_10_LR <- gpDH$geno[ , rownames(map_i)]
	geno_10_BL <- gpBL$geno[ , rownames(map_i)]

	# define groups
	haplotypes <- sort(unique(hap_vec))
	hapGroups <- matrix(NA, nrow = length(haplotypes), ncol = nrow(map_i))
	colnames(hapGroups) <- rownames(map_i)
	rownames(hapGroups) <- haplotypes
	for(hap_i in rownames(hapGroups)){
		exampleID <- names(hap_vec)[which(hap_vec == hap_i)][1]
		hapGroups[hap_i, ] <- geno_10_LR[exampleID, ]
	}
	hapGroups

	# sequences
	seq_groups_10 <- apply(hapGroups, 1, paste0, collapse = "")
	seq_groups_10

	seq_10_LR <- apply(geno_10_LR, 1, paste0, collapse = "")
	seq_10_BL <- apply(geno_10_BL, 1, paste0, collapse = "")
	
	# generate groupList
	groupList <- list()
	prop_LR <- NULL
	prop_BL <- NULL
	for(hap_i in rownames(hapGroups)){
		IDs_LR <- names(seq_10_LR)[which(seq_10_LR == seq_groups_10[hap_i])]
		IDs_BL <- names(seq_10_BL)[which(seq_10_BL == seq_groups_10[hap_i])]
		prop_LR <- c(prop_LR, length(IDs_LR) / length(seq_10_LR))
		prop_BL <- c(prop_BL, length(IDs_BL) / length(seq_10_BL))
		groupList[[hap_i]] <- c(IDs_LR, IDs_BL)
	}
	
	# generate input for 10 SNP window plot
	markers_10 <- colnames(gpDH$geno)[(which(colnames(gpDH$geno) %in% rownames(map_i))[1]) : (which(colnames(gpDH$geno) %in% rownames(map_i))[10])]
	geno_both_10 <- rbind(gpDH$geno[ , markers_10], gpBL$geno[ , markers_10])
	map_10 <- gpBL$map[markers_10, ]
	a_10 <- geno_both_10[IDs_qtl[1],]
		
	
	
###
### plotting
###

# function REF/ALT allele
RefAlt_f <- function(x) {
  y1 <- x[-1]
  y2 <- y1
  y2[which(y1 == x[1])] <- "REF"
  y2[which(y1 != x[1])] <- "ALT"
  return(y2)
}

# generate input data for plotting nucleotide sequences
col_ref <- "gray35"
col_alt <- "gray85"

a_10 <- geno_both_10[IDs_qtl,]
a_10_temp <-  apply(a_10, 1, paste0, collapse = "")
a_10_temp2 <- sort(table(a_10_temp), decreasing = TRUE)
a_10_temp <- names(a_10_temp)[which(a_10_temp == names(a_10_temp2)[1])]
a_10 <- a_10[a_10_temp[1],]
plot_list_seq <- list()
for(hap_i in haplotypes){
  geno_i <- geno_both_10[groupList[[hap_i]], ]
  seq_10_i <- apply(geno_i, 1, paste0, collapse = "")
  freq_10_i <- sort(table(seq_10_i), decreasing = TRUE)
  # order according to frequency
  if(length(freq_10_i) > 1){
    order_vec_i <- NULL
    for(freq_i in names(freq_10_i)){
      order_vec_i <- c(order_vec_i, names(seq_10_i)[which(seq_10_i == freq_i)])
    }
    geno_i <- geno_i[order_vec_i, ]
  }
  
  geno_temp_i <- rbind(a_10, geno_i)
  geno_RefAlt <- apply(geno_temp_i, 2, RefAlt_f)
  geno_RefAlt[which(geno_RefAlt == "REF")] <- col_ref
  geno_RefAlt[which(geno_RefAlt == "ALT")] <- col_alt
  
  plot_list_seq[[hap_i]] <- geno_RefAlt
}

##
## plotting
##
ylim <- c(-0.7,length(haplotypes) + 1.5)
figure_end_x <- 21.5
figure_start_x <- 7.7
xlim <- c(figure_start_x, figure_end_x)
ns <- "ns"
eff_cols <- c("#BB4444","white","#4477AA")
n_col <- 61
range_col <- c(-max(abs(info$HAPeff_sd)), +max(abs(info$HAPeff_sd)))
range_col
cex_text <- 0.95
cex_text2 <- 0.8

png(paste(outfolder, "/HapPlot_Eff", region, "_hapTable.png", sep = ""), width = 1500, height = 2100, res = 300)
par(mar = c(0, 0, 0, 0))
cols <- colorRampPalette(c("#BB4444","white","#4477AA"))(n_col)
breaks_col <- seq(range_col[1], range_col[2], (range_col[2] - range_col[1]) / n_col)
breaks_col[1] <- breaks_col[1] - 0.00000001

plot(x = -10, y = -10,
     xlim = xlim, ylim = ylim,
     bty = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

# add percentage of LR and BL
for(seq_i in 1:length(haplotypes)){
  text(x = 8, y = length(haplotypes) - seq_i + 1, pos = 2,
       labels = LETTERS[seq_i], cex = 1.3)
  temp <- as.character(round(prop_LR[seq_i] * 100, digits = 1))
  if(length(unlist(strsplit(temp, split = ".", fixed = TRUE))) == 1){
	temp <- paste(temp, ".0", sep = "")
  }
  text(x = 9.6, y = length(haplotypes) - seq_i + 1, pos = 2,
       labels = temp, cex = 1.3)
  temp <- as.character(round(prop_BL[seq_i] * 100, digits = 1))
  if(length(unlist(strsplit(temp, split = ".", fixed = TRUE))) == 1){
	temp <- paste(temp, ".0", sep = "")
  }
  text(x = 11.3, y = length(haplotypes) - seq_i + 1, pos = 2,
       labels = temp, cex = 1.3)
  abline(h = length(haplotypes) - seq_i + 0.5, lty = 2)
}

###
### add hap effects
###
x_eff_start <- 11.5
order_ENV <- c("2017.BBG",
               "2017.EIN",
               "2017.GOL",
               "2017.OLI",
               "2017.ROG",
               "2017.TOM",
               "2018.EIN",
               "2018.GOL",
               "2018.KLW",
               "2018.ROG",
               "2018.TOM")
order_ENV <- order_ENV[which(order_ENV %in% Envs)]
order_ENV
steps_eff <- (figure_end_x-x_eff_start) / length(order_ENV)

for(hap_i in 2:length(haplotypes)){
 temp_info <- info[which(info$haplotype == haplotypes[hap_i]),]
 if(nrow(temp_info) > 0){
	 rownames(temp_info) <- temp_info$ENV
	 temp_info <- temp_info[order_ENV, ]
	 x1_temp <- seq(x_eff_start, x_eff_start + steps_eff * nrow(temp_info)- steps_eff, steps_eff)
	 x2_temp <- seq(x_eff_start + steps_eff, x_eff_start + steps_eff * nrow(temp_info), steps_eff)
	 x_num <- temp_info$HAPeff_sd
	 print(summary(x_num))
	 x_col <- x_num
	 for(i in 1:length(cols)){
	   x_col[which((x_num > breaks_col[i]) & (x_num <= breaks_col[i+1]))] <- cols[i]
	 }
	 y_borders <- c(ylim[2] - 1 - (hap_i - 1) - 0.9, ylim[2] - 1 - (hap_i - 1) - 0.1)
	 for(i in 1:length(x_col)){
	   polygon(x = c(x1_temp[i], x1_temp[i], x2_temp[i], x2_temp[i]),
			   y = c(y_borders[1], y_borders[2], y_borders[2], y_borders[1]),
			   col = x_col[i],
			   border = NA)
	   if(temp_info$sig[i] == "ns"){
		 text(x = mean(c(x1_temp[i], x2_temp[i])), y = mean(y_borders),
			  labels = "ns", cex = cex_text2, col = "gray30")
	   }
	 }
 }
}

# add color key
breaks <- c(breaks_col)
colkey_range <- range_col
cols_heatmap <- c(cols)
xstart <- x_eff_start + 2
xend <- figure_end_x - 2
ystart <- length(haplotypes) + 1
yend <- length(haplotypes) + 1.4
color.legend(xl = xstart, yb = ystart, xr = xend, yt = yend, rect.col = cols_heatmap, gradient = "x", legend = "")
points(x = c(xstart, xstart), y = c(ystart, ystart - (yend - ystart)*0.10), type = "l", lwd = 1.5)
points(x = rep(mean(c(xstart, xend)), 2), y = c(ystart, ystart - (yend - ystart)*0.10), type = "l", lwd = 1.5)
points(x = c(xend, xend), y = c(ystart, ystart - (yend - ystart)*0.10), type = "l", lwd = 1)
text(x = c(xstart), y = c(ystart), labels = round(range_col[1], digits = 2), pos = 1, cex = cex_text)
text(x = rep(mean(c(xstart, xend)), 1), y = c(ystart), labels = "0.00", pos = 1, cex = cex_text)
text(x = c(xend), y = c(ystart), labels = round(range_col[length(range_col)], digits = 2), pos = 1, cex = cex_text)

dev.off()

##############
############################################################################################################################################################################################################




