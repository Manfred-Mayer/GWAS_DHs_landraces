###################################################################
###################################################################
####
#### compare haplotype frequencies between DH lines and breeding lines
#### Plot frequencies in 2D plot
#### identify number of haplotypes present in DH but absent in BL and vice versa 
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
library(plot3D)
library(VennDiagram)

# arguments
nSNPs <- 10
steps <- 10

###################################################################

dir.create("comp_BL_DHs")
outFolder <- paste("comp_BL_DHs", sep = "")
outFolder

# load haplotype data
load(paste("geno_InfoList_DH.BL_m", nSNPs, "s", steps, ".RData", sep = ""))
str(geno)

#
# devide into landrace panel and elite panel
#

geno_BL <- geno[which(substr(rownames(geno), 1, 3) == "EL_"), ]
geno_DH <- geno[which(substr(rownames(geno), 1, 3) == "DH_"), ]
str(geno_BL)
str(geno_DH)

# here, count every haplotype
minCount <- 1

#
# check for each haplotype in which landraces it occurs
#

#
# DH lines
#
landraces_DH <- unique(substr(rownames(geno_DH), 1, 5))
landraces_DH

# function for calculating how often a haplotype occurs per landrace
occur_n_perLR_f <- function(x, template_in, template_out){
	carriers <- template_in[which(x == 2)]
	n_perLR <- table(carriers)
	template_out[names(n_perLR)] <- n_perLR
	return(template_out)
}

template_out <- rep(0, length(landraces_DH))
names(template_out) <- landraces_DH
template_in <- substr(rownames(geno_DH), 1, 5)

info_occur_DH <- apply(geno_DH, 2, occur_n_perLR_f, template_in = template_in, template_out = template_out)
str(info_occur_DH)

#
# calculate mean proportion in the 3 DH libraries per haplotype
#
# prop function
LR_prop_f <- function(x, n_LR){
	return(x / n_LR)
}
n_LR <- table(substr(rownames(geno_DH), 1, 5))
n_LR
all.equal(names(n_LR), rownames(info_occur_DH))
prop_perLR_DH <- apply(info_occur_DH, 2, LR_prop_f, n_LR = n_LR)

mean_freq_perLR_DH <- apply(prop_perLR_DH, 2, mean)
summary(mean_freq_perLR_DH)
length(which(mean_freq_perLR_DH == 0)) / length(mean_freq_perLR_DH)
length(which(mean_freq_perLR_DH < 0.01)) / length(mean_freq_perLR_DH)

#
# compare frequency in landrace panel with frequency in elite panel
#
prop_f <- function(x){
	length(which(x==2)) / length(x)
}
freq_BL <- apply(geno_BL , 2, prop_f)
summary(freq_BL)
length(which(freq_BL == 0)) / length(freq_BL)

freq_LR_DH <- apply(geno_DH , 2, prop_f)
summary(freq_LR_DH)
length(which(freq_LR_DH == 0)) / length(freq_LR_DH)


cor_results <- matrix(NA, nrow = 2, ncol = 4)
colnames(cor_results) <- c(	"freq_LR.DH_pearson",
							"mean_freq_perLR.DH_pearson",
							"freq_LR.DH_spearman",
							"mean_freq_perLR.DH_spearman"
							)
rownames(cor_results) <- c("estimate", "p_value")
cor_results

print(all.equal(names(freq_BL), names(freq_LR_DH)))
cor_temp <- cor.test(x = freq_BL, y = freq_LR_DH, method = c("pearson"))
cor_results["estimate", "freq_LR.DH_pearson"] <- cor_temp$estimate
cor_results["p_value", "freq_LR.DH_pearson"] <- cor_temp$p.value
cor_temp <- cor.test(x = freq_BL, y = freq_LR_DH, method = c("spearman"))
cor_results["estimate", "freq_LR.DH_spearman"] <- cor_temp$estimate
cor_results["p_value", "freq_LR.DH_spearman"] <- cor_temp$p.value

print(all.equal(names(freq_BL), names(mean_freq_perLR_DH)))
cor_temp <- cor.test(x = freq_BL, y = mean_freq_perLR_DH, method = c("pearson"))
cor_results["estimate", "mean_freq_perLR.DH_pearson"] <- cor_temp$estimate
cor_results["p_value", "mean_freq_perLR.DH_pearson"] <- cor_temp$p.value
cor_temp <- cor.test(x = freq_BL, y = mean_freq_perLR_DH, method = c("spearman"))
cor_results["estimate", "mean_freq_perLR.DH_spearman"] <- cor_temp$estimate
cor_results["p_value", "mean_freq_perLR.DH_spearman"] <- cor_temp$p.value

cor_results

res_list <- list()
res_list[["freq_BL"]] <- freq_BL
res_list[["freq_LR_DH"]] <- freq_LR_DH
res_list[["mean_freq_perLR_DH"]] <- mean_freq_perLR_DH

save(res_list, cor_results, info_occur_DH, prop_perLR_DH, file = paste(outFolder, "/cor_results_DHvsBL_DHs_ALLhaps_minCount", minCount, ".RData", sep = ""))
write.table(cor_results, file = paste(outFolder, "/cor_results_DHvsBL_ALLhaps_minCount", minCount, ".csv", sep =""), sep = ";" , dec=".", row.names = T, col.names = NA)

#
# number of haplotypes per group
#
present_DH <- names(res_list$freq_LR_DH)[which(res_list$freq_LR_DH > 0)]
n_DH <- length(present_DH)
present_BL <- names(res_list$freq_BL)[which(res_list$freq_BL > 0)]
n_BL <- length(present_BL)
n_DH
n_BL

#
# overlapping
#
n_DH.BL <- length(present_BL[which(present_BL %in% present_DH)])
n_DH.BL

##
## generate Venn Diagram
##

# define colors
col2rgb_MM <- function(x, alpha){
  rgb(x[1],x[2],x[3], alpha = alpha)
}

temp <- c(0.6,0.6,0.9)
BL.2_fill <- col2rgb_MM(temp, alpha = 0.5)

temp <- c(0,0.85,0)
LR.2_fill <- col2rgb_MM(temp, alpha = 0.5)

venn.plot <- draw.pairwise.venn(area1=n_DH, area2=n_BL, cross.area=n_DH.BL, category = c("LR","BL"),cex=1.9,
                                cat.dist=c(rep(-0.03,2)), col = rep("black", 2))
dev.off()

# manually set numbers with thousand separator
venn.plot[[5]][["label"]]
venn.plot[[5]][["label"]] <- "100k"
venn.plot[[6]][["label"]]
venn.plot[[6]][["label"]] <- "94k"
venn.plot[[7]][["label"]]
venn.plot[[7]][["label"]] <- "263k"

# manually set colorsr
BL_col <- BL.2_fill
LR_col <- LR.2_fill
venn.plot[[3]][["gp"]]$fill
venn.plot[[3]][["gp"]]$fill <- BL_col
venn.plot[[3]][["gp"]]$col <- rgb(0,0,0, alpha = 0.2)
venn.plot[[4]][["gp"]]$fill
venn.plot[[4]][["gp"]]$fill <- LR_col
venn.plot[[4]][["gp"]]$col <- rgb(0,0,0, alpha = 0.2)

# manually set cex label
venn.plot[[5]][["gp"]]$cex <- 2.3
venn.plot[[6]][["gp"]]$cex <- 2.3
venn.plot[[7]][["gp"]]$cex <- 2.3
venn.plot[[8]][["gp"]]$cex <- 2.6
venn.plot[[9]][["gp"]]$cex <- 2.6

# manually set font
venn.plot[[5]][["gp"]]$fontfamily <- "sans"
venn.plot[[6]][["gp"]]$fontfamily <- "sans"
venn.plot[[7]][["gp"]]$fontfamily <- "sans"
venn.plot[[8]][["gp"]]$fontfamily <- "sans"
venn.plot[[9]][["gp"]]$fontfamily <- "sans"

tiff(paste(outFolder, "/Venn_DHvsBL_haplotypes.tiff", sep = ""), width=1500, height=1500, res=300)
grid.draw(venn.plot)
dev.off()





# in DH panel
# present
length(present_BL[which(present_BL %in% present_DH)])
length(present_BL[which(present_BL %in% present_DH)]) / length(present_BL)
# absent
length(present_BL[which(present_BL %in% present_DH == FALSE)])
length(present_BL[which(present_BL %in% present_DH == FALSE)]) / length(present_BL)

# in BL panel
# present
length(present_DH[which(present_DH %in% present_BL)])
length(present_DH[which(present_DH %in% present_BL)]) / length(present_DH)
# absent
length(present_DH[which(present_DH %in% present_BL == FALSE)])
length(present_DH[which(present_DH %in% present_BL == FALSE)]) / length(present_DH)

##
## analyze haplotypes absent in elite panel but present in the DH panel
##

# which are absent
absent_BL <- names(res_list$freq_BL)[which((res_list$freq_BL == 0) & (res_list$freq_LR_DH > 0))]
length(absent_BL)
length(present_DH)
length(absent_BL) / length(present_DH)

# what is the average frequency of these haplotypes in the landrace panel
summary(res_list$freq_LR_DH[absent_BL])
summary(res_list$freq_LR_DH)

# how many of those haplotypes are present in 3 landraces
str(info_occur_DH)
info_occur_DH_absentBL <- info_occur_DH[,absent_BL]
str(info_occur_DH_absentBL)
check_3_f <- function(x){
  check <- any(x == 0) == FALSE
  return(check)
}
check_3 <- apply(info_occur_DH_absentBL, 2, check_3_f)
length(which(check_3))

# how many of those haplotypes are present in 2 landraces
check_2_f <- function(x){
  check <- length(which(x == 0)) == 1
  return(check)
}
check_2 <- apply(info_occur_DH_absentBL, 2, check_2_f)
length(which(check_2))

# how many of those haplotypes are present in 1 landrace
check_1_f <- function(x){
  check <- length(which(x == 0)) == 2
  return(check)
}
check_1 <- apply(info_occur_DH_absentBL, 2, check_1_f)
length(which(check_1))

# from those, how many occur in each landrace
check_LR_f <- function(x, LR){
  check <- (length(which(x[-LR] == 0)) == 2) & (x[LR] != 0)
  return(check)
}
check_KE <- apply(info_occur_DH_absentBL, 2, check_LR_f, LR = 1)
check_LL <- apply(info_occur_DH_absentBL, 2, check_LR_f, LR = 2)
check_PE <- apply(info_occur_DH_absentBL, 2, check_LR_f, LR = 3)
length(which(check_KE))
length(which(check_LL))
length(which(check_PE))
length(which(check_KE)) + length(which(check_LL)) + length(which(check_PE))

# check frequencies within those landraces
str(prop_perLR_DH)
freqs_single_KE <- prop_perLR_DH[1, absent_BL[check_KE]]
freqs_single_LL <- prop_perLR_DH[2, absent_BL[check_LL]]
freqs_single_PE <- prop_perLR_DH[3, absent_BL[check_PE]]
summary(freqs_single_KE)
summary(freqs_single_LL)
summary(freqs_single_PE)

# what is the average frequency of those haplotypes within the landraces
prop_perLR_absentAll_allFreqs_list <- apply(prop_perLR_DH[ , absent_BL], 2, function(x) {x[which(x != 0)]})

prop_perLR_absentAll_allFreqs <- unlist(prop_perLR_absentAll_allFreqs_list)
str(prop_perLR_absentAll_allFreqs)
summary(prop_perLR_absentAll_allFreqs)

# what is the maximum frequency of those haplotypes within a landrace
max_per_LR_absentAll_allFreqs <- unlist(lapply(prop_perLR_absentAll_allFreqs_list, max))
str(max_per_LR_absentAll_allFreqs)
summary(max_per_LR_absentAll_allFreqs)


##
## analyze haplotypes absent in DH panel but present in the BL panel
##

# which are absent
absent_DH <- names(res_list$freq_LR_DH)[which((res_list$freq_LR_DH == 0) & (res_list$freq_BL > 0))]
length(absent_DH)
length(present_BL)
length(absent_DH) / length(present_BL)

# what is the average frequency of these haplotypes in the breeding panel
summary(res_list$freq_BL[absent_DH])

#
# generate heatmap for freq_DH vs freq_BL
#

n_cells_LR <- 65

cor_table <- matrix(NA, ncol = n_cells_LR+1, nrow = 65+1)
rownames(cor_table) <- 0:65
colnames(cor_table) <- 0:n_cells_LR


n_BL <- res_list$freq_BL * 65
table(n_BL)

# make bins for LR
lower <- c(-1, seq(0, n_cells_LR - 1, 1))
upper <- seq(0, n_cells_LR, 1)
bins <- cbind(lower, upper)
n_LR_temp <- res_list$freq_LR_DH * n_cells_LR
n_LR <- n_LR_temp
for(i in 1:nrow(bins)){
  n_LR[which((n_LR_temp > bins[i,1]) & (n_LR_temp <= bins[i,2]))] <- bins[i,2]
}
table(n_LR)

for(row_i in rownames(cor_table)){
  for(col_i in colnames(cor_table)){
    print(paste(row_i, col_i, sep = ":"))
    n_i <- length(which((n_BL == as.numeric(row_i)) & (n_LR == as.numeric(col_i))))
    cor_table[row_i, col_i] <- n_i
  }
}
str(cor_table)
cor_table[1:10,1:6]
sum(cor_table)
length(res_list$freq_LR_DH)
length(res_list$freq_BL)
sum(cor_table[1,])
sum(cor_table[,1])

#
# generate plot
#
z <- t(cor_table)
z

color <- z
z[is.na(z)] <- 0
z
color[which(color == 0)] <- NA

# colors
# define color scheme for the kinship plots
cols <- colorRampPalette(c("#ffffff", "#ffe4c8", "#ffdab4", "#ffd0a0", 
                           "#ffc182", "#ffb76e", "#ffad5a", "#ffa346", "#ff9932", 
                           "#ff8f1e", "#ff850a", "#e17100", "#cd6700", "#b95d00", 
                           "#a55300", "#914900", "#7d3f00", "#5f3000", "#4b2600", 
                           "#371c00", "#000000"))(27000)
cols <- cols[length(cols):1]
#cols <- colorRampPalette(c("black", "white"))(24000)
cols <- cols[c(1:4000,seq(4001,7000,5),seq(7001,9000,10),seq(9001,12000,25), seq(12001,15000,50), seq(15001,17000,100),seq(17001,21000,250),seq(21001,26001,500))]
hist(1:length(cols), breaks = 0:length(cols), col = cols, border = NA)
cols <- cols[length(cols):1]
dev.off()

tiff(paste(outFolder, "/freq_DH_vs_freqBL_2D_allHaps_minCount1.tiff", sep = ""), width=2400, height=1910, res=300)
par(mai=c(1.1,1.1,0.1,1.1), las=1, mgp=c(4,1,0))
image2D(z=color, col=cols, border=NA, x=as.numeric(rownames(z)), y=as.numeric(colnames(z)),
        yaxt="n", xaxt="n", xlab="Landrace derived DH lines", ylab="Breeding lines", cex.axis=1.5, cex.lab=1.5, colkey=list(cex.clab=1, cex.axis=1, at = seq(5000, 75000, 5000), length=1.01))
axis(side = 2, at = seq(5,65,10), labels = round(seq(5,65,10) / 65, digits = 2), padj = 0.5, cex.axis=1.5, las = 2)
axis(side = 1, at = seq(5,n_cells_LR,10), labels = round(seq(5,n_cells_LR,10) / n_cells_LR, digits = 2), padj = 0.5, cex.axis=1.5)
dev.off()

###
###################################################################
