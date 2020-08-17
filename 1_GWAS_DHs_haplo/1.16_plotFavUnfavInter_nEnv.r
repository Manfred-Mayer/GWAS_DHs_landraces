###################################################################
###################################################################
####
#### plot the number of environments in which favorable/unfavorable/interacting haplotypes
#### show significant effects
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
TRAIT <- "PH_final"
p_thresh <- 0.01
FDR <- "15"

# graphical parameters
cex <- 0.6
cex.axis <- 1.6
cex.lab <- 1.6

###################################################################################################
###################################################################################################
###################################################################################################

traits <- c("EV_V4", "EV_V6", "PH_V4", "PH_V6",
			"PH_final",
			"FF", "MF",
			"LO", "TILL")

mean_nEnv <- matrix(NA, nrow = length(traits), ncol = 4)
median_nEnv <- matrix(NA, nrow = length(traits), ncol = 4)
colnames(mean_nEnv) <- c("Favorable", "Unfavorable", "Interacting", "All")
colnames(median_nEnv) <- c("Favorable", "Unfavorable", "Interacting", "All")
rownames(mean_nEnv) <- traits
rownames(median_nEnv) <- traits

# for each trait separately
for(TRAIT in traits){

# input folder
infolder <- paste(TRAIT, "/QTLregs/finalQTL", sep = "")
infolder
# output folder
outfolder <- paste("Stability_GWASplots/", TRAIT, sep = "")
outfolder

load(paste("FavUnfav_Stability/", TRAIT, "/haps_all_", TRAIT, ".RData", sep = ""))
if(length(haps_Unfav_info) > 0){
	n_env_max <- sum(haps_Unfav_info[ 1, 2:4])
} else {
	n_env_max <- sum(haps_Fav_info[ 1, 2:4])
}
n_env_FavUnfavInfo <- rbind(haps_Fav_info, haps_Unfav_info, haps_Inter_info)

png(paste(outfolder, "/boxplot_n_env_", TRAIT , ".png", sep =""), width = 2700, height = 900, res = 300)
par(mar = c(3, 5.9, 0, 0) + 0.1, mgp = c(2,0.6,0))
N_env <- c(	n_env_max - n_env_FavUnfavInfo$n_env_neut,
			n_env_max - haps_Fav_info$n_env_neut,
			n_env_max - haps_Unfav_info$n_env_neut,
			n_env_max - haps_Inter_info$n_env_neut)
categ <- c(	rep("All", nrow(n_env_FavUnfavInfo)),
			rep("Favorable", nrow(haps_Fav_info)),
			rep("Unfavorable", nrow(haps_Unfav_info)),
			rep("Changing sign", nrow(haps_Inter_info)))

dat <- data.frame(N_env = N_env, Category = categ)
dat$Category <- factor(dat$Category, levels = c("All", "Interacting", "Unfavorable", "Favorable"))


boxplot(N_env ~ Category, data = dat, xlab = "N sign. environments", ylab = "", col = c("gray30", "palevioletred2", "red", "blue"), cex.axis = 1, cex.lab = 1.1, horizontal = TRUE, las = 1, ylim = c(1,n_env_max))

mean_nEnv[TRAIT, 1] <- mean(n_env_max - haps_Fav_info$n_env_neut)
mean_nEnv[TRAIT, 2] <- mean(n_env_max - haps_Unfav_info$n_env_neut)
mean_nEnv[TRAIT, 3] <- mean(n_env_max - haps_Inter_info$n_env_neut)
mean_nEnv[TRAIT, 4] <- mean(n_env_max - n_env_FavUnfavInfo$n_env_neut)

points(x = mean_nEnv[TRAIT, ], y = c(4,3,2,1), bg = "gray", pch = 23)

dev.off()

median_nEnv[TRAIT, 1] <- median(n_env_max - haps_Fav_info$n_env_neut)
median_nEnv[TRAIT, 2] <- median(n_env_max - haps_Unfav_info$n_env_neut)
median_nEnv[TRAIT, 3] <- median(n_env_max - haps_Inter_info$n_env_neut)
median_nEnv[TRAIT, 4] <- median(n_env_max - n_env_FavUnfavInfo$n_env_neut)

}

outfolder <- paste("Stability_GWASplots", sep = "")
outfolder
write.table(mean_nEnv, file = paste(outfolder, "/mean_nEnv_perTrait.csv", sep = ""), sep = ";", dec = ".", col.names = NA, row.names = TRUE, quote= FALSE)
write.table(median_nEnv, file = paste(outfolder, "/median_nEnv_perTrait.csv", sep = ""), sep = ";", dec = ".", col.names = NA, row.names = TRUE, quote= FALSE)


#
# all early traits together
#
n_env_FavUnfavInfo_all <- NULL
haps_Fav_info_all <- NULL
haps_Unfav_info_all <- NULL
haps_Inter_info_all <- NULL
for(TRAIT in c("EV_V4", "EV_V6", "PH_V4", "PH_V6")){

load(paste("FavUnfav_Stability/", TRAIT, "/haps_all_", TRAIT, ".RData", sep = ""))
if(length(haps_Unfav_info) > 0){
	n_env_max <- sum(haps_Unfav_info[ 1, 2:4])
} else {
	n_env_max <- sum(haps_Fav_info[ 1, 2:4])
}
n_env_FavUnfavInfo <- rbind(haps_Fav_info, haps_Unfav_info, haps_Inter_info)

n_env_FavUnfavInfo_all <- rbind(n_env_FavUnfavInfo_all, n_env_FavUnfavInfo)
haps_Fav_info_all <- rbind(haps_Fav_info_all, haps_Fav_info)
haps_Unfav_info_all <- rbind(haps_Unfav_info_all, haps_Unfav_info)
haps_Inter_info_all <- rbind(haps_Inter_info_all, haps_Inter_info)
}

str(n_env_FavUnfavInfo_all)
str(haps_Fav_info_all)
str(haps_Unfav_info_all)
str(haps_Inter_info_all)

png(paste(outfolder, "/boxplot_n_env_EarlyTraits.png", sep =""), width = 2700, height = 900, res = 300)
par(mar = c(3, 5.9, 0, 0) + 0.1, mgp = c(2,0.6,0))
N_env <- c(	n_env_max - n_env_FavUnfavInfo_all$n_env_neut,
			n_env_max - haps_Fav_info_all$n_env_neut,
			n_env_max - haps_Unfav_info_all$n_env_neut,
			n_env_max - haps_Inter_info_all$n_env_neut)
categ <- c(	rep("All", nrow(n_env_FavUnfavInfo_all)),
			rep("Favorable", nrow(haps_Fav_info_all)),
			rep("Unfavorable", nrow(haps_Unfav_info_all)),
			rep("Interacting", nrow(haps_Inter_info_all)))

dat <- data.frame(N_env = N_env, Category = categ)
dat$Category <- factor(dat$Category, levels = c("All", "Interacting", "Unfavorable", "Favorable"))

boxplot(N_env ~ Category, data = dat, xlab = "N sign. environments", ylab = "", col = c("gray30", "palevioletred2", "red", "blue"), cex.axis = 1, cex.lab = 1.1, horizontal = TRUE, las = 1, ylim = c(1,n_env_max))

allEarly <- c(	mean(n_env_max - haps_Fav_info_all$n_env_neut),
				mean(n_env_max - haps_Unfav_info_all$n_env_neut),
				mean(n_env_max - haps_Inter_info_all$n_env_neut),
				mean(n_env_max - n_env_FavUnfavInfo_all$n_env_neut))
mean_nEnv <- rbind(mean_nEnv, allEarly)

points(x = allEarly, y = c(4,3,2,1), bg = "gray", pch = 23)

dev.off()

allEarly <- c(	median(n_env_max - haps_Fav_info_all$n_env_neut),
				median(n_env_max - haps_Unfav_info_all$n_env_neut),
				median(n_env_max - haps_Inter_info_all$n_env_neut),
				median(n_env_max - n_env_FavUnfavInfo_all$n_env_neut))
median_nEnv <- rbind(median_nEnv, allEarly)

nrow(n_env_FavUnfavInfo_all)
nrow(haps_Fav_info_all)
nrow(haps_Unfav_info_all)
nrow(haps_Inter_info_all)
summary(n_env_max - n_env_FavUnfavInfo_all$n_env_neut)
summary(n_env_max - haps_Fav_info_all$n_env_neut)
summary(n_env_max - haps_Unfav_info_all$n_env_neut)
summary(n_env_max - haps_Inter_info_all$n_env_neut)

write.table(mean_nEnv, file = paste(outfolder, "/mean_nEnv_perTrait.csv", sep = ""), sep = ";", dec = ".", col.names = NA, row.names = TRUE, quote= FALSE)
write.table(median_nEnv, file = paste(outfolder, "/median_nEnv_perTrait.csv", sep = ""), sep = ";", dec = ".", col.names = NA, row.names = TRUE, quote= FALSE)
mean_nEnv
median_nEnv

mean_nEnv["allEarly", ] / 11
median_nEnv["allEarly", ] / 11
