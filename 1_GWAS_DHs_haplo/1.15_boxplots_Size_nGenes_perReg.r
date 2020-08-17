###################################################################
###################################################################
####
#### generate boxplots for
#### the size of the identified trait associated regions
#### the number of annotated genes included in the identified trait associated regions
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

################################################
################################################
################################################
traits <- c("EV_V4", "EV_V6",
            "PH_V4", "PH_V6",
            "PH_final",
            "FF", "MF",
            "LO", "TILL"
)

regions <- NULL
for(TRAIT in traits){
  print(TRAIT)
  regions_i <- read.table(paste(TRAIT, "/QTLregs/finalQTL/finalQTLregs_", TRAIT, ".csv", sep = ""), sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)
  regions_i$trait <- rep(TRAIT, nrow(regions_i))
  print(str(regions_i))
  regions <- rbind(regions, regions_i)
}

str(regions)
head(regions)
tail(regions)

write.table(regions, file = "regions_allTraits.csv", sep = ";", dec = ".", quote = FALSE, col.names = TRUE, row.names = FALSE)


#
# plot size
#
size_kb <- regions$size_qtl_kb
summary(size_kb)

png("boxplot_regions_allTraits_size_kb.png", width=1200, height=1800, res=300)
par(mai=c(0.5,1.1,0.1,0.1), las=1, mgp=c(4,1,0))
boxplot(size_kb, xlab = "", ylab = "Size (kb)", col = "gray75", cex.lab = 1.2, cex.axis = 1.1)
title(xlab="Trait assoc. regions", mgp=c(1,1,0), cex.lab = 1.2)
dev.off()


#
# plot n-genes
#
n_genes <- regions$n_genes
summary(n_genes)
length(which(n_genes == 0))
length(n_genes)
length(which(n_genes == 0)) / length(n_genes)
length(which(n_genes == 1)) / length(n_genes)
length(which(n_genes == 2)) / length(n_genes)
length(which(n_genes == 3)) / length(n_genes)

png("boxplot_regions_allTraits_n_genes.png", width=1200, height=1800, res=300)
par(mai=c(0.5,1.1,0.1,0.1), las=1, mgp=c(3,1,0))
boxplot(n_genes, xlab = "", ylab = "Number of genes", col = "gray75", cex.lab = 1.2, cex.axis = 1.1)
title(xlab="Trait assoc. regions", mgp=c(1,1,0), cex.lab = 1.2)
dev.off()

length(which(n_genes > 10)) / length(n_genes)
length(which(n_genes > 50)) / length(n_genes)
length(which(n_genes > 100)) / length(n_genes)


