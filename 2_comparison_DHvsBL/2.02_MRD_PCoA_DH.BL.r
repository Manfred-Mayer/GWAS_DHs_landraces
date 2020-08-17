###################################################################
###################################################################
####
#### PCoA based on modified Rogers distances
#### 1) Calculate MRD
#### 2) Calculate PCoA
#### 3) Plot PCoA
####
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
library(ape)


###################################################################
###
### 1) Calculate MRD
###

# load genotypic data
load("gpDH.BL.RData")
geno <- gpDH.BL$geno
str(geno)

# filter for only polymorphic markers
# function for filtering for polymorphic markers (output = TRUE if polymorphic, FALSE if monomorphic)
polyFilter_f <- function(x) {
	  sd(x, na.rm=TRUE) != 0
	  }
markers_poly <- apply(geno, 2, polyFilter_f)
markers_poly <- colnames(geno)[which(markers_poly)]
markers_mono <- which(colnames(geno) %in% markers_poly == FALSE)
length(markers_poly)
length(markers_mono)
if(length(markers_mono) != 0){
geno <- geno[ , -markers_mono]
}
str(geno)

# function MDR, input=geno of gpData object
MRD_ind_f <- function(x)	{
		m <- ncol(x)	# number of SNPs
		d <- dist(x, method = "euclidean", diag = TRUE, upper = TRUE)
		D <- d/sqrt(4 * m)
		return(D)
		}

# calculate MRD
	D <- MRD_ind_f(geno)
	print(str(D))
	save(D, file = paste("MRD_DH.BL.RData", sep=""))

###
###################################################################


###################################################################
###
### 2) Calculate PCoA
###

D_sub <- as.matrix(D)
str(D_sub)

set.seed(4848)
# sample n_sam representative individuals from each landrace
n_sam <- 22
ID_KE <- sample(rownames(D_sub)[which(substr(rownames(D_sub), 4, 5) == "KE")], size = n_sam)
ID_PE <- sample(rownames(D_sub)[which(substr(rownames(D_sub), 4, 5) == "PE")], size = n_sam)
ID_LL <- sample(rownames(D_sub)[which(substr(rownames(D_sub), 4, 5) == "LL")], size = n_sam)
# make sure PE0075 is part of PE (PE0075 is the sequenced reference line; we used to highlight this line as well in the plot)
if("DH_PE0075" %in% ID_PE == FALSE){
	ID_PE[sample(n_sam, size = 1)] <- "DH_PE0075"
}
ex_KE <- rownames(D_sub)[which(substr(rownames(D_sub), 4, 5) == "KE")]
ex_KE <- ex_KE[which(ex_KE %in% ID_KE == FALSE)]
ex_PE <- rownames(D_sub)[which(substr(rownames(D_sub), 4, 5) == "PE")]
ex_PE <- ex_PE[which(ex_PE %in% ID_PE == FALSE)]
ex_LL <- rownames(D_sub)[which(substr(rownames(D_sub), 4, 5) == "LL")]
ex_LL <- ex_LL[which(ex_LL %in% ID_LL == FALSE)]
ID_KE
ID_PE
ID_LL
D_sub <- D_sub[which(rownames(D_sub) %in% c(ex_KE, ex_PE, ex_LL) == FALSE), which(colnames(D_sub) %in% c(ex_KE, ex_PE, ex_LL) == FALSE)]
str(D_sub)
D <- as.dist(D_sub, diag=TRUE, upper=TRUE)
print(str(D))

# PCoA
pcoaR <- pcoa(D)
class(pcoaR)
str(pcoaR)

###
###################################################################



###################################################################
###
### 3) Plot PCoA
###

# define colors
col2rgb_MM <- function(x, alpha){
  rgb(x[1],x[2],x[3], alpha = alpha)
}

temp <- c(0.6,0.6,0.9)
BL.2_edge <- col2rgb_MM(temp, alpha = 1)
BL.2_fill <- col2rgb_MM(temp, alpha = 0.25)

temp <- c(0,0,0.6)
BL.1_edge <- col2rgb_MM(temp, alpha = 1)
BL.1_fill <- col2rgb_MM(temp, alpha = 0.25)

temp <- c(0,0.85,0)
LR.2_edge <- col2rgb_MM(temp, alpha = 1)
LR.2_fill <- col2rgb_MM(temp, alpha = 0.25)

temp <- c(0,0.4,0)
LR.1_edge <- col2rgb_MM(temp, alpha = 1)
LR.1_fill <- col2rgb_MM(temp, alpha = 0.25)

	# colors
	# define colors
	# colors
	DH <- c(LR.2_edge, LR.2_edge, LR.2_edge)
	names(DH) <- c('DH_KE','DH_LL','DH_PE')
	EL <- c(BL.2_edge, BL.2_edge)
	names(EL) <- c('EL_NF','EL_nN')
	col_vec <- c(DH, EL)
	
	group_vec <- substr(attr(D, "Labels"),1,5)
	print(table(group_vec))
	
	name_vec <- attr(pcoaR$vectors, "dimnames")[[1]]
	name_vec <- gsub("DH_", "", name_vec)
	name_vec <- gsub("EL_NF_", "", name_vec)
	name_vec <- gsub("EL_nN_", "", name_vec)
	name_vec
	
	# prominent lines
	cands <- c("DK105", "EP1", "F7", "F2")
	cands
	
		# define vector with colours per family
		# fill color
		f_cols <- col_vec[group_vec]
		names(f_cols) <- name_vec
		f_cols
		
		# border color
		b_cols <- rep("black",length(f_cols))
		b_lwd <- rep(0.4, length(f_cols))
		names(b_lwd) <- name_vec
		b_lwd[cands] <- 1.3

		# vector for plotting type
		pchs <- group_vec
		pchs[pchs %in% names(EL)] <- 24
		pchs[pchs == 'DH_KE'] <- 21
		pchs[pchs == 'DH_LL'] <- 22
		pchs[pchs == 'DH_PE'] <- 23
		names(pchs) <- name_vec
		pchs[cands] <- 25		
		pchs <- as.numeric(pchs)

		# cex
		cex <- rep(0.9, length(pchs))
		cex[which(pchs == 25)] <- 1.1

		# vectors for legend
		leg <- c("KE", "LL", "PE", "Breeding lines", "Founder lines")
		l_pchs <- c(21, 22, 23, 24, 25)
		l_f_cols <- c(LR.2_edge, LR.2_edge, LR.2_edge, BL.2_edge, BL.2_edge)
		l_b_cols <- rep("black", 5)

s <- sample(seq(1, length(pchs),1))
label <- name_vec
label[which(pchs != 25)] <- NA
#f_cols[which(f_cols == "yellow")] <- "darkblue"

png(file = paste("PCoA_DH.BL_PCo12_sample", n_sam, ".png", sep=""), width=2100, height=2100, res=300)
par(mai=c(1,1,0,0), mgp = c(3.3,1,0))
plot(x=10, y=10, xlab = paste("PCo 1 (",round(pcoaR$values$Relative_eig[1] *100, digits=1),"%)", sep = ""),
	ylab = paste("PCo 2 (",round(pcoaR$values$Relative_eig[2] *100, digits=1),"%)", sep = ""), 
	las=1, cex.axis=1.2, cex.lab=1.5, cex=0.8, ylim=c(-0.215,0.215),xlim=c(-0.215,0.215))
points(pcoaR$vectors[s,c(1,2)], pch=pchs[s], col=b_cols[s], bg=f_cols[s], cex = cex[s], lwd = b_lwd[s])
legend("topright", ncol=1, pch=as.numeric(l_pchs), col=l_b_cols, pt.bg=l_f_cols, pt.lwd = c(rep(0.4, 4), 1.2), legend=leg, bty="n", cex=1.2)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
text(x = pcoaR$vectors[s,1], y = pcoaR$vectors[s,2]-0.01, labels = label[s], cex = 1, pos = 4, offset = 0.4)
dev.off()


###
###################################################################

