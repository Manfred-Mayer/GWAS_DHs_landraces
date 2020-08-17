###################################################################
###################################################################
####
#### calculate r2 and r2 decay distances according to Hill and Robertson (1968) and Hill and Weir (1988), respectively
#### for DHs of KE, LL, PE separately as well as overall
#### for BLs
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

# load packages
library(synbreed)

# arguments
nSNPs <- 10
steps <- 10
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
group

# function for filtering polymorphic markers (output = TRUE if polymorphic, FALSE if monomorphic)
polyFilter_f <- function(x) {
	  sd(x, na.rm=TRUE) != 0
	  }

# function for non-linear regression
smooth.fit <- function(overallDist, overallr2, n) {
        nonlinearoverall <- nls(overallr2 ~ ((10 + p * overallDist))/((2 +
            p * overallDist) * (11 + p * overallDist)) * (1 +
            ((3 + p * overallDist) * (12 + 12 * p + p^2 * overallDist^2))/(n *
            (2 + p * overallDist) * (11 + p * overallDist))),
            start = list(p = 0.001), contro=list(maxiter=1000))
        p <- coef(nonlinearoverall)
		return(p)	        
		} 

# function for calculating the distance for which the fitted curve crosses the r_threshold (here 0.2)
fitcurve2 <- function(x,p,n) {
   -1*0.2 + ((10 + p*x)) / ((2+p*x) * (11 + p*x) ) *
     ( 1 + ( (3+ p*x) * (12 + 12 * p + p^2*x^2)) / ( n*(2+p*x) * (11 + p*x)))	
	 }

###################################################################
###

# output folder
dir.create(paste("diversityParams", sep = ""))
dir.create(paste("diversityParams/LD", sep = ""))
outfolder <- paste("diversityParams/LD", sep = "")
outfolder

# load genotypic data of DH and BL lines
load("gpDH.BL.RData")
str(gpDH.BL)

# filter for respective germplasm group
geno <- gpDH.BL$geno
geno <- geno[which(substr(rownames(geno), 1, nchar(group)) %in% group), ]
print(rownames(geno))
print(dim(geno))
geno[geno==2] <- 1

n <- nrow(geno)
print(paste("n =", n))

map <- gpDH.BL$map

rm(gpDH.BL)
gc()

pos <- map$pos
names(pos) <- rownames(map)

LDall_1Mb <- NULL
ALLdist_all <- NULL
mean_r2_all <- NULL

for (CHR in 10:1)	
	{
	LD_chr_1Mb <- NULL
	Alldist_chr_parts <- mean_r2_chr_parts <- NULL
	ALLdist <- mean_r2_chr <- NULL
	print(paste("Chr",CHR,sep=""))
	# gametes in rows, markers in columns
	hap <- geno[, map$chr == CHR]
	# Only polymorphic markers
	# per LR
	n_all <- ncol(hap)
	hap <- hap[ ,apply(hap, 2, polyFilter_f)]
	n_poly <- ncol(hap)
	print(paste(n_all - n_poly, "monomorphic markers removed", sep=" "))
	print(ncol(hap))
	# because of memory limitations of R, subdividing of the corr-matrix is required
	rx <- ncol(hap) # number of markers
	# devide matrix into 10 parts:
	
														#################################
														#		#		#		#		#
														# part1	# part2	# part3	# part4	#
														#		#		#		#		#
														#################################
														#		#		#		#		#
														# 		# part5	# part6	# part7	#
														#		#		#		#		#
														#################################
														#		#		#		#		#
														#		# 		# part8	# part9	#
														#		#		#		#		#
														#################################
														#		#		#		#		#
														# 		# 		# 		# part10#
														#		#		#		#		#
														#################################
						
	# create vector for break points between the parts
	breaks <- cbind(c(1, 1, 1, 1, (floor(rx*(1/4))+1), (floor(rx*(1/4))+1), (floor(rx*(1/4))+1), (floor(rx*(2/4))+1), (floor(rx*(2/4))+1), (floor(rx*(3/4))+1)),
					c(floor(rx*(1/4)), floor(rx*(1/4)), floor(rx*(1/4)), floor(rx*(1/4)), floor(rx*(2/4)), floor(rx*(2/4)), floor(rx*(2/4)), floor(rx*(3/4)), floor(rx*(3/4)), floor(rx*(4/4))),
					c(1,(floor(rx*(1/4)+1)), (floor(rx*(2/4)+1)), (floor(rx*(3/4)+1)), (floor(rx*(1/4)+1)), (floor(rx*(2/4)+1)), (floor(rx*(3/4)+1)), (floor(rx*(2/4)+1)), (floor(rx*(3/4)+1)), (floor(rx*(3/4)+1))),
					c(floor(rx*(1/4)), floor(rx*(2/4)), floor(rx*(3/4)), floor(rx*(4/4)), floor(rx*(2/4)), floor(rx*(3/4)), floor(rx*(4/4)), floor(rx*(3/4)), floor(rx*(4/4)), floor(rx*(4/4))))
	rownames(breaks) <- paste("part", c(1:10), sep="")
	colnames(breaks) <- c("x1", "x2", "y1", "y2")
	print(breaks)
	# calculation for each part seperately (10 parts of the matrix)
	for (part in 1:10)	{
	if (part %in% c(1,5,8,10)) {
		cols <- cbind(colnames(hap)[sequence(1:(breaks[part, 2]-breaks[part, 1]+1))+(breaks[part, 1]-1)],colnames(hap)[rep(breaks[part, 3]:breaks[part, 4], 1:(breaks[part, 2]-breaks[part, 1]+1))])
		} else {
			if((breaks[part, 2]-breaks[part, 1]) == (breaks[part, 4]-breaks[part, 3])) {
			cols <- cbind(colnames(hap)[rep(breaks[part, 1]:breaks[part, 2],(breaks[part, 2]-breaks[part, 1]+1))],colnames(hap)[rep(breaks[part, 3]:breaks[part, 4], each=(breaks[part, 2]-breaks[part, 1]+1))])
			} else {
			if((breaks[part, 2]-breaks[part, 1]) < (breaks[part, 4]-breaks[part, 3])) {
				cols <- cbind(colnames(hap)[rep(breaks[part, 1]:breaks[part, 2],(breaks[part, 2]-breaks[part, 1]+1))],colnames(hap)[rep(breaks[part, 3]:(breaks[part, 4]-1), each=(breaks[part, 2]-breaks[part, 1]+1))])
				cols <- rbind(cols,cbind(colnames(hap)[breaks[part, 1]:breaks[part, 2]], colnames(hap)[rep(breaks[part, 4],(breaks[part, 2]-breaks[part, 1]+1))]))
				} else {
				cols <- cbind(colnames(hap)[rep(breaks[part, 1]:breaks[part, 2],(breaks[part, 2]-breaks[part, 1]))],colnames(hap)[rep(breaks[part, 3]:(breaks[part, 4]), each=(breaks[part, 2]-breaks[part, 1]+1))])
				}
			}			
		}				
	# calculate r^2
	m <- x <- y <- NULL
	if (part %in% c(1,5,8,10)) {
	# calculate distances between markers
	distance <- pos[cols[,2]]-pos[cols[,1]]
	m <- hap[,colnames(hap)[breaks[paste("part",part,sep=""),1]:breaks[paste("part",part,sep=""),2]]]
	print(paste("part_",part," m  from ",colnames(m)[1], " to ", colnames(m)[ncol(m)], sep=""))
	print(paste("part_",part," m  from ",colnames(m)[1], " to ", colnames(m)[ncol(m)], sep=""))
	corrs <- cor(m, use="pairwise.complete.obs")
	corrs1<- as.numeric(corrs[upper.tri(corrs,diag=TRUE)])
	corrs2<- corrs1^2
	# merge distance and corresponding r^2 and remove entries for same SNPs
	corrs3 <- cbind(distance,corrs1,corrs2)
	corrs3 <- corrs3[which(corrs3[,1]!=0),]
	cols <- cols[which(cols[,1]!=cols[,2]),]
	} else {
	cols <- cols[which(cols[,1]!=cols[,2]),]
	distance <- pos[cols[,2]]-pos[cols[,1]]
	x <- hap[,colnames(hap)[breaks[paste("part",part,sep=""),1]:breaks[paste("part",part,sep=""),2]]]
	y <- hap[,colnames(hap)[breaks[paste("part",part,sep=""),3]:breaks[paste("part",part,sep=""),4]]]
	print(paste("part_",part," x  from ",colnames(x)[1], " to ", colnames(x)[ncol(x)], sep=""))
	print(paste("part_",part," y  from ",colnames(y)[1], " to ", colnames(y)[ncol(y)], sep=""))
	corrs <- cor(x=x, y=y, use="pairwise.complete.obs")
	corrs1<- as.numeric(corrs)
	corrs2<- corrs1^2
	corrs3 <- cbind(distance,corrs1,corrs2)
	}
	# create data.frame with SNP_names and corresponding distance and r^2 value
	SNP_A <- cols[,1]
	SNP_B <- cols[,2]
	distance <- corrs3[,1]
	r <- corrs3[,2]
	r2 <- corrs3[,3]
	LD <- data.frame(SNP_A, SNP_B, distance, r, r2, stringsAsFactors=FALSE)
	LD_chr_1Mb <- rbind(LD_chr_1Mb, LD[LD$distance <= 1000000, ])
	Alldist_chr_parts <- c(Alldist_chr_parts, as.numeric(nrow(LD)))
	mean_r2_chr_parts <- c(mean_r2_chr_parts, mean(LD$r2))
	rm(cols, distance, corrs, corrs1, corrs2, corrs3, SNP_A, SNP_B, r, r2)		
	rm(LD)
	}	
	save(LD_chr_1Mb, file=paste(outfolder, "/LD_NoMAF_Chr",CHR,"_", group, ".RData",sep=""))
	ALLdist <- sum(Alldist_chr_parts)
	mean_r2_chr <- sum((mean_r2_chr_parts*Alldist_chr_parts)/sum(Alldist_chr_parts))
	print(ALLdist)
	print(mean_r2_chr)
	ALLdist_all <- c(ALLdist_all, ALLdist)
	mean_r2_all <- c(mean_r2_all, mean_r2_chr)
	LDall_1Mb <- rbind(LDall_1Mb, LD_chr_1Mb)

	LD_chr_1Mb <- LD_chr_1Mb[, c(1,2,3,5)]
	t <- try(p <- smooth.fit(LD_chr_1Mb$distance, LD_chr_1Mb$r2, n))
	if (inherits(t, "try-error")){
	p <- NA
	}
	print(p)
	if (is.na(p)){
	r2decay <- NA
	} else {
	r2decay <- try(uniroot(f=fitcurve2,n=n, interval=c(0, 10000000000), p=p))
	if (inherits(r2decay, "try-error")){
	r2decay <- NA
	}
	}
	print(r2decay)
	r2.decay.kb <- r2decay[[1]]/1000
	# create summary output
	n_SNPs <- nrow(LD_chr_1Mb)
	min.r2 <- min(LD_chr_1Mb$r2)
	max.r2 <- max(LD_chr_1Mb$r2)
	mean.r2 <- mean(LD_chr_1Mb$r2)
	median.r2 <- median(LD_chr_1Mb$r2)
	overview <- cbind(CHR, n_SNPs, min.r2, max.r2, mean.r2, median.r2, p, r2.decay.kb)
	write.table(overview, file=paste(outfolder, "/Overview_LD_nlr_1Mb_Chr",CHR,"_", group, ".txt", sep=""), row.names=FALSE)
	rm(p, n_SNPs, min.r2, max.r2, mean.r2, median.r2, overview, r2.decay.kb)
}

# for whole genome	
	# calculate decay distance
	LDall_1Mb <- LDall_1Mb[, c(1,2,3,5)]
	t <- try(p <- smooth.fit(LDall_1Mb$distance, LDall_1Mb$r2, n))
	if (inherits(t, "try-error")){
	p <- NA
	}
	print(p)
	if (is.na(p)){
	r2decay <- NA
	} else {
	r2decay <- try(uniroot(f=fitcurve2,n=n, interval=c(0, 10000000000), p=p))
	if (inherits(r2decay, "try-error")){
	r2decay <- NA
	}
	}
	print(r2decay)
	r2.decay.kb <- r2decay[[1]]/1000
	# create summary output
	n_SNPs <- nrow(LDall_1Mb)
	min.r2 <- min(LDall_1Mb$r2)
	max.r2 <- max(LDall_1Mb$r2)
	mean.r2 <- mean(LDall_1Mb$r2)
	median.r2 <- median(LDall_1Mb$r2)
	overview <- cbind(n_SNPs, min.r2, max.r2, mean.r2, median.r2, p, r2.decay.kb)
	write.table(overview, file=paste(outfolder, "/Overview_LD_nlr_wholeG_1Mb_", group, ".txt", sep=""), row.names=FALSE)
	rm(p, n_SNPs, min.r2, max.r2, mean.r2, median.r2, overview, r2.decay.kb)
	
	mean_r2_chr_vec <- mean_r2_all
	print(mean_r2_all)
	mean_r2_all <- mean_r2_all*ALLdist_all
	print(mean_r2_all)
	mean_r2_all <- mean_r2_all/(sum(ALLdist_all))
	print(mean_r2_all)
	mean_r2_all <- sum(mean_r2_all)
	print(mean_r2_all)
	mean_r2_chr_vec <- c(mean_r2_chr_vec, mean_r2_all)
	names(mean_r2_chr_vec) <- c(paste("Chr", 1:10, sep=""), "wholeG")
	write.csv(mean_r2_chr_vec, file= paste(outfolder, "/mean_r2_perChr_wholeG_allPairs_", group, ".csv", sep=""))


###
###################################################################
