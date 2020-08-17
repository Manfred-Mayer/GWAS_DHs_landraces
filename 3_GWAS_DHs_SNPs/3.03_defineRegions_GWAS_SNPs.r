###################################################################
###################################################################
####
#### for each environment:
#### based on the significant SNPs identified in GWAS,
#### define trait associated regions, grouping SNPs together which
#### - are within 1Mb of each other
#### - have an r2 >= 0.8
####
#### Manfred Mayer (Technical University of Munich, Plant Breeding)
#### manfred.mayer@tum.de
####
#### date: 29.06.2020
###################################################################
###################################################################

# general settings
options(stringsAsFactors=FALSE)
options(scipen=99)
options(warn = 1)
set.seed(212)

# arguments
TRAIT <- "PH_final"

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
TRAIT

# thresholds for defining regions
dist_bp <- 1000000
r2_thresh <- 0.8 # this was set to be equal to the haplotype based approach, however, for SNP-based analyses you might consider lowering it

###################################################################
###
### 
###

# input folder
infolder <- paste("GWAS_SNPs/", TRAIT, "/summary_GWAS", sep ="")
infolder

# output folder
dir.create(paste("GWAS_SNPs/", TRAIT, "/QTLregs", sep = ""))
outfolder <- paste("GWAS_SNPs/", TRAIT, "/QTLregs", sep = "")
outfolder

# load annotated gene information for B73v4
# it's a dataframe with columns: name, chr, start_bp, end_bp, strand
# and one row per gene
load(paste("gene_pos_Zea_mays.B73_RefGen_v4.43_FGS.RData", sep = ""))
str(gene_pos)
# consider only genes mapped on one of the 10 chromosomes
gene_pos <- gene_pos[which(gene_pos$chr %in% c(1,2,3,4,5,6,7,8,9,10)), ]
str(gene_pos)
# also include regions ...kb upstream of genes
extGenes <- 5
for(i in 1: nrow(gene_pos)){
	if(gene_pos[i, 5] == "+"){
		gene_pos[i, 3] <- gene_pos[i, 3] - (extGenes * 1000)
	}
	if(gene_pos[i, 5] == "-"){
		gene_pos[i, 4] <- gene_pos[i, 4] + (extGenes * 1000)
	}
}
gene_pos$chr <- as.numeric(gene_pos$chr)
names(gene_pos)[1] <- "Gene_ID"
str(gene_pos)
head(gene_pos)
	
# load data
# geno and map
load(paste("geno_InfoList_DH_SNPs.RData", sep = ""))
str(InfoList)
str(geno)

# environments
Envs <- read.table(paste0("GWAS_SNPs/", TRAIT, "/", list.files(path = paste0("GWAS_SNPs/", TRAIT), pattern = ".phenotypes_order.txt")), stringsAsFactors = FALSE)
Envs <- Envs[,1]
Envs

# run for each environment (and for the across-environment-BLUEs)
for (Env_i in 1:length(Envs)){
ENV <- Envs[Env_i]
print(ENV)

		# define exactly the geno which was used for this gwas run
		Y_temp <- read.table(paste0("GWAS_SNPs/", TRAIT, "/", list.files(path = paste0("GWAS_SNPs/", TRAIT), pattern = ".pheno.txt")), stringsAsFactors = FALSE)
		colnames(Y_temp) <- Envs
		Y_temp <- Y_temp[ , ENV]
		GENOTYPE_temp <- read.table(paste0("GWAS_SNPs/", TRAIT, "/", list.files(path = paste0("GWAS_SNPs/", TRAIT), pattern = ".geno_order.txt")), stringsAsFactors = FALSE)
		GENOTYPE_temp <- GENOTYPE_temp[,1]
		GENOTYPE_temp <- GENOTYPE_temp[which(is.na(Y_temp) == FALSE)]
		geno_gwas <- geno[GENOTYPE_temp, ]

	# load significant haplotypes
	if(file.exists(paste(infolder, "/Hits_", TRAIT, "_", ENV, "_FDR0.15.csv", sep =""))){
	sigList <- read.table(paste(infolder, "/Hits_", TRAIT, "_", ENV, "_FDR0.15.csv", sep =""), sep = ";", header = TRUE, row.names = 1)
	print(str(sigList))
	
	if(is.na(sigList[1,1]) == FALSE){
		
		# define for each haplotype (+- ld based distance)
		reg_list_merge <- NULL
		gene_list_merge <- list()
		
		# 0) per chromosome
		regions <- NULL
		for(CHR in 1:10){
			print(paste("CHR",CHR))
			sigList_chr <- sigList[which(sigList$Chr == CHR), ]

			if(nrow(sigList_chr) > 0){
				  # 1) define regions around every haplotype
				  for(m_i in rownames(sigList_chr)){
					geno_focus <- geno_gwas[ , m_i]
					# consider significant haplotypes with maximum distance of "dist_bp" to haplotype m_i
					focus_start <- sigList_chr[m_i, 2]
					focus_end <- sigList_chr[m_i, 3]
					right_of <- sigList_chr[, 2] - focus_end
					left_of <- focus_start - sigList_chr[, 3]
					withinDistMarkers <- rownames(sigList_chr)[which((abs(right_of) <= dist_bp) | abs(left_of) <= dist_bp)]
					if(length(withinDistMarkers) == 1){
						LD_marker <- m_i						
					} else {
						geno_gwas_sig <- geno_gwas[ , withinDistMarkers]
						# calculate r2
						R <- cor(x = geno_gwas_sig, y = geno_focus,  use = "pairwise.complete.obs")
						R2 <- R^2
						LD_marker <- rownames(R2)[which(R2 >= r2_thresh)]
					}	
						# for all haplotypes which are in LD with haplotype m_i, consider max and min position
						Start_bp <- min(sigList_chr[LD_marker, 2])
						End_bp <- max(sigList_chr[LD_marker, 3])
						reg <- data.frame(	CHR,
											Start_bp = Start_bp,
											End_bp = End_bp,
											Marker = m_i,
											allele_freq = sigList_chr[m_i, "allele_freq"],
											pvalue = sigList_chr[m_i, "pvalue"],
											effect = sigList_chr[m_i, "effect"])
						regions <- rbind(regions, reg)
					}
				}
			}
		# regions is now data frame which states for each significant haplotype m_i the start and end point of the region with all significant haplotypes in LD with m_i
		
		###
		### merge overlapping regions
		### 
		
		list_all <- regions
		list_all <- list_all[order(list_all[,1], list_all[,2], list_all[,3]), ]
		cons_list <- NULL
		for(CHR in 1:10){
			list_all_chr <- list_all[which(list_all[,1] == CHR), ]
			cons_list_chr <- NULL
			if(length(list_all_chr) > 0){
				repeat{
					if(nrow(list_all_chr) == 0){
						break
					}
					if(nrow(list_all_chr) == 1){
						temp <- list_all_chr
						cons_list_chr <- rbind(cons_list_chr, temp)
						rm(temp)				
						break
					}
					i <- 1
					temp <- list_all_chr[i,]
					repeat{
						i <- i + 1
						if(temp[1,3] >= list_all_chr[i,2]){
							temp[1,3] <- max(c(temp[1,3], list_all_chr[i,3]))
						} else {
							i <- i - 1
							break
						}
						if(i == nrow(list_all_chr)){
							break
						}					
					}
					cons_list_chr <- rbind(cons_list_chr, temp)
					list_all_chr <- list_all_chr[-c(1:i), ]
					rm(temp, i)
				}
			}
			cons_list <- rbind(cons_list, cons_list_chr)
		}
		cons_list <- cons_list[, 1:3]
		print(cons_list)
		
		###
		### for merged regions, identify mostSigMarker (the haplotype with lowest p-value in that region), n_SNPs_reg, ...
		###
		for(reg_i in 1:nrow(cons_list)){
							sigList_reg <- sigList[which(sigList$Chr == cons_list[reg_i, 1]), ]
							sigList_reg <- sigList_reg[which(sigList_reg$Start_bp >= cons_list[reg_i, 2]), ]
							sigList_reg <- sigList_reg[which(sigList_reg$End_bp <= cons_list[reg_i, 3]), ]
							#
							p_values <- sort(sigList_reg$pvalue)
							p_values <- unique(p_values)
							# 1st the one with the lowest p-value (if multiple, take the one in the middle)
							wind1 <- sigList_reg[which(sigList_reg$pvalue == p_values[1]), ]
							if(nrow(wind1) > 1){
								wind1 <- wind1[round(mean(1:nrow(wind1))), ]
							}
							mostSigMarker <- rownames(wind1)
							Chr <- wind1$Chr[1]
							start_bp <- cons_list[reg_i, 2]
							end_bp <- cons_list[reg_i, 3]
							size_kb <- (end_bp - start_bp + 1) / 1000
							allele_freq <- wind1$allele_freq[1]
							eff <- wind1$effect[1]
							eff_prop <- wind1$effect_prop_sd[1]
							min_p <- wind1$pvalue[1]
							# how many annotated genes lie within that region
							geneIDs <- gene_pos$Gene_ID[which((gene_pos$chr == Chr) & (gene_pos$start_bp <= end_bp) & (gene_pos$end_bp >= start_bp))]
							if(length(geneIDs) > 0){
								gene_list_merge[[length(gene_list_merge)+1]] <- geneIDs
							} else {
								gene_list_merge[[length(gene_list_merge)+1]] <- NA
							}
							n_genes <- length(geneIDs)
							micro_reg <- data.frame(	Chr,
														start_bp,
														end_bp,
														size_kb,
														mostSigMarker,
														allele_freq,
														min_p,
														eff,
														eff_prop,
														n_genes
														)
							reg_list_merge <- rbind(reg_list_merge, micro_reg)	
		}
		
		# list of genes within each of the identified regions
		gene_list_merge <- gene_list_merge[order(reg_list_merge[,1], reg_list_merge[,2], reg_list_merge[,3])]
		# list of regions
		colnames(reg_list_merge) <- c("Chr", "Start_bp", "End_bp", "Size_kb", "mostSigMarker", "allele_freq", "pvalue_log10", "eff_abs", "eff_prop", "n_genes")
		reg_list_merge <- reg_list_merge[order(reg_list_merge[,1], reg_list_merge[,2], reg_list_merge[,3]), ]
		reg_list_merge[ , "pvalue_log10"] <- -log10(reg_list_merge[ , "pvalue_log10"])
		rownames(reg_list_merge) <- paste(TRAIT, ".", ENV, "_", 1:nrow(reg_list_merge), sep = "")
		print(reg_list_merge)		
		names(gene_list_merge) <- rownames(reg_list_merge)
		# outputs
		save(gene_list_merge, file = paste(outfolder, "/gene_list_merge_", TRAIT, ".", ENV, "_FDR0.15.RData", sep = ""))
		save(reg_list_merge, file = paste(outfolder, "/reg_list_merge_", TRAIT, ".", ENV, "_FDR0.15.RData", sep = ""))
		write.table(reg_list_merge, file = paste(outfolder, "/reg_list_merge_", TRAIT, ".", ENV, "_FDR0.15.txt", sep = ""), row.names=FALSE, quote = FALSE)

}

}

}

###
###################################################################
