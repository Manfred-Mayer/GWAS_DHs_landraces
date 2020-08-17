###################################################################
###################################################################
####
#### test the focus haplotypes of one trait (TRAIT1) for having a significant effect also on anothor trait (TRAIT2)
#### in a bivariate model
####
#### Manfred Mayer (Technical University of Munich, Plant Breeding)
#### manfred.mayer@tum.de
####
#### date: 25.05.2020
###################################################################
###################################################################

# general settings
options(stringsAsFactors=FALSE)
options(scipen=99)
options(warn = 1)
set.seed(212)

# load packages
library(synbreed)
library(asreml)

# graphical parameters
cex <- 0.6
cex.axis <- 1.6
cex.lab <- 1.6

# arguments
nSNPs <- 10
steps <- 10
TRAIT1 <- "PH_V4"
TRAIT2 <- "PH_final"
FDR <- "0.15"
p_thresh <- 0.01

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
TRAIT1
TRAIT2

# functions

# calculate a positive semidefinite matrix
possemiD <- function(matrix, epsilon = 0.0001){
  eigenval <- eigen(matrix)
  eigenval$values[eigenval$values < epsilon] <- epsilon
  n <- length(eigenval$values)
  Eigen <- matrix(nrow = n, ncol = n, data = 0)
  diag(Eigen) <- eigenval$values
  kmat2 <- (eigenval$vectors %*% Eigen) %*% t(eigenval$vectors)
  rownames(kmat2) <- rownames(matrix)
  colnames(kmat2) <- colnames(matrix)
  return(kmat2)
}

# pin function (ASReml) http://www.homepages.ed.ac.uk/iwhite/asreml/
pin <- function (object, transform) 
{
  pframe <- as.list(object$gammas)
  #res.pos <- which(names(pframe) == "R!variance")
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), 
                 pframe)
  X <- as.vector(attr(tvalue, "gradient"))
  #X[res.pos] <- 0
  #X[object$gammas.type == 1] <- 0
  tname <- if (length(transform) == 3) 
    transform[[2]]
  else ""
  n <- length(pframe)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  k <- 1 + (i > j)
  Vmat <- object$ai
  se <- sqrt(sum(Vmat * X[i] * X[j] * k))
  data.frame(row.names = tname, Estimate = tvalue, SE = se)
}

# input folder
infolder <- paste(TRAIT1, "/QTLregs/finalQTL", sep = "")
infolder

###################################################################
###
### 
###

# geno and map
load(paste("geno_InfoList_DH_m", nSNPs, "s", steps, ".RData", sep = ""))
str(InfoList)
str(geno)

##
## load final set of haplotypes for TRAIT1
##
# load haps of TRAIT1
load(paste(infolder, "/final_HapSet_", TRAIT1, ".RData", sep =""))
finalHapSet <- unique(info_qtls_sig_95_final$Marker)
finalHapSet
# if in biallelic case, haplotypes were added, we have to do it here also for geno and InfoList
if(any(finalHapSet %in% colnames(geno) == FALSE)){
  temp_x <- finalHapSet
  temp_x <- temp_x[which(temp_x %in% colnames(geno) == FALSE)]
  for(temp_xi in temp_x){
    temp_xix <- geno[ , which(substr(colnames(geno), 1, 13) == substr(temp_xi, 1, 13))]
    temp_xix[which(temp_xix == 0)] <- 99
    temp_xix[which(temp_xix == 2)] <- 0
    temp_xix[which(temp_xix == 99)] <- 2
    geno <- cbind(geno, temp_xix)
    colnames(geno)[ncol(geno)] <- temp_xi
    
    InfoList_xi <- InfoList[which(substr(rownames(InfoList), 1, 13) == substr(temp_xi, 1, 13)), ]
    InfoList <- rbind(InfoList, InfoList_xi)
    rownames(InfoList)[nrow(InfoList)] <- temp_xi
    print(tail(InfoList, n = 10))
  }
}

dir.create(paste("TraitCorBivar", sep = ""))
dir.create(paste("TraitCorBivar/", TRAIT1, ".", TRAIT2, sep = ""))
outfolder <- paste("TraitCorBivar/", TRAIT1, ".", TRAIT2, sep = "")
outfolder

##
## load input data for bivariate model
## generate input objects for bivariate model
##
#
# phenotypic data
#
load("pheno_perEnv_list_DHs.RData")
pheno_list <- pheno_perEnv_list
str(pheno_list)
pheno_TRAIT1 <- pheno_list[[TRAIT1]]
str(pheno_TRAIT1)
pheno_TRAIT2 <- pheno_list[[TRAIT2]]
str(pheno_TRAIT2)
# filter for DHs with genotypic info
pheno_TRAIT1 <- pheno_TRAIT1[which(pheno_TRAIT1$Geno %in% rownames(geno)), ]
str(pheno_TRAIT1)
head(pheno_TRAIT1)
pheno_TRAIT2 <- pheno_TRAIT2[which(pheno_TRAIT2$Geno %in% rownames(geno)), ]
str(pheno_TRAIT2)
head(pheno_TRAIT2)
# stack environments
geno_dummy <- unique(c(pheno_TRAIT1$Geno, pheno_TRAIT2$Geno))
env_dummy <- unique(c(colnames(pheno_TRAIT1)[-1], colnames(pheno_TRAIT2)[-1]))
pheno_dummy <- matrix(NA, nrow = length(geno_dummy) * length(env_dummy), ncol = 2)
colnames(pheno_dummy) <- c(TRAIT1, TRAIT2)
i <- 0
for(env_i in env_dummy){
  print(env_i)
  for(geno_i in geno_dummy){
    i <- i + 1
    x1 <- pheno_TRAIT1[which(pheno_TRAIT1$Geno == geno_i), env_i]
    x2 <- pheno_TRAIT2[which(pheno_TRAIT2$Geno == geno_i), env_i]
    if(is.null(x1)){
      x1 <- NA
    }
    if(is.null(x2)){
      x2 <- NA
    }
    pheno_dummy[i, TRAIT1] <- x1
    pheno_dummy[i, TRAIT2] <- x2
  }
}
# landrace covariates
LR_cov <- rep(NA, length(geno_dummy))
LR_cov[which(substr(geno_dummy, 1, 5) == "DH_KE")] <- "KE"
LR_cov[which(substr(geno_dummy, 1, 5) == "DH_LL")] <- "LL"
LR_cov[which(substr(geno_dummy, 1, 5) == "DH_PE")] <- "PE"
print(str(LR_cov))
print(table(LR_cov))

pheno <- data.frame(GENOTYPE = rep(geno_dummy, length(env_dummy)),
                    LR = rep(LR_cov, length(env_dummy)),
                    ENV = rep(env_dummy, rep(length(geno_dummy), length(env_dummy))),
                    pheno_dummy
)
str(pheno)
head(pheno)

#
# kinship matrix
#
# filter kinship for overlapping individuals
load("Kin_AB_DHs_SNPs.RData")
Kinship <- Kinship[geno_dummy, geno_dummy]
print(all.equal(rownames(Kinship), geno_dummy))
print(all.equal(colnames(Kinship), geno_dummy))
str(Kinship)
G <- possemiD(Kinship, epsilon = 0.0001)
G_inv <- solve(G)

#
# establish bivariate model without markers as fixed effects (null_model)
#
str(pheno)
pheno$GENOTYPE <- factor(pheno$GENOTYPE)
pheno$LR <- factor(pheno$LR)
pheno$ENV <- factor(pheno$ENV)
colnames(pheno)[4] <- "y1"
colnames(pheno)[5] <- "y2"
str(pheno)
head(pheno)
pheno <- pheno[which((is.na(pheno[,"y1"]) & is.na(pheno[,"y2"])) == FALSE),]
str(pheno)
head(pheno)

# Null model
model_0 <- asreml(fixed = cbind(y1,y2) ~ trait + trait:LR + trait:ENV, 
                  random   = ~giv(GENOTYPE):us(trait),
                  ginverse = list(GENOTYPE = G_inv),
                  rcov=~units:us(trait),
                  maxiter = 50, data = pheno,
                  workspace=40e08,
                  start.values=TRUE)

iv.null <- model_0$gammas.table

model_0 <- asreml(fixed = cbind(y1,y2) ~ trait + trait:LR + trait:ENV, 
                  random   = ~giv(GENOTYPE):us(trait),
                  ginverse = list(GENOTYPE = G_inv),
                  rcov=~units:us(trait),
                  maxiter = 50, data = pheno,
                  workspace=40e08)

iv.null$Value  <- summary(model_0)$varcomp$component

model_0 <- asreml(fixed = cbind(y1,y2) ~ trait + trait:LR + trait:ENV, 
                  random   = ~giv(GENOTYPE):us(trait),
                  ginverse = list(GENOTYPE = G_inv),
                  rcov=~units:us(trait),
                  maxiter = 50, data = pheno,
                  workspace=40e08,
                  G.param=iv.null, 
                  R.param=iv.null)

varcomp_G <-  summary(model_0)$varcomp$component[1:3]
varcomp_R <-  summary(model_0)$varcomp$component[5:7]
varcomp_G
varcomp_R
# calculate genomic correlation
genomic_cor_y1y2 <- pin(model_0, genomicCor ~ V2/sqrt(V1*V3))
genomic_cor_y1y2

#
# now, include markers as fixed effects
#
finalHapSet <- sort(finalHapSet)
finalHapSet 

overall_p_main <- NULL
Gcomp <- NULL
Rcomp <- NULL
genomicCor_y1y2_est <- NULL
genomicCor_y1y2_se <- NULL
info_CI <- matrix(NA, nrow = length(finalHapSet), ncol = 6)
rownames(info_CI) <- finalHapSet
colnames(info_CI) <- c("eff_y1", "lowCI_y1", "upCI_y1",
                       "eff_y2", "lowCI_y2", "upCI_y2"
)
ori <- NULL
info_list_CI <- list()
for(mk_i in finalHapSet){
  print(mk_i)
  
  ori <- c(ori, TRAIT1)

  haplo <- geno[as.character(pheno$GENOTYPE), mk_i]
  pheno_i <- pheno
  pheno_i$haplo <- haplo
  
  model_haplo <- asreml(fixed = cbind(y1,y2) ~ trait + trait:ENV + trait:LR + trait:haplo, 
                        random   = ~giv(GENOTYPE):us(trait, init = varcomp_G),
                        ginverse = list(GENOTYPE = G_inv),
                        rcov=~units:us(trait, init = varcomp_R),
                        maxiter = 50, data = pheno_i,
                        workspace=40e08)
  
  wald.obj <- wald(model_haplo, denDF = "default", workspace=40e08)
  print(wald.obj$Wald)
  overall_p_main <- c(overall_p_main, wald.obj$Wald["trait:haplo", "Pr"])
  Gcomp <- rbind(Gcomp, summary(model_haplo)$varcomp$component[1:3])
  Rcomp <- rbind(Rcomp, summary(model_haplo)$varcomp$component[5:7])
  genomicCor_y1y2_i <- pin(model_haplo, genomicCor ~ V2/sqrt(V1*V3))
  genomicCor_y1y2_est <- c(genomicCor_y1y2_est, genomicCor_y1y2_i$Estimate[1])
  genomicCor_y1y2_se <- c(genomicCor_y1y2_se, genomicCor_y1y2_i$SE[1])
  
  eff <- model_haplo$coefficients$fixed
  SEeff <- sqrt(model_haplo$vcoeff$fixed)
  lowerCI95 <- eff - qnorm(0.975, mean = 0, sd = 1) * SEeff
  upperCI95 <- eff + qnorm(0.975, mean = 0, sd = 1) * SEeff
  info_CI[mk_i, ] <- c( eff["trait_y1:haplo"], lowerCI95["trait_y1:haplo"], upperCI95["trait_y1:haplo"],
                        eff["trait_y2:haplo"], lowerCI95["trait_y2:haplo"], upperCI95["trait_y2:haplo"]
  )
  print(info_CI[mk_i, ])
  info_list_CI_i <- list(eff, SEeff, lowerCI95, upperCI95)
  info_list_CI[[mk_i]] <- info_list_CI_i
}

colnames(Gcomp) <- c(paste0("VARg_", TRAIT1), paste0("COVg"), paste0("VARg_", TRAIT2))
colnames(Rcomp) <- c(paste0("VARr_", TRAIT1), paste0("COVr"), paste0("VARr_", TRAIT2))
VAR_y1_expl <- (varcomp_G[1] - Gcomp[,1]) / varcomp_G[1]
COV_y12_expl <- (varcomp_G[2] - Gcomp[,2]) / varcomp_G[2]
VAR_y2_expl <- (varcomp_G[3] - Gcomp[,3]) / varcomp_G[3]
VARCOV_expl <- cbind(VAR_y1_expl, COV_y12_expl, VAR_y2_expl)
colnames(VARCOV_expl) <- c(paste0("VARexpl_", TRAIT1), paste0("COVexpl"), paste0("VARexpl_", TRAIT2))

info <- data.frame( haplotype = finalHapSet,
                    Q_trait = ori,
                    p_main = overall_p_main,
                    info_CI,
                    Gcomp,
                    Rcomp,
                    VARCOV_expl,
                    genomicCor_y1y2_est,
                    genomicCor_y1y2_se
)
sig_y1 <- rep(0, nrow(info))
sig_y2 <- rep(0, nrow(info))
sig_y1[which(((info$lowCI_y1 > 0) & (info$upCI_y1 > 0)) | ((info$lowCI_y1 < 0) & (info$upCI_y1 < 0)))] <- 1
sig_y2[which(((info$lowCI_y2 > 0) & (info$upCI_y2 > 0)) | ((info$lowCI_y2 < 0) & (info$upCI_y2 < 0)))] <- 1

info$sig_y1 <- sig_y1
info$sig_y2 <- sig_y2
info

write.table(info, file = paste(outfolder, "/info_cor_", TRAIT1, "_", TRAIT2, "_main.csv", sep =""), sep = ";", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)
save(info_list_CI, genomic_cor_y1y2, varcomp_G, varcomp_R, file = paste(outfolder, "/info_list_CI_", TRAIT1, "_", TRAIT2, "_main.RData", sep =""))

###
###################################################################
