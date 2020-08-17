###################################################################
###################################################################
####
#### calculate average physical and genetic window size
####
#### Manfred Mayer (Technical University of Munich, Plant Breeding)
#### manfred.mayer@tum.de
####
#### date: 03.07.2020
###################################################################
###################################################################

# general settings
options(stringsAsFactors=FALSE)
options(scipen=99)
options(warn = 1)
set.seed(212)

##############
# load phys and  genet positions of haplotype windows (generated in script 2.06)
load("Info_Lists_genet_phys.RData")
str(Info_Lists)
phys <- NULL
genet <- NULL
for(CHR in 1:10){
	phys <- rbind(phys, Info_Lists[[CHR]][["phys_InfoList"]])
	genet <- rbind(genet, Info_Lists[[CHR]][["genet_InfoList"]])
}
str(phys)
str(genet)

head(phys)
head(genet)

summary(phys[ ,"Size_bp"])
summary(genet[ ,"Size_cM"])

##############
############################################################################################################################################################################################################




