# R-scripts accompanying the article 'Discovery of beneficial haplotypes for complex traits in maize landraces' (Mayer et al., 2020)<br/>

### Reference<br/>
Mayer M, Hölker AC, González-Segovia E, Bauer E, Presterl T, Ouzunova M, Melchinger AE, Schön CC (2020) Discovery of beneficial haplotypes for complex traits in maize landraces. Nat. Commun. *issue and page not available yet*. *doi not available yet*<br/>

## System requirements<br/>

### Hardware and OS<br/>
The scripts were tested on Linux operating system but they should also be compatible with MAC and, except for the GEMMA application, also for Windows operating systems. We recommend a computer with the following specs:<br/>
RAM: 32+ GB<br/>
CPU: 4+ cores, 3.3+ GHz/core<br/>

### Software<br/>
R version 3.6.0<br/>
GEMMA 0.98.1<br/>
R packages: 'ape' version 5.3; 'asreml' version 3.0; 'ggplot2' version 3.2.0; 'plot3D' version 1.3; 'synbreed' version 0.12-9; 'VennDiagram' version 1.6.20; 'zoo' version 1.8-6<br/>

## Pipeline<br/>
### 1. Haplotype based GWAS in landrace derived DH lines<br/>
Folder: *1_GWAS_DHs_haplo*<br/>

**1.01_dataPrep.r**: Prepares genotypic data; output is a gpData object (synbreed package).<br/>

**1.02_haplotypeConstruction.r**: Constructs haplotypes to be used in GWAS; main output is a binary coded haplotype matrix.<br/>

**1.03_input_GWAS.r**: Generates all input files needed for the analysis with gemma, including an sh-script with the respective gemma commands; run this sh-script before step 1.04.<br/>

**1.04_summary_GWAS.r**: Summarizes gemma output; generates Manhattan plot, QQ plot and list of significant haplotypes for each GWAS run.<br/>

**1.05_defineRegions_GWAS.r**: Defines trait associated regions by assigning significant haplotypes to one region if they are within 1 Mb of each other and in high linkage disequilibrium (*r*<sup>2</sup> ≥ 0.8); output is a list of regions for each environment-specific GWAS run, with a specified focus haplotype for each region.<br/>

**1.06_backwardElimination_GWAS.r**: Performs backward elimination of the haplotypes from step 1.05, using a multi-locus multi-environment model; for all haplotypes passing the backward elimination, environment-specific effect estimates are calculated; the script further outputs the proportion of genotypic variance explained by the haplotypes in a simultaneous fit.<br/>

**1.07_recheckFinalHapSet.r**: In rare cases, multiple haplotypes within the same window might be significant in GWAS but only one of them was retained in the backward elimination due to collinearity; here the respective other haplotypes are tested and if significant added to the set of focus haplotypes.<br/>

**1.08_finalRegionsTable.r**: Generates final list of trait associated genomic regions, based on the final set of focus haplotypes (obtained in 1.06) and the corresponding genomic regions calculated in step 1.05.<br/>

**1.09_effectDirections_and_Stability.r**: Distinguishes between favorable, unfavorable and interacting haplotypes.<br/>

**1.10_stability_over_landraces.r**: Estimates landrace-specific effects for the set of focus haplotypes and compares effect direction and significance between landraces.<br/>

**1.11_modelWholeWindow.r**: For every 10-SNP window including a focus haplotype, this script calculates the proportion of genetic variance explained by the whole window (haplotypes coded as categorical variable) and the effect of each alternative haplotype within the respective window relative to the focus haplotype.<br/>

**1.12_bivarModel.r**: Analyzes in a bivariate model if the focus haplotypes of one trait also have a significant effect on other traits.<br/>

**1.13_plotRegions.r**: Plots the identified trait-associated regions (as defined in 1.08) for all traits.<br/>

**1.14_plotEffects_perEnvironment.r**: Plots for a given trait the position of the identified focus haplotypes as well as their effect direction and size in each environment.<br/>

**1.15_boxplots_Size_nGenes_perReg.r**: Generates boxplots for the size of identified trait- associated regions (in kb) as well as for the number of annotated genes within those regions.<br/>

**1.16_plotFavUnfavInter_nEnv.r**: Generates boxplots for the number of environments for which focus haplotypes of different categories (favorable, unfavorable, interacting, all) show significant effects.<br/>

**1.17_barplots_effectOnMultiTraits.r**: Uses results from 1.12 to generate barplots with the proportion of genotypic variance explained by haplotypes with effects on multiple traits.<br/>

### 2. Comparison of landrace derived DH lines and breeding lines<br/>
Folder: *2_comparison_DHvsBL*<br/>

**2.01_dataPrep.r**: Generates merged file with genotypic data from DH and breeding lines (BL).<br/>

**2.02_MRD_PCoA_DH.BL.r**: Calculates Modified Rogers Distances (MRD) between all lines included in 2.01; subsamples 22 DH lines per landrace; performs Principal Coordinate Analysis (PCoA); generates PCoA plot.<br/>

**2.03_haplotypeConstruction.r**: Constructs haplotypes for the combined set of DH and BL, separately for each chromosome.<br/>

**2.04_haplotype_mergeCHR.r**: Merges the haplotype files from each chromosome to one file for the whole genome.<br/>

**2.05_haplotypeCompareAll_DHvsBL.r**: Compares haplotype frequencies between DH and BL; calculates correlation of haplotype frequencies between DH and BL; generates heatmap for the haplotype frequencies; generates Venn diagram showing the number of haplotypes overlapping and private for BL and DH.<br/>

**2.06_FreqFocusHaps_inBL.r**: Calculates the frequencies of identified favorable or unfavorable focus haplotypes as well as of 500 random haplotypes in the set of 65 BL; performs Mann-Whitney test if frequencies are significantly different; generates density plots.<br/>

**2.07_CompareFocusHap_pheno_DHvsBL.r**: For a given focus haplotype, the script compares the phenotypic performance between DH lines carrying the focus haplotype and breeding lines not carrying the focus haplotype; performs permutation test; generates density plots for across environment BLUEs; generates boxplots per environment.<br/>

**2.08_CompareFocusHap_effectsPlot.r**: For a given window, the script calculates the frequencies of all haplotypes within that window for DHs and BLs; plots the effects of all alternative haplotypes relative to the focus haplotype (per environment).<br/>

**2.09_calculateAverage_cM.r**: Calculates the average physical and genetic window size.<br/>

### 3. SNP based GWAS in landrace derived DH lines<br/>
Folder: *3_GWAS_DHs_SNPs*<br/>

**3.01_input_GWAS_SNPs.r**: Generates all input files needed for the SNP-based analysis with gemma, including an sh-script with the respective gemma commands; run this sh-script before step 3.02.<br/>

**3.02_summary_GWAS_SNPs.r**: Summarizes gemma output; generates Manhattan plot, QQ plot and list of significant SNPs for each GWAS run.<br/>

**3.03_defineRegions_GWAS_SNPs.r**: Defines trait associated regions by assigning significant SNPs to one region if they are within 1 Mb of each other and in high linkage disequilibrium (*r*<sup>2</sup> ≥ 0.8); output is a list of regions for each environment-specific GWAS run, with a specified focus SNP for each region.<br/>

**3.04_backwardElimination_GWAS_SNPs.r**: Performs backward elimination of the SNPs from step 3.03, using a multi-locus multi-environment model; for all SNPs passing the backward elimination, environment-specific effect estimates are calculated; the script further outputs the proportion of genotypic variance explained by the SNPs in a simultaneous fit.<br/>

**3.05_finalRegionsTable_SNPs.r**: Generates final list of trait associated genomic regions, based on the final set of focus SNPs (obtained in 3.04) and the corresponding genomic regions calculated in step 3.03.<br/>

### 4. Calculation of diversity parameters<br/>
Folder: *4_DiversityMeasures*<br/>

**4.01_PIC_DH.BL_SNPs.r**: Calculates Polymorphism Information Content based on SNPs (PIC<sub>SNP</sub>) separately for DH lines from KE, LL and PE, for the combined set of DH lines and for the set of breeding lines.<br/>

**4.02_PIC_DH.BL_haplotypes.r**: Calculates Polymorphism Information Content based on haplotypes (PIC<sub>HAP</sub>) separately for DH lines from KE, LL and PE, for the combined set of DH lines and for the set of breeding lines.<br/>

**4.03_expHet_DH.BL_SNPs.r**: Calculates expected heterozygosity based on SNPs (*H*<sub>SNP</sub>) separately for DH lines from KE, LL and PE, for the combined set of DH lines and for the set of breeding lines.<br/>

**4.04_expHet_DH.BL_haplotypes.r**:  Calculates expected heterozygosity based on haplotypes (*H*<sub>HAP</sub>) separately for DH lines from KE, LL and PE, for the combined set of DH lines and for the set of breeding lines.<br/>

**4.05_nR_DH.BL.r**: Calculates minimum number of historical recombination events (nR) separately for DH lines from KE, LL and PE, for the combined set of DH lines and for the set of breeding lines.<br/>

**4.06_LD_DH.BL.r**: Calculates pairwise *r*<sup>2</sup> between SNPs within 1 Mb of each other as well as *r*<sup>2</sup>  decay distance for DH lines from KE, LL and PE, for the combined set of DH lines and for the set of breeding lines.<br/>
