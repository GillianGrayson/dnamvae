rm(list=ls())


if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DMRcate")
BiocManager::install("methylumi")
BiocManager::install("ChAMP")
BiocManager::install("minfi")
BiocManager::install("minfiData")
BiocManager::install("wateRmelon")
BiocManager::install("shinyMethyl")
BiocManager::install("FlowSorted.Blood.EPIC")
BiocManager::install("FlowSorted.Blood.450k")
BiocManager::install("FlowSorted.DLPFC.450k")

library(ChAMP)
library("xlsx")
library(openxlsx)
library(minfi)
library(minfiData)
library(shinyMethyl)
library(wateRmelon)
library(FlowSorted.Blood.EPIC)
library(FlowSorted.Blood.450k)
library(FlowSorted.DLPFC.450k)

path <- "E:/YandexDisk/Work/pydnameth/datasets/GPL13534/GSE87648/raw/idat"
chip_type = "450K"
detPcut = 0.01

setwd(path)

myLoad = champ.load(directory = path,
           arraytype = chip_type, # Choose microarray type is "450K" or "EPIC".(default = "450K")
           method = "minfi", # Method to load data, "ChAMP" method is newly provided by ChAMP group, while "minfi" is old minfi way.(default = "ChAMP")
           methValue = "B", # Indicates whether you prefer m-values M or beta-values B. (default = "B")
           autoimpute = TRUE, # If after filtering (or not do filtering) there are NA values in it, should impute.knn(k=3) should be done for the rest NA?
           filterDetP = TRUE, # If filter = TRUE, then probes above the detPcut will be filtered out.(default = TRUE)
           ProbeCutoff = 0.1, # The NA ratio threshhold for probes. Probes with above proportion of NA will be removed.
           SampleCutoff = 0.1, # The failed p value (or NA) threshhold for samples. Samples with above proportion of failed p value (NA) will be removed.
           detPcut = detPcut, # The detection p-value threshold. Probes about this cutoff will be filtered out. (default = 0.01)
           filterBeads = TRUE, # If filterBeads=TRUE, probes with a beadcount less than 3 will be removed depending on the beadCutoff value.(default = TRUE)
           beadCutoff = 0.05, # The beadCutoff represents the fraction of samples that must have a beadcount less than 3 before the probe is removed.(default = 0.05)
           filterNoCG = TRUE, # If filterNoCG=TRUE, non-cg probes are removed.(default = TRUE)
           filterSNPs = TRUE, # If filterSNPs=TRUE, probes in which the probed CpG falls near a SNP as defined in Nordlund et al are removed.(default = TRUE)
           filterMultiHit = TRUE, # If filterMultiHit=TRUE, probes in which the probe aligns to multiple locations with bwa as defined in Nordlund et al are removed.(default = TRUE)
           filterXY = TRUE, # If filterXY=TRUE, probes from X and Y chromosomes are removed.(default = TRUE)
           force=TRUE
)

passed_cpgs_origin = rownames(myLoad$beta)
RGset <- myLoad$rgSet
Mset <- myLoad$mset
observables = myLoad$pd

rm(myLoad)

# detectionP ==================================================================
detP_before <- detectionP(RGset)

passed_cpgs = intersect(passed_cpgs_origin, rownames(detP_before))
detP_after <- detP_before[passed_cpgs,]

failed_before <- summary(detP_before>detPcut)
failed_before <- data.frame(failed_before)
failed_before <- failed_before[seq(3, nrow(failed_before), 3), 2:3]
colnames(failed_before) <- c("Sample",paste("Number of probes (detP_before>", detPcut, ")",sep=''))
write.csv(failed_before, "failed_before.csv", row.names = FALSE)

failed_after <- summary(detP_after>detPcut)
failed_after <- data.frame(failed_after)
failed_after <- failed_after[seq(3, nrow(failed_after), 3), 2:3]
colnames(failed_after) <- c("Sample",paste("Number of probes (failed_after>", detPcut, ")",sep=''))
write.csv(failed_after, "failed_after.csv", row.names = FALSE)

# QC ==========================================================================
h = 5

pdf("NEGATIVE.pdf", width = 10, height = h)
controlStripPlot(RGset, controls="NEGATIVE")
dev.off()

pdf("BISULFITE CONVERSION I.pdf", width = 10, height = h)
controlStripPlot(RGset, controls="BISULFITE CONVERSION I")
dev.off()

pdf("BISULFITE CONVERSION II.pdf", width = 10, height = h)
controlStripPlot(RGset, controls="BISULFITE CONVERSION II")
dev.off()

pdf("EXTENSION.pdf", width = 10, height = h)
controlStripPlot(RGset, controls="EXTENSION")
dev.off()

pdf("HYBRIDIZATION.pdf", width = 10, height = h)
controlStripPlot(RGset, controls="HYBRIDIZATION")
dev.off()

pdf("NON-POLYMORPHIC.pdf", width = 10, height = h)
controlStripPlot(RGset, controls="NON-POLYMORPHIC")
dev.off()

pdf("SPECIFICITY I.pdf", width = 10, height = h)
controlStripPlot(RGset, controls="SPECIFICITY I")
dev.off()

pdf("SPECIFICITY II.pdf", width = 10, height = h)
controlStripPlot(RGset, controls="SPECIFICITY II")
dev.off()

pdf("TARGET REMOVAL.pdf", width = 10, height = h)
controlStripPlot(RGset, controls="TARGET REMOVAL")
dev.off()

qcReport(RGset,
         pdf = "qcReport.pdf", 
         controls = c("NEGATIVE",
                      "BISULFITE CONVERSION I",
                      "BISULFITE CONVERSION II",
                      "EXTENSION",
                      "HYBRIDIZATION",
                      "NON-POLYMORPHIC",
                      "SPECIFICITY I",
                      "SPECIFICITY II",
                      "TARGET REMOVAL"
         )
)

qc <- getQC(Mset)
pdf("QCplot.pdf")
plotQC(qc)
dev.off()

green <- getGreen(RGset)
red <- getRed(RGset)

pdf("red_intensities.pdf", width=15, height=5)
par(mfrow=c(1,1))
boxplot(red,outline=F)
dev.off()

pdf("green_intensities.pdf", width=15, height=5)
par(mfrow=c(1,1))
boxplot(green,outline=F)
dev.off()

# Raw =========================================================================
raw = preprocessRaw(RGset)
beta_raw <- getBeta(raw)

passed_cpgs = intersect(passed_cpgs_origin, rownames(beta_raw))
beta_raw_filtered <- beta_raw[passed_cpgs,]

pdf("density_raw.pdf")
densityPlot(beta_raw_filtered, sampGroups = observables$Sample_Group)
dev.off()

pdf("boxplot_raw.pdf",width=15,height=5)
par(mfrow=c(1,1))
boxplot(beta_raw_filtered, outline=F, main="No Normalization")
dev.off()

# funnorm =====================================================================
funnorm <- preprocessFunnorm(RGset)
beta_funnorm <- getBeta(funnorm)

passed_cpgs = intersect(passed_cpgs_origin, rownames(beta_funnorm))
beta_funnorm_filtered <- beta_funnorm[passed_cpgs,]

pdf("boxplots_funnorm.pdf",width=15,height=5)
par(mfrow=c(1,1))
boxplot(beta_funnorm_filtered,outline=F,main="Funnorm Normalization")
dev.off()

pdf("density_funnorm.pdf")
densityPlot(beta_funnorm_filtered, sampGroups = observables$Sample_Group)
dev.off()

beta_funnorm_filtered_df <- data.frame(row.names(beta_funnorm_filtered),beta_funnorm_filtered)
colnames(beta_funnorm_filtered_df)[1] <- "CpG"
write.table(beta_funnorm_filtered_df,file="beta_funnorm_filtered.txt",row.names=F,sep="\t",quote=F)

# quantile ====================================================================
quantile <- preprocessQuantile(RGset)
beta_quantile <- getBeta(quantile)

passed_cpgs = intersect(passed_cpgs_origin, rownames(beta_quantile))
beta_quantile_filtered <- beta_quantile[passed_cpgs,]

pdf("boxplots_quantile.pdf",width=15,height=5)
par(mfrow=c(1,1))
boxplot(beta_quantile_filtered,outline=F,main="Quantile Normalization")
dev.off()

pdf("density_quantile.pdf")
densityPlot(beta_quantile_filtered, sampGroups = observables$Sample_Group)
dev.off()

beta_quantile_filtered_df <- data.frame(row.names(beta_quantile_filtered),beta_quantile_filtered)
colnames(beta_quantile_filtered_df)[1] <- "IlmnID"
write.table(beta_quantile_filtered_df,file="beta_quantile_filtered",row.names=F,sep="\t",quote=F)

# BMIQ ========================================================================
myNorm <- champ.norm(beta=myLoad$beta,
                     rgSet=myLoad$rgSet,
                     mset=myLoad$mset,
                     method="BMIQ",
                     arraytype=chip_type)

pdf("boxplots_BMIQ.pdf",width=15,height=5)
par(mfrow=c(1,1))
boxplot(myNorm,outline=F,main="Quantile Normalization")
dev.off()

pdf("density_BMIQ.pdf")
densityPlot(myNorm, sampGroups = observables$Sample_Group)
dev.off()

beta_filtered_normalized = data.frame(row.names(myNorm), myNorm)
colnames(beta_filtered_normalized)[1] <- "IlmnID"
write.table(beta_filtered_normalized,file="beta_BMIQ_filtered.txt",col.name=TRUE, row.names=FALSE,sep="\t",quote=F)

# Cells =======================================================================
cellType = "Blood" # Should be one of "Blood", "CordBloodCombined", "CordBlood", "CordBloodNorway", "CordTissueAndBlood", or "DLPFC"
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu") # c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu") for "Blood" and c("NeuN_neg", "NeuN_pos") for "DLPFC"

refPlatform = "IlluminaHumanMethylation450k"
refPlatform = "IlluminaHumanMethylationEPIC"

cell_counts <- estimateCellCounts2(RGset,
                                   compositeCellType = cellType,
                                   cellTypes = cellTypes,
                                   processMethod = "preprocessFunnorm", 
                                   referencePlatform=refPlatform,
                                   meanPlot = TRUE,
                                   verbose = TRUE,
                                   lessThanOne = FALSE)

write.csv(cell_counts, "cell_counts.csv")
# =============================================================================

