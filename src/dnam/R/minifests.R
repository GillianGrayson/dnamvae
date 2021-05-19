rm(list=ls())

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylation27kmanifest")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("IlluminaHumanMethylation27kanno.ilmn12.hg19")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

library(minfi)
library(openxlsx)
library(IlluminaHumanMethylation27kmanifest)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation27kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

path <- "E:/YandexDisk/Work/dnamvae"

ann27k <- getAnnotation(IlluminaHumanMethylation27kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

fn_save = paste(path, "/ann27k.xlsx",  sep='')
write.xlsx(ann27k, fn_save, sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
fn_save = paste(path, "/ann450k.xlsx",  sep='')
write.xlsx(ann450k, fn_save, sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
fn_save = paste(path, "/ann850k.xlsx",  sep='')
write.xlsx(ann850k, fn_save, sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
