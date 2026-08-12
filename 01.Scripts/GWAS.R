list.of.packages <- c("openxlsx","GWASpoly", "tidyverse","BGLR","data.table","sommer","gtools","adegenet")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]                       
if(length(new.packages)) {install.packages(new.packages)}

invisible(lapply(list.of.packages, library, character.only = TRUE))      

rm(list = ls())

gc()

setwd("E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/Frontiers_Paper/3.4.Revisiones")
# source("codigo CV/CV_paper.R")

phenofile <- "GWAS_BINOMIAL_DSI.csv"
genofile <- "base_genotipica_121168_SNPs_519_genotipos_ACGT.csv"

## GWAS Analysis
data <- read.GWASpoly(ploidy=10,pheno.file=phenofile,geno.file=genofile,format="ACGT",n.traits=11,delim=",")

data = set.K(data = data, n.core = 1, LOCO = FALSE)


params = set.params(fixed = c("Q1","Q2","Q3", "Q4"),
                    fixed.type = c("numeric","numeric","numeric", "numeric"), 
                    n.PC = 0,
                    MAF = NULL,
                    geno.freq = 0.95,
                    P3D = T)

DGWAS <- GWASpoly(data = data, models = c("additive","general","1-dom","2-dom","3-dom","4-dom","5-dom"),
                  traits=c("DSI_BLUP_general"),params = params, 
                  n.core = 1, quiet = F )

Diatraea = DGWAS
Diatraea <- set.threshold(Diatraea,  method="FDR", level = 0.05)

# save(Diatraea, file = "GWASPoly_GWAS_DIATRAEA_DSIndex.RData")
load("GWAS_BLUP_DSI_121166_SNPs_4Q.RData")

QTL <- get.QTL(Diatraea)
