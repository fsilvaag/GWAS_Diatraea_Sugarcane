list.of.packages <- c("openxlsx","GWASpoly", "tidyverse","data.table","gtools","ggrepel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]                       
if(length(new.packages)) {install.packages(new.packages)}

invisible(lapply(list.of.packages, library, character.only = TRUE))      

rm(list = ls())

gc()

library(data.table)
# phenofile <- "//192.168.153.238/biodata10/rstudio/fernando/02.GWAS/Sanidad/Diatraea/Diatraea__solo_220_121k_combinado.csv"
phenofile <- "E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/01.Abril2026/Result/BLUES_DSIndex.csv"
genofile <- "//192.168.153.238/biodata10/rstudio/fernando/02.GWAS/papers_2025_sanidad/base_genotipica_121168_SNPs_519_genotipos_ACGT_lista_gwas.csv"

# 
# source("//192.168.153.238/biodata10/rstudio/fernando/20.Funciones/READ_GWAS_SNP_Remover1.R")
# # 
# # # n <- count.fields(phenofile, sep = ",")[1] - 6 # 4 = Name + Q1 + Q2 + Q3
# # 
# # 
# data <- GWAS_CSD(ploidy = 10,pheno.file = phenofile, geno.file = genofile,
# #                  format = "ACGT",
# #                  n.traits = 34,delim = ",",
#                  snp_remover = NULL,colIDy = 1 ,imp = "mode")

library(GWASpoly)
data <- read.GWASpoly(ploidy=10,pheno.file=phenofile,geno.file=genofile,
                      format="ACGT",n.traits=8,delim=",")



data = set.K(data = data, n.core = 10, LOCO = FALSE)


params = set.params(fixed = c("Q1","Q2","Q3","Q4"),
                    fixed.type = c("numeric","numeric","numeric","numeric"), 
                    n.PC = 0,
                    MAF = NULL,
                    geno.freq = 0.95,
                    P3D = T)

DGwas <- GWASpoly(data = data, models = c("additive","general","1-dom","2-dom","3-dom","4-dom","5-dom"),
                                                 traits=colnames(data@pheno)[-1],params = params, n.core = 1, quiet = F )

# save(DGwas, file = "/biodata10/rstudio/fernando/02.GWAS/Sanidad/Diatraea/05.Enero2026/GWASPoly_Diatraea__solo_220_121166_BLUes_BLUPs_Combinado_individual.RData")
load("E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/Frontiers_Paper/3.3.GWAS/GWASPoly_Diatraea_Discovery_Population_220genotypes_121166_DSIndex_BLUPs_Combinado.RData")

# D_FDR<- set.threshold(DGwas,  method="FDR")

# QTL_FDR<- get.QTL(D_FDR)

nm = "DSI_Comb_BLUP"
fdr = NULL
for (i in 1:length(nm)) {
  a <- DGwas@scores[[nm[i]]]
  a1 = apply(a, 2, function(dfad) p.adjust(10^-dfad, "fdr"))
  a1 = as.data.frame(a1)
  a1$marker = rownames(a1)
  a1 <- a1 %>% remove_rownames()
  a2 <- a1 %>% pivot_longer(-'marker', names_to = "Modelo", values_to = "pj") %>% 
    mutate(Trait = nm[i])
  fdr = rbind(fdr,a2)
  rm(a)
  rm(a1)
  rm(a2)
  
}

fdr1_1 <- fdr %>% filter(pj <= 0.05) 

## Marker effect for General Model
source("GWASpoly_General_efectos.R")

gg = GWASpoly_general(data= DGwas,models = c("general"),
                      traits="DSI_Comb_BLUP",params=params,n.core=1,quiet=F)

# save(gg, file = "E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/Frontiers_Paper/3.3.GWAS/Marker_Effects_General_Model.RData")

## Manhattan Plots
source("E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/01.Abril2026/codes for association analysis/Manhattan_plot_FSA_CSD_qValue.R")

mkr = c('5:4170778_T/G','6:18496137_G/A','6:33057028_A/G','9:9208699_T/C','contig_16578:78638_A/G','contig_3381:85998_A/C')

tiff(filename = "E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/Frontiers_Paper/3.3.GWAS/manhattan_plot_model_general_Validated.tif",
     width = 19,height = 9,units = "in",res = 600)
manhattan.plot_FSA(data = DGwas, 
  map = NULL,
  traits = "DSI_Comb_BLUP",        
  models = "general",  mkr = mkr,          
  chrom.hold = NULL, solid = TRUE,
  umbral = 0.05,umb2 = NULL, nme = "Contigs",
  plot = F, xtes= 15, hj = 0.5, lp = "none",
  at = 15, axt = 14, scx = "Chromosome",
  scm = c("#EC5f67", "#FAC863",'#99C794','#6699CC','#C594C5'),
  ma = NULL, mi = NULL, us = 1, uc = "darkgreen",
  snpsize = 4,snppadd = 0.2,bpad = 0.5,segcol = "grey25")
dev.off()

mkr = c('5:4170778_T/G')

tiff(filename = "E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/Frontiers_Paper/3.3.GWAS/manhattan_plot_model_additive_Validated.tif",
     width = 19,height = 9,units = "in",res = 600)
manhattan.plot_FSA(data = DGwas, 
                   map = NULL,
                   traits = "DSI_Comb_BLUP",        
                   models = "additive",  mkr = mkr,          
                   chrom.hold = NULL, solid = TRUE,
                   umbral = 0.05,umb2 = NULL, nme = "Contigs",
                   plot = F, xtes= 15, hj = 0.5, lp = "none",
                   at = 15, axt = 14, scx = "Chromosome",
                   scm = c("#EC5f67", "#FAC863",'#99C794','#6699CC','#C594C5'),
                   ma = NULL, mi = NULL, us = 1, uc = "darkgreen",
                   snpsize = 4,snppadd = 0.2,bpad = 0.5,segcol = "grey25")
dev.off()
