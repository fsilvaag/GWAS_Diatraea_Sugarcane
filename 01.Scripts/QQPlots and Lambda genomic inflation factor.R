list.of.packages <- c("openxlsx","GWASpoly", "tidyverse","BGLR","data.table","sommer","gtools","adegenet")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]                       
if(length(new.packages)) {install.packages(new.packages)}

invisible(lapply(list.of.packages, library, character.only = TRUE))      

rm(list = ls()[!ls() %in% c("DGWAS","Diatraea")])

gc()

#### 0. GWAS Results ####

load("GWAS_BLUP_DSI_121166_SNPs_4Q.RData")
QTL <- get.QTL(Diatraea)

#### 1. Inflation Factor ####

markers <- rownames(Diatraea@scores[["DSI_BLUP_general"]])
trait_name <- "DSI_BLUP_general"
modelo <- "general"

# EStimate the degrees of freedom for the general model for each marker.
# This are equal to the number of classes or allele dossages
df_vector <- sapply(markers,function(marker) {
  S <- try(GWASpoly:::.design.score(Diatraea@geno[, marker],
                                    model = modelo,
                                    ploidy = Diatraea@ploidy,
                                    min.MAF = 0,
                                    max.geno.freq = 1),silent = T)
  
  
  if (inherits(S, "try-error")) {
    S <- NULL
  } else {
    ncol(S)  
  }
})

save(df_vector,file = "Degrees_of_freedon_General_Model.RData")

df = bind_rows(df_vector)
df = t(as.data.frame(df))
colnames(df) <- "DF"
df <- as.data.frame(df) |> rownames_to_column("Marker")

scores <- data.frame(Marker = rownames(Diatraea@scores[[trait_name]]),
                     general = Diatraea@scores[[trait_name]][, modelo])
scores <- left_join(scores, df) |> mutate(p_values = 10^(-general),
                                          chisq = ifelse(!is.na(p_values),
                                                         qchisq(p_values, df = DF, lower.tail = F)
                                                         ,NA),
                                          chisq_expected = qchisq(0.5,df = DF))



lambda_gc <- median(scores$chisq / scores$chisq_expected, na.rm = T) 

round(lambda_gc,2)

scores <- scores |> mutate(chisq_corregido = scores$chisq / lambda_gc,
                           p_corregidos = pchisq(chisq_corregido,df = DF,
                                                 lower.tail = FALSE),
                           scores_corregidos = -log10(p_corregidos),
                           FDR_Corregidos = p.adjust(p_corregidos,method = "BH"),
                           FDR_orignal = p.adjust(p_values,method = "BH"))

modelo <- colnames(Diatraea@scores$DSI_BLUP_general)[-2]
tb= data.frame(Modelo = as.character(), Lambda = as.numeric())
for (v in 1:length(modelo)) {
  log_p_values <- na.omit(Diatraea@scores[[trait_name]][, modelo[v]])
  p_values <- 10^(-log_p_values)
  chisq <- qchisq(1 - p_values, df = 1)
  lambda_gc5 <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1) 
  tb = rbind(tb, data.frame(Modelo = modelo[v], Lambda = lambda_gc5))
  
}
tb <- rbind(tb, data.frame(Modelo = "general", Lambda = lambda_gc))

Diatraea@scores[[trait_name]][, modelo] <- scores$scores_corregidos
scores |> filter(FDR_orignal<= 0.05)


#### 2. QQ Plots ####
source("qqPlod_CSD.R")
d <- Diatraea@scores$DSI_BLUP_general
d <- d %>% rownames_to_column("M")
d$SNP <- d$M
d <- d %>% separate(M, into = c("Chromosome","B"), sep = "\\:") |> 
  separate(B, into = c("Position","M"), sep = "_") |> 
  separate(M, into = c("Ref","Alt"), sep = "\\/")
d <- transform(d, Position = as.numeric(Position))
scores <- d
a1 <- qqCDS(data = d,modelos = c("5-dom-alt"))
a1


library(ggplot2)
# setwd("E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/Frontiers_Paper/3.4.Revisiones/")
ggsave(
  filename = "QQPLOT_DSI_1_dom_ref.tiff",
  plot = a1,
  device = "tiff",
  width = 16,
  height = 8,
  units = "in",
  dpi = 500,
  compression = "lzw"
)


source("Manhattan_plot_FSA_24072026.R")


PL = manhattan.plot_FSA(data = Diatraea,traits="DSI_BLUP_general",models="5-dom-alt", 
                   chrom.hold = NULL, solid= TRUE, 
                   umbral = 0.05, 
                   nme = "CS",map = NULL,
                   mkr = c("3:18796998_T/A","9:9208699_T/C"),
                   plot = F, xtes= 15, hj = 0.5, lp = "none",
                   at = 20, axt = 20, scx = "Chromosome",
                   scm = c("#EC5f67", "#FAC863",'#99C794','#6699CC','#C594C5'),
                   ma = NULL, mi = NULL, us = 1, uc = "darkgreen",
                   snpsize = 5, snppadd = 0.6, bpad = 0.6, segcol = "grey25",
                   nombre_SNP=c("CC-SNP1","CC-SNP3"))


PL
ggsave(
  filename = "3.5.Figures/Manhattan_DSI_5domalt.tiff",
  plot = PL$DSI_BLUP_general,
  device = "tiff",
  width = 16,
  height = 8,
  units = "in",
  dpi = 300,
  compression = "lzw"
)



library(magick)
# setwd("E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/Frontiers_Paper")
img1 <- image_read("3.5.Figures/Manhattan_DSI_general.tiff")
img2 <- image_read("3.3.GWAS/QQPLOT_DSI_general.tiff")

img1 <- image_trim(image_background(img1, color = "white", flatten = TRUE))
img2 <- image_trim(image_background(img2, color = "white", flatten = TRUE))
img1_resized <- image_resize(img1, "2000x")
img2_resized <- image_resize(img2, "2000x")

img1_con_espacio <- image_border(img1_resized, color = "white", geometry = "0x50")
img1_con_espacio <- image_crop(img1_con_espacio, geometry = "2000x")

img_final <- image_append(c(img1_con_espacio, img2_resized), stack = TRUE)
img_final <- image_border(img_final, color = "white", geometry = "40x0")
# img_final <- image_trim(img_final)

image_write(
  img_final, 
  path = "3.5.Figures/Figure_3_general_Hor.tiff", 
  format = "tiff", 
  density = "300", 
  compression = "lzw"
)
