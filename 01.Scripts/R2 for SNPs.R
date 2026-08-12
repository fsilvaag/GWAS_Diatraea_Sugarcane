list.of.packages <- c("openxlsx","GWASpoly", "tidyverse","BGLR","data.table","sommer","gtools","adegenet")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]                       
if(length(new.packages)) {install.packages(new.packages)}

invisible(lapply(list.of.packages, library, character.only = TRUE))      

rm(list = ls())

gc()

#### 0. GWAS Results ####

load("GWAS_BLUP_DSI_121166_SNPs_4Q.RData")
QTL <- get.QTL(Diatraea)

qtl_9_gen <- QTL[,c("Model","Marker")]


#### 1. Maddala's 1983 R2 ####
source("Fit_R2_CSD.R")
resultado_r2 <- fit.QTL_CSD(data = Diatraea,
                            trait = "DSI_BLUP_general",
                            qtl = qtl_9_gen,
                            fixed = NULL)

round(100*resultado_r2$R2_M,2)

#### 2. Cross-Validation R2 ####
source("CV_paper.R")
#### 2.1 general - 10:7656049_T/C ####
seed(456)
prepared_9_general <- prepare_gwas_data(
  data = Diatraea,
  trait = "DSI_BLUP_general",
  qtl = qtl_9_gen[6,],
  fixed = Diatraea@fixed[,c(1:4)])

cv_9_dom <- cross_validate_qtl(
  prepared = prepared_9_general,
  k = 5,
  repeats = 50,
  method = "ML",
  seed = 2026
)


tab =cv_9_dom %>%
  filter(metric == "R2_pred") %>%
  # group_by(Model, Marker) %>% 
  summarise(Full_median = median(full, na.rm = TRUE) * 100,
    Full_P2.5 = quantile(full, 0.025, na.rm = TRUE) * 100,
    Full_P97.5 = quantile(full, 0.975, na.rm = TRUE) * 100,
    Delta_median = median(increment, na.rm = TRUE) * 100,
    Delta_P2.5 = quantile(increment, 0.025, na.rm = TRUE) * 100,
    Delta_P97.5 = quantile(increment, 0.975, na.rm = TRUE) * 100,
    Positive_folds = mean(increment > 0, na.rm = TRUE) * 100,
    .groups = "drop") %>%
  mutate(Full_R2_CV = sprintf(
      "%.2f (%.2f to %.2f)",
      Full_median,Full_P2.5,Full_P97.5),
      Delta_R2_CV = sprintf("%.2f (%.2f to %.2f)",
                            Delta_median,Delta_P2.5,Delta_P97.5),
    Positive_folds = round(Positive_folds, 1))



rmse_cv <- cv_9_dom |> filter(metric == "RMSE") |> 
  summarise(RMSE_reducido_media = mean(reduced, na.rm = TRUE),
            RMSE_reducido_SD = sd(reduced, na.rm = TRUE),
            RMSE_completo_media = mean(full, na.rm = TRUE),
            RMSE_completo_SD = sd(full, na.rm = TRUE),
            Reduccion_RMSE_media = mean(reduced - full,na.rm = TRUE),
            Reduccion_RMSE_SD = sd(reduced - full,na.rm = TRUE),
            Proporcion_mejora = mean(full < reduced,na.rm = TRUE))

rmse_mae <- cv_9_dom |> filter(metric == "MAE") |> 
  summarise(MAE_reducido_media = mean(reduced, na.rm = TRUE),
            MAE_reducido_SD = sd(reduced, na.rm = TRUE),
            MAE_completo_media = mean(full, na.rm = TRUE),
            MAE_completo_SD = sd(full, na.rm = TRUE),
            Reduccion_MAE_media = mean(reduced - full,na.rm = TRUE),
            Reduccion_MAE_SD = sd(reduced - full,na.rm = TRUE),
            Proporcion_mejora = mean(full < reduced,na.rm = TRUE))

rmse_cv
rmse_mae



