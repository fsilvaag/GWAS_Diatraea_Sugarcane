list.of.packages <- c("openxlsx","lme4", "tidyverse","data.table","gtools","Matrix",
                      "purrr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]                       
if(length(new.packages)) {install.packages(new.packages)}

invisible(lapply(list.of.packages, library, character.only = TRUE))      

rm(list = ls())

gc()

setwd("E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/Frontiers_Paper/3.4.Revisiones")
dt = read.xlsx("../../01.Abril2026/Raw_Phenotypic_Data_Base_Diatraea.xlsx")

dt = dt |> transform(Genotype = as.factor(Genotype), eval = as.factor(eval))
dt <- dt |> mutate(ind1 = (nlarvas + npupas + nexuvias))
dtdisco <- dt |> filter(Population =='Discovery') |> droplevels()

#---------------------------------------------------#
# SInd
#---------------------------------------------------#
calcular_h2_cullis_binomial <- function(datos) {
  
  datos <- datos |> filter(!is.na(ind1),
                           !is.na(larvas),
                           !is.na(Genotype),
                           larvas > 0,
                           ind1 >= 0,
                           ind1 <= larvas) |>
    mutate(Genotype = droplevels(factor(Genotype)))
  
  evaluaciones <- unique(datos$eval)
  
  evaluacion <- evaluaciones[1]
  
  modelo <- glmer(cbind(ind1, larvas - ind1) ~ 1 + (1 | Genotype),
    data = datos,
    family = binomial(link = "logit"),
    nAGQ = 1,
    control = glmerControl(optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)))
  
  mensajes <- modelo@optinfo$conv$lme4$messages
  
  converged <- is.null(mensajes)
  singular  <- isSingular(modelo, tol = 1e-5)
  
  vc <- as.data.frame(VarCorr(modelo))
  
  sigma_g2 <- vc |> filter(grp == "Genotype",
      var1 == "(Intercept)") |> pull(vcov)
  
  if (length(sigma_g2) != 1 || sigma_g2 <= 0) {
    stop("No se obtuvo una varianza genotípica positiva.")
  }
  
  X <- getME(modelo, "X")
  Z <- getME(modelo, "Z")
  
  X <- Matrix(X, sparse = TRUE)
  Z <- Matrix(Z, sparse = TRUE)
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  p_hat <- fitted(modelo)
  
  epsilon <- 1e-10
  
  p_hat <- pmin(pmax(p_hat, epsilon),1 - epsilon)
  
  # ------------------------------------------------------------
  # Matriz de pesos binomiales PIRLS
  #
  # Para datos binomiales agrupados:
  # W_i = larvas_i * p_i * (1 - p_i)
  # ------------------------------------------------------------
  
  w <- datos$larvas * p_hat * (1 - p_hat)
  
  W <- Diagonal(n = n,x = w)
  
  # ------------------------------------------------------------
  # Matriz G para los efectos genotípicos
  #
  # G = sigma_g^2 * I
  # G^-1 = I / sigma_g^2
  # ------------------------------------------------------------
  
  G_inv <- Diagonal(n = q,x = rep(1 / sigma_g2, q))
  
  # ------------------------------------------------------------
  # 8. Matriz de ecuaciones mixtas penalizadas
  #
  #          X'WX        X'WZ
  # C = 
  #          Z'WX   Z'WZ + G^-1
  # ------------------------------------------------------------
  
  C11 <- crossprod(X, W %*% X)
  C12 <- crossprod(X, W %*% Z)
  C21 <- t(C12)
  C22 <- crossprod(Z, W %*% Z) + G_inv
  
  C <- rbind(cbind(C11, C12),cbind(C21, C22))
  
  C_inv <- tryCatch(solve(C),error = function(e) {
      message("Se utilizará una inversa generalizada.")
      Matrix(MASS::ginv(as.matrix(C)), sparse = FALSE)
    }
  )
  
  # ------------------------------------------------------------
  # BLUP de Genotype
  # ------------------------------------------------------------
  
  indices_u <- (p + 1):(p + q)
  
  C22_g <- as.matrix(C_inv[indices_u, indices_u, drop = FALSE])
  
  # Orden de los niveles de Genotype en Z
  niveles_Genotype <- levels(getME(modelo, "flist")$Genotype)
  
  if (length(niveles_Genotype) == nrow(C22_g)) {
    rownames(C22_g) <- niveles_Genotype
    colnames(C22_g) <- niveles_Genotype
  }
  
  # ------------------------------------------------------------
  # BLUP de los Genotypes
  # ------------------------------------------------------------
  
  tabla_blup <- ranef(modelo)$Genotype |>
    as.data.frame() |>
    rownames_to_column("Genotype") |>
    rename(BLUP = `(Intercept)`)
  
  # ------------------------------------------------------------
  # 12. Varianza de cada diferencia entre BLUP
  #
  # Var(e_i - e_j) = C_ii + C_jj - 2*C_ij
  # ------------------------------------------------------------
  
  pares <- combn(seq_len(q), 2)
  
  tabla_diferencias <- map_dfr(seq_len(ncol(pares)), function(k) {
      i <- pares[1, k]
      j <- pares[2, k]
      
      var_diferencia <- C22_g[i, i] + C22_g[j, j] - 2 * C22_g[i, j]
      
      tibble(eval = evaluacion,
             Genotype_1 = niveles_Genotype[i],
             Genotype_2 = niveles_Genotype[j],
             BLUP_1 = tabla_blup$BLUP[match(niveles_Genotype[i], tabla_blup$Genotype)],
             BLUP_2 = tabla_blup$BLUP[match(niveles_Genotype[j], tabla_blup$Genotype)],
             Difference_BLUP = BLUP_1 - BLUP_2,
             VarianceDifference = var_diferencia,
             SE_Difference = sqrt(max(var_diferencia, 0)))
    }
  )
  
  # ------------------------------------------------------------
  # Media de las varianzas de diferencias
  # ------------------------------------------------------------
  
  vbar_blup <- mean(tabla_diferencias$VarianceDifference,na.rm = TRUE)
  
  # ------------------------------------------------------------
  # Heredabilidad tipo Cullis
  #
  # H2 = 1 - vbar_BLUP / (2*sigma_g^2)
  # ------------------------------------------------------------
  
  H2_Cullis <- 1 - vbar_blup / (2 * sigma_g2)
  
  resultado <- tibble(eval = evaluacion,
                      N_observations = nrow(datos),
                      N_genotypes = q,
                      N_pairs = nrow(tabla_diferencias),
                      GeneticVariance = sigma_g2,
                      MeanVarianceDifference = vbar_blup,
                      H2_Cullis = H2_Cullis,
                      Singular = singular,
                      Converged = converged,
                      ConvergenceMessage = if (converged) {
                        NA_character_ } else {
                          paste(mensajes, collapse = "; ") 
                          }
                      )
    list(result = resultado,
         pairwise = tabla_diferencias,
         blup = tabla_blup,
         C22_genotype = C22_g,
         model = modelo)
}


nm = levels(dtdisco$eval)


# res_eval_1 <- calcular_h2_cullis_binomial(dt_eval1)
# 
# res_eval_1$result
# 
# 
resultados_por_eval <- split(dtdisco, dtdisco$eval) |>
  purrr::map(calcular_h2_cullis_binomial)

tabla_h2_cullis <- resultados_por_eval |>
  map_dfr("result")


#---------------------------------------------------#
# BI
#---------------------------------------------------#
calcular_h2_cullis_binomial_BI <- function(datos) {
  
  datos <- datos |> filter(!is.na(TE),
                           !is.na(ED),
                           !is.na(Genotype)) |>
    mutate(Genotype = droplevels(factor(Genotype)))
  
  evaluaciones <- unique(datos$eval)
  
  evaluacion <- evaluaciones[1]
  
  modelo <- glmer(cbind(ED, TE - ED) ~ 1 + (1 | Genotype),
                  data = datos,
                  family = binomial(link = "logit"),
                  nAGQ = 1,
                  control = glmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 200000)))
  
  mensajes <- modelo@optinfo$conv$lme4$messages
  
  converged <- is.null(mensajes)
  singular  <- isSingular(modelo, tol = 1e-5)
  
  vc <- as.data.frame(VarCorr(modelo))
  
  sigma_g2 <- vc |> filter(grp == "Genotype",
                           var1 == "(Intercept)") |> pull(vcov)
  
  if (length(sigma_g2) != 1 || sigma_g2 <= 0) {
    stop("No se obtuvo una varianza genotípica positiva.")
  }
  
  X <- getME(modelo, "X")
  Z <- getME(modelo, "Z")
  
  X <- Matrix(X, sparse = TRUE)
  Z <- Matrix(Z, sparse = TRUE)
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  p_hat <- fitted(modelo)
  
  epsilon <- 1e-10
  
  p_hat <- pmin(pmax(p_hat, epsilon),1 - epsilon)
  
  # ------------------------------------------------------------
  # Matriz de pesos binomiales PIRLS
  #
  # Para datos binomiales agrupados:
  # W_i = larvas_i * p_i * (1 - p_i)
  # ------------------------------------------------------------
  
  w <- datos$larvas * p_hat * (1 - p_hat)
  
  W <- Diagonal(n = n,x = w)
  
  # ------------------------------------------------------------
  # Matriz G para los efectos genotípicos
  #
  # G = sigma_g^2 * I
  # G^-1 = I / sigma_g^2
  # ------------------------------------------------------------
  
  G_inv <- Diagonal(n = q,x = rep(1 / sigma_g2, q))
  
  # ------------------------------------------------------------
  # 8. Matriz de ecuaciones mixtas penalizadas
  #
  #          X'WX        X'WZ
  # C = 
  #          Z'WX   Z'WZ + G^-1
  # ------------------------------------------------------------
  
  C11 <- crossprod(X, W %*% X)
  C12 <- crossprod(X, W %*% Z)
  C21 <- t(C12)
  C22 <- crossprod(Z, W %*% Z) + G_inv
  
  C <- rbind(cbind(C11, C12),cbind(C21, C22))
  
  C_inv <- tryCatch(solve(C),error = function(e) {
    message("Se utilizará una inversa generalizada.")
    Matrix(MASS::ginv(as.matrix(C)), sparse = FALSE)
  }
  )
  
  # ------------------------------------------------------------
  # BLUP de Genotype
  # ------------------------------------------------------------
  
  indices_u <- (p + 1):(p + q)
  
  C22_g <- as.matrix(C_inv[indices_u, indices_u, drop = FALSE])
  
  # Orden de los niveles de Genotype en Z
  niveles_Genotype <- levels(getME(modelo, "flist")$Genotype)
  
  if (length(niveles_Genotype) == nrow(C22_g)) {
    rownames(C22_g) <- niveles_Genotype
    colnames(C22_g) <- niveles_Genotype
  }
  
  # ------------------------------------------------------------
  # BLUP de los Genotypes
  # ------------------------------------------------------------
  
  tabla_blup <- ranef(modelo)$Genotype |>
    as.data.frame() |>
    rownames_to_column("Genotype") |>
    rename(BLUP = `(Intercept)`)
  
  # ------------------------------------------------------------
  # 12. Varianza de cada diferencia entre BLUP
  #
  # Var(e_i - e_j) = C_ii + C_jj - 2*C_ij
  # ------------------------------------------------------------
  
  pares <- combn(seq_len(q), 2)
  
  tabla_diferencias <- map_dfr(seq_len(ncol(pares)), function(k) {
    i <- pares[1, k]
    j <- pares[2, k]
    
    var_diferencia <- C22_g[i, i] + C22_g[j, j] - 2 * C22_g[i, j]
    
    tibble(eval = evaluacion,
           Genotype_1 = niveles_Genotype[i],
           Genotype_2 = niveles_Genotype[j],
           BLUP_1 = tabla_blup$BLUP[match(niveles_Genotype[i], tabla_blup$Genotype)],
           BLUP_2 = tabla_blup$BLUP[match(niveles_Genotype[j], tabla_blup$Genotype)],
           Difference_BLUP = BLUP_1 - BLUP_2,
           VarianceDifference = var_diferencia,
           SE_Difference = sqrt(max(var_diferencia, 0)))
  }
  )
  
  # ------------------------------------------------------------
  # Media de las varianzas de diferencias
  # ------------------------------------------------------------
  
  vbar_blup <- mean(tabla_diferencias$VarianceDifference,na.rm = TRUE)
  
  # ------------------------------------------------------------
  # Heredabilidad tipo Cullis
  #
  # H2 = 1 - vbar_BLUP / (2*sigma_g^2)
  # ------------------------------------------------------------
  
  H2_Cullis <- 1 - vbar_blup / (2 * sigma_g2)
  
  resultado <- tibble(eval = evaluacion,
                      N_observations = nrow(datos),
                      N_genotypes = q,
                      N_pairs = nrow(tabla_diferencias),
                      GeneticVariance = sigma_g2,
                      MeanVarianceDifference = vbar_blup,
                      H2_Cullis = H2_Cullis,
                      Singular = singular,
                      Converged = converged,
                      ConvergenceMessage = if (converged) {
                        NA_character_ } else {
                          paste(mensajes, collapse = "; ") 
                        }
  )
  list(result = resultado,
       pairwise = tabla_diferencias,
       blup = tabla_blup,
       C22_genotype = C22_g,
       model = modelo)
}

resultados_por_eval_BI <- split(dtdisco, dtdisco$eval) |>
  purrr::map(calcular_h2_cullis_binomial_BI)

tabla_h2_cullis_BI <- resultados_por_eval_BI |>
  map_dfr("result")

### Modelo combinado
modelo <- glmer(cbind(ind1,larvas-ind1) ~ eval + (1|Genotype:eval), data = dtdisco, 
                family = binomial(link = "logit"),
                control = glmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 200000)))

summary(modelo)


#### DSIndex
library(glmmTMB)

datos_BI <- dtdisco %>%
  transmute(Genotype,
    eval,
    mesa,
    num,
    Trait = "BI",
    Exitos = ED,
    Fracasos = TE - ED
  )

datos_SInd <- dtdisco %>%
  transmute(
    Genotype,
    eval,
    mesa,num,
    Trait = "SInd",
    Exitos = ind1,
    Fracasos = larvas - ind1
  )

dt_long <- bind_rows(datos_BI, datos_SInd) %>%
  mutate(
    Genotype = factor(Genotype),
    eval = factor(eval),
    Trait = factor(Trait, levels = c("BI", "SInd"))
  )

datos_eval <- droplevels(subset(dt_long, eval == "Diatraea_SAM_1"))

base_eval = datos_eval
  
  
  
# ------------------------------------------------
# 2. Función para una evaluación
# ------------------------------------------------

analizar_evaluacion <- function(base_eval) {
  
  base_eval <- droplevels(base_eval)
  
  nombre_eval <- as.character(unique(base_eval$eval))
  
  n_Genotypes <- nlevels(base_eval$Genotype)
  
  modelo <- glmmTMB(cbind(Exitos, Fracasos) ~
      0 + Trait + us(0 + Trait | Genotype),
    family = binomial(link = "logit"),
    data = base_eval,
    control = glmmTMBControl(
      optimizer = optim,
      optArgs = list(method = "BFGS")
    )
  )
  
  convergio <- isTRUE(modelo$sdr$pdHess)
  
  mensaje_convergencia <- modelo$fit$message
  
  vc <- VarCorr(modelo)
  G_aux <- vc$cond$Genotype
  sd_g <- attr(G_aux,"stddev")
  cor_g <- attr(G_aux,"correlation")
  
  # Revisar nombres
  nombres_efectos <- names(sd_g)
  nombre_BI <- grep("BI",nombres_efectos,value = TRUE)[1]
  nombre_SInd <- grep("SInd",nombres_efectos,value = TRUE)[1]
  
  sd_BI <- unname(sd_g[nombre_BI])
  sd_SInd <- unname(sd_g[nombre_SInd])
  rG <- unname(cor_g[nombre_BI, nombre_SInd])
  
  varG_BI <- sd_BI^2
  varG_SInd <- sd_SInd^2
  covG_BI_SInd <- (rG *sd_BI *sd_SInd)
  
  G <- matrix(c(varG_BI,covG_BI_SInd,covG_BI_SInd,varG_SInd),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("BI", "SInd"),c("BI", "SInd")))
  
  pesos <- matrix(c(0.5, 0.5),ncol = 1)
  
  varG_DS <- as.numeric(t(pesos) %*%G %*%pesos)
  
  re <- ranef(modelo,condVar = TRUE)$cond$Genotype
  
  condVar <- attr(re,"condVar")
  nombre_re_BI <- grep("BI",colnames(re),value = TRUE)[1]
  nombre_re_SInd <- grep("SInd",colnames(re),value = TRUE)[1]
  
  indice_BI <- match(nombre_re_BI,colnames(re))
  indice_SInd <- match(nombre_re_SInd,colnames(re))
  Genotypes <- rownames(re)
  n <- nrow(re)
  
  BLUP_BI <- re[[nombre_re_BI]]
  BLUP_SInd <- re[[nombre_re_SInd]]
  
  BLUP_DS <- (0.5 * BLUP_BI +0.5 * BLUP_SInd)
  PEV_DS <- numeric(n)
  PEV_BI <- numeric(n)
  PEV_SInd <- numeric(n)
  PECov_BI_SInd <- numeric(n)
  
  for (i in seq_len(n)) {
    Vi_completa <- condVar[, , i]
    Vi <- Vi_completa[c(indice_BI, indice_SInd),c(indice_BI, indice_SInd),drop = FALSE]
    
    PEV_BI[i] <- Vi[1, 1]
    PEV_SInd[i] <- Vi[2, 2]
    PECov_BI_SInd[i] <- 0
    
    PEV_DS[i] <- 0.25 * PEV_BI[i] + 0.25 * PEV_SInd[i]
  }
  
  
  pares <- combn(seq_len(n),2)
  
  PEV_dif_DS <- apply(pares,2,function(z) {
      
      i <- z[1]
      j <- z[2]
      
      PEV_DS[i] + PEV_DS[j]
    })
  
  mean_PEV_dif_DS <- mean(PEV_dif_DS,na.rm = TRUE)
  
  if (is.na(varG_DS) || varG_DS <= 0) {
    
    H2_DS_Cullis <- NA_real_
    
  } else {
    
    H2_DS_Cullis <- 1 - mean_PEV_dif_DS / (2 * varG_DS)
  }
  
  PEV_dif_BI <- apply(pares,2,function(z) {
      i <- z[1]
      j <- z[2]
      PEV_BI[i] + PEV_BI[j]
    })
  
  mean_PEV_dif_BI <- mean(PEV_dif_BI,na.rm = TRUE)
  
  H2_BI_Cullis <- if (is.na(varG_BI) || varG_BI <= 0) {
    NA_real_
  } else {
    1 - mean_PEV_dif_BI / (2 * varG_BI)
  }
  
  PEV_dif_SInd <- apply(pares,2,function(z) {
      i <- z[1]
      j <- z[2]
      PEV_SInd[i] + PEV_SInd[j]
    })
  
  mean_PEV_dif_SInd <- mean(PEV_dif_SInd,na.rm = TRUE)
  
  H2_SInd_Cullis <- if (is.na(varG_SInd) || varG_SInd <= 0) {
    NA_real_
  } else {
    1 - mean_PEV_dif_SInd /(2 * varG_SInd)
  }
  
  tabla_blups <- tibble(
    eval = nombre_eval,
    Genotype = Genotypes,
    BLUP_BI_logit = BLUP_BI,
    BLUP_SInd_logit = BLUP_SInd,
    BLUP_DS_logit = BLUP_DS,
    PEV_BI = PEV_BI,
    PEV_SInd = PEV_SInd,
    PECov_BI_SInd = PECov_BI_SInd,
    PEV_DS = PEV_DS
  )
  
  resumen <- tibble(
    Evaluacion = nombre_eval,
    n_Genotypes = n,
    convergio = convergio,
    varG_BI = varG_BI,
    varG_SInd = varG_SInd,
    covG_BI_SInd = covG_BI_SInd,
    corG_BI_SInd = rG,
    varG_DS = varG_DS,
    mean_PEV_dif_BI = mean_PEV_dif_BI,
    mean_PEV_dif_SInd = mean_PEV_dif_SInd,
    mean_PEV_dif_DS = mean_PEV_dif_DS,
    H2_Cullis_BI = H2_BI_Cullis,
    H2_Cullis_SInd = H2_SInd_Cullis,
    H2_Cullis_DS = H2_DS_Cullis
  )
  
  list(
    modelo = modelo,
    matriz_G = G,
    blups = tabla_blups,
    resumen = resumen,
    mensaje_convergencia = mensaje_convergencia
  )
}

# Separar por evaluación

resultados_por_eval <- split(dt_long, dt_long$eval) |>
  purrr::map(analizar_evaluacion)


tabla_h2_cullis_DSI <- resultados_por_eval |>
  map_dfr("resumen")



fwrite(tabla_h2_cullis_DSI, "Heredabilidades_DSIndex.csv")
fwrite(tabla_h2_cullis, "Heredabilidades_SInd.csv")
fwrite(tabla_h2_cullis_BI, "Heredabilidades_BI.csv")



