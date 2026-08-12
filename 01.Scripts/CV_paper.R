library(GWASpoly)
library(rrBLUP)

# ============================================================
# 1. Construir el diseño del marcador según el modelo GWASpoly
# ============================================================

make_marker_design <- function(data, marker, model) {
  
  stopifnot(marker %in% data@map$Marker)
  
  S <- GWASpoly:::.design.score(
    data@geno[, marker],
    model = model,
    ploidy = data@ploidy,
    min.MAF = 0,
    max.geno.freq = 1
  )
  
  if (is.null(S)) {
    stop("No fue posible construir el diseño para el marcador: ", marker)
  }
  
  S <- as.matrix(S)
  rownames(S) <- rownames(data@geno)
  
  return(S)
}


# ============================================================
# 2. Construir matriz de efectos fijos
# ============================================================

make_fixed_design <- function(data, rows, fixed = NULL,efecto = NULL) {
  
  X <- matrix(1,
    nrow = length(rows),
    ncol = 1
  )
  
  colnames(X) <- "(Intercept)"
  
  if (!is.null(fixed)) {
    
    for (i in seq_len(ncol(fixed))) {
      
      efe <- fixed[,i]
      tipo <- efecto[i]
      
      xx <- as.numeric(data@fixed[rows, i])
      X <- cbind(X, xx)
      colnames(X)[ncol(X)] <- colnames(data@fixed)[i]
      
      }
    }
  
  
  return(X)
}


# ============================================================
# 3. Extraer matriz K
# ============================================================

get_K_matrix <- function(data, marker = NULL) {
  
  if (length(data@K) == 1) {
    return(as.matrix(data@K[[1]]))
  }
  
  if (is.null(marker)) {
    stop("Debe indicar un marcador cuando se usa LOCO.")
  }
  
  chrom <- as.character(
    data@map$Chrom[match(marker, data@map$Marker)]
  )
  
  exclude_chr <- match(
    chrom,
    levels(data@map$Chrom)
  )
  
  K <- GWASpoly:::makeLOCO(
    data@K,
    exclude = exclude_chr
  )
  
  return(as.matrix(K))
}


# ============================================================
# 4. Preparar datos para validación cruzada
# ============================================================

prepare_gwas_data <- function(data,
                              trait,
                              qtl,
                              fixed = NULL,
                              K = NULL,
                              efecto = NULL) {
  
  stopifnot(inherits(data, "GWASpoly.K"))
  
  if (!trait %in% colnames(data@pheno)) {
    stop("El trait no aparece en data@pheno.")
  }
  
  rows <- which(!is.na(data@pheno[, trait])) # Identify NA data
  
  y <- as.numeric(data@pheno[rows, trait]) # Remove all NA data
  
  gid <- as.character(data@pheno[rows, 1]) # Bring the genotype names
  
  geno_ids <- rownames(data@geno)# Bring the genotype names
  
  idx_geno <- match(gid,geno_ids) # Check that all genotypes has phenotypic and genotypic data
  
  if (anyNA(idx_geno)) {
    stop(
      "Some genotypes from the phenofile are not present in the genotypic matrix."
    )
  }
  
  if (anyDuplicated(gid)) {
    stop(
      "There are some duplicated genotypes in data@pheno. There should be only one row per genotype."
    )
  }
  Xbase <- make_fixed_design(data = data,rows = rows,fixed = fixed) 
  # Creates a X design matrix with fixed effects (e.g., Population structure) and intercept
  # for all genotypes that has data points (rows)
  
  if (is.null(K)) {
    K <- get_K_matrix(data = data,marker = qtl$Marker[1])
  }
  
  K <- as.matrix(K)
  
  # Ordenar K según los individuos del fenotipo
  if (!is.null(rownames(K)) && all(gid %in% rownames(K))) {
    
    K <- K[gid, gid, drop = FALSE]
    } else {
      K <- K[idx_geno, idx_geno, drop = FALSE]
      rownames(K) <- gid
      colnames(K) <- gid
  }
  
  # Calcular la matriz genotipica de acuerdo con el modelo genetico
  diseños <- lapply(seq_len(nrow(qtl)),function(i) {
      
      S <- make_marker_design(data = data,marker = qtl$Marker[i],
                              model = qtl$Model[i])
      S[gid, , drop = FALSE]
    }
  )
  
  Dosis <- do.call(cbind,diseños) # Matriz de dosis alelica segun mdelo genetico
  
  nombres_qtl <- unlist(lapply(seq_len(nrow(qtl)),function(i) {
    paste0(qtl$Marker[i],"_",qtl$Model[i],"_",1+seq_len(ncol(diseños[[i]])))
    }))
  
  colnames(Dosis) <- nombres_qtl
  md = qtl$Model[1]
  return(list(y = y,gid = gid,Xbase = Xbase,Dosis = Dosis,K = K,rows = rows, Modelo = md))
}


library(rrBLUP)

# ============================================================
# 5. Crear particiones de validación cruzada
# ============================================================

make_cv_folds <- function(n,k = 5,repeats = 50,seed = 123) {
  # n =  Numero total de individuos de la poblacion 
  set.seed(seed)
  folds <- vector("list", repeats)
  for (r in seq_len(repeats)) {
    folds[[r]] <- sample(
      rep(seq_len(k), length.out = n)
    )
  }
  
  return(folds)
}


# ============================================================
# 6. Ajustar el modelo en entrenamiento y predecir prueba
# ============================================================

fit_predict_mixed <- function(y,X,K,train,test,method = "ML", Dosis = NULL) {
  
  y_train <- y[train]
  X_train <- X[train, , drop = FALSE]
  X_test  <- X[test, , drop = FALSE]
  K_train <- K[train, train, drop = FALSE]
  K_cross <- K[test, train, drop = FALSE]
  # Eliminar columnas sin variación en entrenamiento,
  # excepto el intercepto
  keep <- apply(X_train,2,function(z) {length(unique(z[!is.na(z)])) > 1})
  keep[1] <- TRUE
  
  X_train <- X_train[, keep, drop = FALSE]
  X_test  <- X_test[, keep, drop = FALSE]
  
  # Eliminar columnas linealmente dependientes
  
  qr_x <- qr(X_train)
  columnas_independientes <- sort(qr_x$pivot[seq_len(qr_x$rank)])
  
  X_train <- X_train[,columnas_independientes,drop = FALSE]
  
  X_test <- X_test[,columnas_independientes,drop = FALSE]
  
  fit <- rrBLUP::mixed.solve(y = y_train,X = X_train,Z = diag(length(train)),
    K = K_train,method = method)
 
  LogL <- (fit$LL)
  if (!purrr::is_empty(Dosis)) {
    df = ncol(Dosis) # Columnas de la matriz de dosis
  } 
  n = length(y_train)
  
 
  beta <- as.numeric(fit$beta)
  
  lambda <- fit$Ve / fit$Vu
  residual_fijo <- as.numeric(y_train - X_train %*% beta)
  
  A <- K_train + diag(lambda, nrow(K_train))
  
  alpha <- tryCatch(solve(A, residual_fijo),error = function(e) {
    qr.solve(A, residual_fijo)
    })
  
  u_test <- as.numeric(K_cross %*% alpha)
  
  pred_test <- as.numeric(X_test %*% beta + u_test  )
  
  return(list(observed = y[test],predicted = pred_test,fit = fit, n = n,df =df, LogL=LogL))
}


# ============================================================
# 7. Calcular métricas de predicción
# ============================================================

prediction_metrics <- function(observed,predicted) {
  
  ok <- complete.cases(observed,predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  
  if (length(observed) < 3) {
    return(c(cor2 = NA_real_,R2_pred = NA_real_,RMSE = NA_real_,
        MAE = NA_real_))
    }
  
  cor2 <- suppressWarnings(cor(observed, predicted)^2)
  
  SSE <- sum((observed - predicted)^2)
  
  SST <- sum((observed - mean(observed))^2)
  
  R2_pred <- ifelse(SST > 0,1 - SSE / SST,NA_real_)
  
  RMSE <- sqrt(mean((observed - predicted)^2))
  
  MAE <- mean(abs(observed - predicted))
  
  return(c(cor2 = cor2,R2_pred = R2_pred,RMSE = RMSE,MAE = MAE))
}


# ============================================================
# 8. Comparar modelo reducido y completo
# ============================================================

evaluate_qtl_outsample <- function(prepared,train,test,method = "ML") {
  
  y <- prepared$y
  K <- prepared$K
  
  # Modelo sin marcador
  X_reduced <- prepared$Xbase
  
  # Modelo con marcador
  X_full <- cbind(prepared$Xbase,prepared$Dosis)
  Dosis = prepared$Dosis
  
  reduced <- fit_predict_mixed(y = y,
    X = X_reduced,
    K = K,
    train = train,
    test = test,
    method = method, Dosis = NULL)
  
  full <- fit_predict_mixed(y = y,
    X = X_full,
    K = K,
    train = train,
    test = test,
    method = method, Dosis = Dosis)
  
  metric_reduced <- prediction_metrics(
    observed = reduced$observed,
    predicted = reduced$predicted
  )
  
  metric_full <- prediction_metrics(
    observed = full$observed,
    predicted = full$predicted
  )
  
  resultado <- data.frame(
    metric = names(metric_reduced),
    reduced = as.numeric(metric_reduced),
    full = as.numeric(metric_full),
    stringsAsFactors = FALSE
  )
  
  
  deviance <- 2 * (full$LogL - reduced$LogL)
  pval <- pchisq(q = deviance, df = full$df, lower.tail = FALSE)
  R2_M <- 1 - exp(-deviance/full$n)
  R2_MF <- 1 - (full$LogL / reduced$LogL)
  #'R_M = Pseudo-R2 de Madala (1983) y Cox and Snell (1989)
  #'R_MF = Pseudo-R2 de McFadden (1974)
  
  resultado$increment <- (resultado$full -resultado$reduced)
  
  return(list(results = resultado,
              observed = reduced$observed,
              prediction_reduced = reduced$predicted,
              prediction_full = full$predicted,
              R2_Madala= R2_M, 
              R2_McFadden = R2_MF,
              Pvalue = pval, 
              Deviance = deviance ))
}


# ============================================================
# 5. Validación cruzada repetida
# ============================================================

cross_validate_qtl <- function(prepared,
                               k = 5,
                               repeats = 50,
                               method = "ML",
                               seed = 123,
                               verbose = TRUE) {
  
  n <- length(prepared$y)
  
  folds <- make_cv_folds(n = n,k = k,repeats = repeats,seed = seed)
  
  results <- vector("list",repeats * k)
  
  contador <- 1
  
  for (r in seq_len(repeats)) {
    
    if (verbose) {message("Repetición ",r," de ",repeats)}
    
    for (f in seq_len(k)) {
      test <- which(folds[[r]] == f)
      train <- setdiff(seq_len(n),test)
      
      ans <- tryCatch(evaluate_qtl_outsample(prepared = prepared,
          train = train,test = test,method = method),
          error = function(e) {
          
          warning("Error en repetición ",r,", fold ",f,": ",conditionMessage(e))
          return(NULL)
        })
      
      if (!is.null(ans)) { 
        tmp <- ans$results
        tmp$repea <- r
        tmp$fold <- f
        tmp$n_train <- length(train)
        tmp$n_test <- length(test)
        # tmp$R2_Madala= ans$R2_Madala
        # tmp$R2_McFadden = ans$R2_McFadden
        # tmp$Pvalue = ans$Pvalue
        # tmp$Deviance = ans$Deviance 
        
        results[[contador]] <- tmp
        
        contador <- contador + 1
      }
    }
  }
  
  results <- results[
    !vapply(
      results,
      is.null,
      logical(1)
    )
  ]
  
  if (length(results) == 0) {
    stop(
      "No se completó ninguna partición. Revise los mensajes de error."
    )
  }
  
  resultado_final <- do.call(rbind,results)
  
  rownames(resultado_final) <- NULL
  
  return(resultado_final)
}



