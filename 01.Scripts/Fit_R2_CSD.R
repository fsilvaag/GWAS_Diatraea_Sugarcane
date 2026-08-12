#'Funcion para calcular el Pseudo-R2 de Madala (1983), Cox and Snell (1989) como 
#'propuesto en el articulo Giselmar et al. (2018) [fsa_656]
#'
#'@data = base datos en formato GWASpoly.thresh
#'@trait = Trait que se quiere analizar
#'@qtl = Dataframe obtenida del get.QTL de GWASpoly
#'@fixed = dataframe con la coleccion de efectos fijos
#'
#'@Resultados
#'R_M = Pseudo-R2 de Madala (1983) y Cox and Snell (1989)
#'R_MF = Pseudo-R2 de McFadden (1974)



.make.full <-
  function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
library(rrBLUP)

fit.QTL_CSD = function (data, trait, qtl, fixed = NULL) 
{
  stopifnot(inherits(data, "GWASpoly.K"))
  stopifnot(is.element(trait, names(data@scores)))
  stopifnot(qtl$Model %in% c("additive", "general", paste(1:(data@ploidy/2), 
                                                          "dom-ref", sep = "-"), paste(1:(data@ploidy/2), "dom-alt", 
                                                                                       sep = "-")))
  stopifnot(qtl$Marker %in% data@map$Marker)
  not.miss <- which(!is.na(data@pheno[, trait]))
  y <- data@pheno[not.miss, trait]
  pheno.gid <- data@pheno[not.miss, 1]
  geno.gid <- rownames(data@geno)
  n.gid <- length(geno.gid)
  n <- length(y)
  Z <- matrix(0, n, n.gid)
  Z[cbind(1:n, match(pheno.gid, geno.gid))] <- 1
  X <- matrix(1, n, 1)
  if (!is.null(fixed)) {
    for (i in 1:nrow(fixed)) {
      if (fixed$Type[i] == "factor") {
        xx <- factor(data@fixed[not.miss, fixed$Effect[i]])
        if (length(levels(xx)) > 1) {
          X <- cbind(X, model.matrix(~x, data.frame(x = xx))[, 
                                                             -1])
        }
      }
      else {
        X <- cbind(X, data@fixed[not.miss, fixed$Effect[i]])
      }
    }
  }
  chrom <- as.character(data@map$Chrom[match(qtl$Marker, data@map$Marker)])
  if (length(data@K) > 1) {
    if (length(setdiff(levels(data@map$Chrom), chrom)) == 
        0) {
      stop("LOCO model cannot be used because there are QTL on every chromosome. Run set.K with LOCO=FALSE")
    }
    K <- makeLOCO(data@K, exclude = match(chrom, levels(data@map$Chrom)))
  }
  else {
    K <- data@K[[1]]
  }
  n.qtl <- nrow(qtl)
  S <- vector("list", n.qtl)
  df <- integer(n.qtl)
  X0 <- X
  
  reduced.model <- mixed.solve(y = y, X = .make.full(X0), Z = Z, 
                            K = K, method = "ML")
  pval <- R2_M <- R2_MF <- numeric(n.qtl)
  for (i in 1:n.qtl) {
    X1 <- X0
    S[[i]] <- .design.score(data@geno[, qtl$Marker[i]], 
                            model = qtl$Model[i], ploidy = data@ploidy, min.MAF = 0, 
                            max.geno.freq = 1)
    stopifnot(!is.null(S[[i]]))
    df[i] <- ncol(S[[i]])
    X1 <- cbind(X1, Z %*% S[[i]])
    
    
    #if (n.qtl > 1) {
    #  for (j in setdiff(1:n.qtl, i)) {
    ##    X <- cbind(X, Z %*% S[[j]])
    #  }
    #}
    full.model <- mixed.solve(y = y, X = .make.full(X1), 
                                 Z = Z, K = K, method = "ML")
    deviance <- 2 * (full.model$LL - reduced.model$LL)
    pval[i] <- pchisq(q = deviance, df = df[i], lower.tail = FALSE)
    R2_M[i] <- 1 - exp(-deviance/n)
    R2_MF[i] <- 1 - (full.model$LL / reduced.model$LL)
  }
  return(data.frame(data@map[match(qtl$Marker, data@map$Marker), 
                             1:3], Model = qtl$Model, R2_M = R2_M, R2_MF = R2_MF, pval = pval))
}
