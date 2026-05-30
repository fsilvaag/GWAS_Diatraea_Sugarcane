list.of.packages <- c("openxlsx", "tidyverse","data.table","sommer","gtools","AGHmatrix")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]                       
if(length(new.packages)) {install.packages(new.packages)}

invisible(lapply(list.of.packages, library, character.only = TRUE))      
gc()
rm(list = ls())
####1. Input data ####

# genotypic <- fread("E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/01.Abril2026/Raw data for GWAS/base_genotipica_121166_SNPs_519_genotipos_ACGT_lista_gwas.csv", data.table = F)
genotypic <- fread("E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/01.Abril2026/Raw data for GWAS/base_genotipica_121166_SNPs_519_genotipos_Dosis_Alternativo_lista_LD.csv", data.table = F)
Markers <- fread("E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/01.Abril2026/Raw data for GWAS/Associated_Markers_DSI_DiatraeaSpp.csv", data.table = F)
Phenotypic <- fread("E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/01.Abril2026/Raw data for GWAS/BLUPs_Discovery_and_Validation_Populations.csv", data.table = F)
# map = fread("E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/01.Abril2026/Raw data for GWAS/map_file_121166_SNPs_519_genotipos.csv", data.table = F)

####2.Data Manipulation ####
#Keep only associated markers and Validation Population
geno <- genotypic %>% filter(Name %in% Phenotypic$Name[Phenotypic$Population == "Validation"]) %>% filter(!Name %in% c("S3","S58")) %>% 
  dplyr::select(c("Name",Markers$Marker))

#Check if markers are polymorphic
sort(apply(geno[,-1], 2, function(s) length(unique(s))))

Pheno = Phenotypic[ Phenotypic$Population == "Validation" & !Phenotypic$Name %in% c("S3","S58"),]

#Removing monomorphic marker
Markers <- Markers %>% filter(Marker != "contig_46852:10466_C/T")

# Kinship
K = AGHmatrix::Gmatrix(SNPmatrix = as.matrix(genotypic %>% remove_rownames() %>% 
                                               column_to_rownames("Name")),ploidy = 10)
K1 = K[rownames(K) %in% Pheno$Name, colnames(K) %in% Pheno$Name]

####3. Convert to genetic model ####
genMar <- function(Mi, model) {
  #'Mi = Marker i
  #'model = Modelo associated with the marker i
  #'ploidy = Level of ploidy of the marker i
  
  n <- length(Mi)
  Mi <- round(Mi)
  q <- mean(Mi, na.rm = TRUE) / ploidy # Allele Frequency
  p = 1 - q
  
  # 1. Transformation according to the model
  if (grepl("dom", model, fixed = TRUE)) {
    tmp <- strsplit(model, "-", fixed = TRUE)[[1]]
    dom.order <- as.integer(tmp[1])
    Mi <- if (tmp[3] == "alt") as.numeric(Mi >= dom.order) else as.numeric(Mi > ploidy - dom.order)
    
  } else if (model %in% c("general")) {
    X <- model.matrix(~ 0+factor(Mi))
    colnames(X) <- gsub("factor|\\(|Mi|\\)","",colnames(X))
    colnames(X) <- paste0("Dosis_",colnames(X))
    return(if (ncol(X) == 0) NULL else X)
  } else if (model %in% "additive") {
    Mi = Mi
  }
  
  return(matrix(Mi))
}
ploidy = 10
res = list()
for (mk in 1:nrow(Markers)) {
  if (Markers$Model[mk] == "general") {
    dat = left_join(Pheno, geno[,c("Name", Markers$Marker[mk])])
    dat <- dat %>% bind_cols(genMar(Mi = dat[,colnames(dat) %in% Markers$Marker[mk]],
               model = Markers$Model[mk]))
    var_x <- colnames(dat)[grep("Dosis", colnames(dat))]
    Q = as.matrix(dat[,c("Q1","Q3")])
    Q1 = svd(Q)
    r <- max(which(Q1$d > 1e-08))
    Q2 = as.data.frame(as.matrix(Q1$u[, 1:r])) %>% mutate(Name = dat$Name)
    
    dat = left_join(dat,Q2)
    
    dat = dat %>%mutate(across(c(var_x,"Name"), as.factor))
    
    modelos <- lapply(var_x, function(gen) {
      form <- as.formula(paste("BLUP_DSI ~", gen, "+ V1 "))
      ft = mmes(fixed  = form,
        random = ~ vsm(ism(Name), Gu = K1),
        data   = dat,
        verbose = FALSE
      )
      wald <- wald.test(b = ft$b,
        Sigma = ft$Ci[1:nrow(ft$b), 1:nrow(ft$b)],
        Terms = 2   # The Marker
      ) 
      sm = summary(ft)$betas %>% 
        rownames_to_column("FV") %>% 
        full_join(data.frame(FV=Markers$Marker[mk],
                             chi2=wald[["result"]][["chi2"]][["chi2"]] ,
                             df=wald[["result"]][["chi2"]][["df"]], 
                             pval=wald[["result"]][["chi2"]][["P"]])) %>% 
        mutate(Fvalue = t.value^2)
      
      return(list(Modelo = form, fit = ft, Wald = wald, Summary = sm, MK = paste0(Markers$Model[mk],"--",Markers$Marker[mk])))
    })
    names(modelos) <- var_x
    
  } else if (grepl("dom", Markers$Model[mk], fixed = TRUE)) {
     
    dat = left_join(Pheno, geno[,c("Name", Markers$Marker[mk])])
    dat <- dat %>% mutate(Dom = as.numeric(genMar(Mi = dat[,colnames(dat) %in% Markers$Marker[mk]],
                       model = Markers$Model[mk])))
    Q = as.matrix(dat[,c("Q1","Q3")])
    Q1 = svd(Q)
    r <- max(which(Q1$d > 1e-08))
    Q2 = as.data.frame(as.matrix(Q1$u[, 1:r])) %>% mutate(Name = dat$Name)
    
    dat = left_join(dat,Q2)
    dat = dat %>%mutate(across(c("Name","Dom"), as.factor))
    
    
    form <- as.formula(("BLUP_DSI ~ Dom + V1 "))
    ft = mmes(fixed  = form,
                random = ~ vsm(ism(Name), Gu = K1 ),
                data   = dat,
                verbose = FALSE, dateWarning = F
      )
    
    wald <- wald.test(b = ft$b,
                        Sigma = ft$Ci[1:nrow(ft$b), 1:nrow(ft$b)],
                        Terms = 2   # The Marker
      )
      sm = summary(ft)$betas %>%
        rownames_to_column("FV") %>%
        full_join(data.frame(FV=Markers$Marker[mk],
                             chi2=wald[["result"]][["chi2"]][["chi2"]] ,
                             df=wald[["result"]][["chi2"]][["df"]],
                             pval=wald[["result"]][["chi2"]][["P"]])) %>%
        mutate(Fvalue = t.value^2)

      modelos = list(Modelo = form, fit = ft, Wald = wald, Summary = sm, MK = paste0(Markers$Model[mk],"--",Markers$Marker[mk]))
  } else if(Markers$Model[mk] == "additive") {
    dat = left_join(Pheno, geno[,c("Name", Markers$Marker[mk])])
    Q = as.matrix(dat[,c("Q1","Q3")])
    Q1 = svd(Q)
    r <- max(which(Q1$d > 1e-08))
    Q2 = as.data.frame(as.matrix(Q1$u[, 1:r])) %>% mutate(Name = dat$Name)
    
    dat = left_join(dat,Q2)
    dat = dat %>%mutate(across(c("Name"), as.factor))
    colnames(dat)[grep(Markers$Marker[mk], colnames(dat))] <- "SNP"
    
    form <- as.formula(paste0("BLUP_DSI ~ SNP + V1 "))
    ft = mmes(fixed  = form,
              random = ~ vsm(ism(Name), Gu = K1 ),
              data   = dat,
              verbose = FALSE, dateWarning = F
    )
    wald <- wald.test(b = ft$b,
                      Sigma = ft$Ci[1:nrow(ft$b), 1:nrow(ft$b)],
                      Terms = 2   # The Marker
    )
    sm = summary(ft)$betas %>%
      rownames_to_column("FV") %>%
      full_join(data.frame(FV=Markers$Marker[mk],
                           chi2=wald[["result"]][["chi2"]][["chi2"]] ,
                           df=wald[["result"]][["chi2"]][["df"]],
                           pval=wald[["result"]][["chi2"]][["P"]])) %>%
      mutate(Fvalue = t.value^2)
    
    modelos = list(Modelo = form, fit = ft, Wald = wald, Summary = sm, MK = paste0(Markers$Model[mk],"--",Markers$Marker[mk]))
    
  }
  
   res[[paste0("SNP",Markers$Marker[mk],"--",Markers$Model[mk])]] = modelos
   
  }

gen = res[grepl("--general", names(res))]
lista_plana <- unlist(gen, recursive = FALSE)
solo_summaries <- lapply(lista_plana, function(x) x$Summary)
total_summary <- rbindlist(solo_summaries, 
                           use.names = TRUE, 
                           fill = TRUE, 
                           idcol = "Origen") %>% 
  separate(Origen, into = c("Marcador", "Dosis"), sep = "\\.") %>% 
  separate(Marcador, into = c("SNP", "Modelo"), sep = "\\--")

nogen = res[!grepl("--general", names(res))]
solo_summaries <- lapply(nogen, function(x) x$Summary)

total_summary_otros <- rbindlist(solo_summaries, 
                                 use.names = TRUE, 
                                 fill = TRUE, 
                                 idcol = "Marcador") %>% 
  separate(Marcador, into = c("SNP","Modelo"), sep = "\\--") %>% 
  mutate(Dosis = NA)

All = rbind(total_summary, total_summary_otros)

fwrite(All, "E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/01.Abril2026/Result/Marker_Validation_for_Dsaccharalis_MLM.csv")
all1 <- All %>% filter(!is.na(chi2))
all1 = all1[ all1$SNP %in% all1$SNP[all1$pval <= 0.05],]
fwrite(all1, "E:/Unidades compartidas/SAM_Diatraea (articulo Socolen)/01.Abril2026/Result/Marker_Validation_for_Dsaccharalis_MLM_Validated.csv")

