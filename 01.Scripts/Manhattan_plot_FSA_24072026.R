#' Create Manhattan plot
#' 
#' Create Manhattan plot
#' 
#' Results for the ref and alt versions of the dominance model are combined. If \code{data} is the output from \code{\link{set.threshold}}, then the threshold is displayed as a horizontal dashed line when \code{models} contains a single model. Because the threshold varies between models, it is not drawn when multiple models are included. Although the ref and alt versions of each dominance model are slightly different (as seen with \code{\link{qq.plot}}), they are treated as a single model for the Manhattan plot, and the average threshold is shown.
#' 
#' @param data Variable of class \code{GWASpoly.fitted}
#' @param traits Vector of trait names (by default, all traits plotted)
#' @param models Vector of model names (by default, all models plotted)
#' @param chrom.hold Vector of chromosomes to hold in the analysis. This is intended when organism  has alpha-numeric chromosomes (by default it is NULL)
#' @param solid a boolean to indicate if the points of the manhatan should be solid or not (Default TRUE)
#' @param FDR a boolean to indicate if Manhattan plot will follow the FDR Threshold or not
#' @param map Archivo con la informacion de los marcadores, cromosoma, posicion, y demas.
#'            Para archivos de GWASpoly, este se toma del archivo, para los demas hay que proveerlo.
#' @param umb2 Segundo Umbral para plotear
#' @param nme a name to plot in the graph when there are contigs or scaffolds.
#'            if null the names of the genotypic data are used as chromosomes
#'            else the plot will use this name 
#' @param mkr Listado de nombres de los marcadores que se quieren plotear en el Manhattan
#' @param plot Boleano para indicar si se quiere graficar todos los modelos en la misma grafica (TRUE)
#'        o si se quiere graficar un solo modelo por grafica (FALSE). Por defecto es FALSE
#' @param xtes Tamano de letra del eje x de la grafica. Por defecto esta en 15        
#' @param hj Es el valor del hjust del grafico en el plot title de ggplos. Por defecto es 0.5 x
#' @param lp la posicion de la leyenda en el ggplot. Por defect es "none"
#' @param at Es el tamano de la letra de los titulos del eje x. Por defecto es 15
#' @param axt Es el tamano de la letra de los ejes del grafico. Por defecto es 14
#' @param scx Titulo del eje X del grafico. POr defecto es "Chromosome"
#' @param scm Es el listado de colores para los puntos de la grafica. Por defecto son
#' c("#EC5f67", "#FAC863",'#99C794','#6699CC','#C594C5'
#' @param ma Es el maximo valor del eje Y para graficar en el ggplot. Por defecto es nulo
#'            para que el mismo programa escoja el maximo. Se puede especificar para 
#'            que todos los graficos tengan el mismo maximo 
#' @param mi Es el minimo valor del eje Y para graficar en el ggplot. Por defecto es nulo
#'            para que el mismo programa escoja el maximo. Se puede especificar para 
#'            que todos los graficos tengan el mismo maximo 
#' @param us Tamano o grosor de la linea del umbral. Por defecto es 1
#' @param uc Color de la linea de umbral. Por defecto es "darkgreen"
#' @param NOMBRE Nombre que se le quiere asignar a la combinacion de los contigs y scaffolds
#'               para las graficas
#' @param nombre_SNP Nombre que queremos a priori darle a los SNPs que queremos resaltar en el grafico
#' @param threshol un boleano de TRUE para cuando se quiere utilizar el threshold de GWASpoly
#'                  o FALSE para cuando se quiere utilizar el FDR teorico.
#' 
#' 
#' @return ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#

get_x <- function(map, chrom, position){
  
  asd = as.factor( gtools::mixedsort(unique(map[,chrom])))
  colnames(map)[chrom] <- "Chr"
  map <- map %>% 
    mutate(Chr = factor(Chr), 
           Chr = as.numeric(fct_relevel(Chr,paste0(asd)))) %>% 
    arrange(Chr)
  a <- tapply(as.numeric(paste0(map[,position])), map[,"Chr"], max)
  n <- length(a)
  m <- tapply(map[,position], map[,"Chr"], length)
  b <- c(0, cumsum(a)[-n])
  # b <- c(0,apply(array(1:(n-1)),1,function(k) {
  #   sum(a[1:k])
  #   }))
  map$x <- map[,position] + rep(b,times=m)
  return(map)
}



manhattan.plot_FSA <- function(data,traits=NULL,models=NULL, 
                               chrom.hold = NULL, solid= TRUE, 
                               umbral = 0.05, 
                               nme = NULL,map = NULL,
                               mkr = NULL,
                               plot = F, xtes= 15, hj = 0.5, lp = "none",
                               at = 15, axt = 14, scx = "Chromosome",
                               scm = c("#EC5f67", "#FAC863",'#99C794','#6699CC','#C594C5'),
                               ma = NULL, mi = NULL, us = 1, uc = "darkgreen",
                               snpsize = 3, snppadd = 0.6, bpad = 0.6, segcol = "grey25",nombre_SNP=NULL){
  
  render_plot <- function(df, mark, md) {
    labels_x <- gtools::mixedsort(unique(df$Chrom))
    if (!is.null(nme)) {
      labels_x <- c(labels_x[-max(as.numeric(unique(df$Chrom)))], nme)
    }
    breaks_x <- (tapply(df$x, df$Chrom, max) + tapply(df$x, df$Chrom, min)) / 2
    y_max <- if(!is.null(ma)) ma else max(df$Q, na.rm = TRUE) + 0.5
    
    if(is.null(mi)){
      mi <- min(asd$Q, na.rm = T)
    }
    
    ural = data@threshold[1,colnames(data@threshold) %in% md]
    
    p <- ggplot(df, aes(x=x, y=Q, colour=factor(Chrom), shape=model)) +
      geom_point() +
      geom_hline(yintercept = ural, 
                 col = uc, linewidth = us) +
      scale_x_continuous(expand = c(0,0), 
                         name=scx, 
                         breaks=breaks_x, 
                         labels=labels_x) +
      scale_y_continuous(limits = c(0, y_max), 
                         expand = c(0,0)) +
      scale_color_manual(values = rep(scm, length.out = length(labels_x))) +
      labs(y = expression(bold(paste("-log"[10],"(p)")))) +
      theme_bw() + 
      theme(text = element_text(size=xtes), 
            panel.grid = element_blank(),
            plot.title = element_text(hjust = hj), 
            legend.position = lp,
            axis.title = element_text(face = "bold", size = at),
            axis.text = element_text(size = axt)) +
      {if (!is.null(mark)) ggrepel::geom_text_repel(data = mark,  
                                                    aes(x = x, y = score, label = Name),
                                                    inherit.aes = FALSE,
                                                    size = snpsize,
                                                    fontface = "bold",
                                                    box.padding = bpad,
                                                    point.padding = snppadd,
                                                    min.segment.length = 0,
                                                    segment.color = segcol,
                                                    max.overlaps = Inf,
                                                    seed = 123)
        }  +
      guides(colour = "none") 
    
    return(p)
  }
  
  if(inherits(data, "GWASpoly.fitted")) {
    map1 <- data@map
    traits <- traits %||% names(scores)
    all_models <- colnames(data@scores[[1]])
  } else {
    if(is.null(map)) stop("Error: Si 'data' no es objeto GWASpoly, se debe incluir el objeto 'map'.")
    map_base <- map
    traits <- traits %||% names(scores_all)
    all_models <- names(scores_all[[traits[1]]])
  }
  
  n.trait <- length(traits)
  
  if (is.null(models)) {
    stop("Model to be plotted has to be specified")
  } else {
    models <- unlist(lapply(as.list(models),function(x){all_models[grep(x,all_models,fixed=T)]}))
    stopifnot(all(is.element(models,all_models)))
  }
  n.model <- length(models)
  
  if (!is.factor(map1[,2])) {
    map1[,2] <- as.factor(map1[,2])
  }
  
  if (is.null(chrom.hold)) {
    cat("There are Contig, Scaffolds or alpha-numeric coded chromosomes. All of them will be coded as c+1, for c = last chromosome\n")
    chr <- droplevels(unique(map1[,2])[!unique(map1[,2]) %in% c(paste0(1:1000))])
    chr1 <- as.numeric(droplevels(unique(map1[,2])[unique(map1[,2]) %in% c(paste0(1:1000))]))
    
    map1 <- transform(map1,Chr = with(map1,ifelse(map1[,2] %in% chr,as.character(max(chr1)+1),as.character(paste0(map1[,2])))))
    
    map1 <- map1[order(map1[,"Chrom"],map1[,"Position"]),]
    qasd <- nrow(subset(map1, !map1$Chr %in% chr1))
    map1$Pos <- ifelse(!map1$Chr %in% chr1, seq(from = 1, to =1000*qasd, by =1000), map1$Position)
    map_plotme <- transform(map1, Chr = as.numeric(Chr))
  } else { 
    cat("The user specified the chromosomes to plot\n")
    map1 <- droplevels(map1 %>% filter(map1[,2] %in% chrom.hold))
    map1$Chr <- as.numeric(map1[,2])
    map1$Pos = map1[,3]
    map_plotme <- map1
  }
  
  #### Crear el archivo para graficar el qqplot
  plotme <- NULL
  pos = grep(paste(c("^Pos$"), collapse = "|"), colnames(map_plotme))
  crom = grep(paste(c("^Chr$"), collapse = "|"), colnames(map_plotme))
  x <- get_x(map = map_plotme, chrom = crom, position = pos)
  
  for (k in 1:n.trait) {
    ###
    asd = as.factor( gtools::mixedsort(unique(map_plotme[,"Chr"])))
    scores <- as.data.frame(data@scores[[traits[k]]]) %>% 
      rownames_to_column("Marker") %>%
      filter(Marker %in% map_plotme$Marker) %>% 
      mutate(Chrom = factor(map_plotme[,"Chr"]), Chrom = fct_relevel(Chrom,paste0(asd))) %>% 
      arrange(Chrom) %>% left_join(x[,c("Marker","x")], by = "Marker") %>% 
      mutate(color = factor(ifelse(as.integer(factor(map_plotme[,2]))%%2==1,1,0)))
    
    tmp <- as.data.frame(pivot_longer(data=scores,cols=match(models,colnames(scores)),
                                      names_to="model",values_to="score", 
                                      values_drop_na=TRUE)) %>%
      dplyr::select(c("Marker","Chrom","x","color","model","score"))
    tmp$trait <- traits[k]
    plotme <- bind_rows(plotme,tmp)
    rm(asd)
  }
  
  plotme <- transform(plotme, model = factor(model),Chr = as.numeric(Chrom),
                      trait = factor(trait))
  
  if (!is.null(mkr)) {
    mmmmm = plotme %>% filter(Marker %in% mkr) %>% 
      dplyr::select(c("Marker","Chr","x","score"))
      
    if (is.null(nombre_SNP)) {
      mmmmm <- mmmmm  %>% mutate(Name = paste0("CC-SNP", seq(1,nrow(.))))
    } else {
          mmmmm <- mmmmm |> mutate(Name = nombre_SNP)
        }
             
  } else {
    mmmmm = NULL
  }
  
  
  Manhatan <- list()
  for (jkq in 1:n.trait) {
    tr = traits[jkq]
    for (qwe in 1:n.model) {
      asd <- droplevels(subset(plotme, plotme$trait == tr & model == models[qwe]))
      asd$Q <- asd$score
      
      p <-  render_plot(df = asd, mark = mmmmm, md = models[qwe])
      Manhatan[[paste0(traits[jkq])]][[paste0(gsub("\\-","_",models[qwe]))]] <- p
    }
  }
  return(Manhatan)
}

# data = Diatraea
# map = NULL
# traits = "DSI_BLUP_general"
# models = "general"
# chrom.hold = NULL; solid = TRUE
# umbral = 0.05; nme = "CS";umb2 = NULL
# plot = F; xtes= 15; hj = 0.5; lp = "none"
# at = 15; axt = 14; scx = "Chromosome"
# scm = c("#EC5f67", "#FAC863",'#99C794','#6699CC','#C594C5')
# ma = NULL; mi = NULL; us = 1; uc = "darkgreen"
# snpsize = 3; snppadd = 0.2
# NOMBRE = "CS"
# mkr = c("10:7656049_T/C","4:31216826_T/C","9:9208699_T/C")
# nombre_SNP = c("SNP2","SNP3","SNP4")
# snpsize = 3; snppadd = 0.6; bpad = 0.6; segcol = "grey25"
