list.of.packages <- c("tidyverse","data.table","sommer","gtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {install.packages(new.packages)}
invisible(lapply(list.of.packages, library, character.only = TRUE))

rm(list = ls()[!ls() %in% c("M", "path","res","ld")])
d = fread("base_genotipica_121166_SNPs_519_genotipos_Dosis_Alternativo_lista_LD.csv", data.table = F)
map = fread("map_file_121166_SNPs_519_genotipos.csv", data.table = F)

map = map[,1:3]
colnames(map) = c("Locus","LG","Position")
mm = map %>% group_by(LG) %>% summarise(N = sum(!is.na(Locus)))
# map1 = map[ !map$LG %in% mm$LG[mm$N <= 50],]
map1 = subset(map, LG %in% c("1","2","3","4","5","6","7","8","9","10")) %>%
arrange(LG, Position)


d1 = (d[,colnames(d) %in% c("Name", map1$Locus)] %>% remove_rownames() %>%
  column_to_rownames("Name"))
d1 = d1[ !rownames(d1) %in% rownames(d1)[grep("SCC",rownames(d1))],]
dim(map1)
dim(d1)

d1 <- d1[,match(colnames(d1), map1$Locus)]

res3 <- LD.decay( M, map,silent=FALSE,unlinked=T,gamma=0.95)
res <- LD.decay( d1, map1,silent=FALSE,unlinked=FALSE,gamma=0.95)
save(res, file = "/biodatas/biodata10/rstudio/fernando/LD_Decay_sommer_unlinked_FALSE_gamma_095_chromosomes_1_10_with_519_Genotypes.RData")

# load("/biodatas:/02.GWAS/01.BD_Genotipicas/LD_Decay/01.Mayo2026/01.Crom_1_10_220Genoti/LD_Decay_sommer_unlinked_FALSE_gamma_095_chromosomes_1_10.RData")

cr = gtools::mixedsort(names(res$by.LG)[-1])
ld = res
for (i in 1:length(cr)) {
  ld[["by.LG"]][[cr[i]]] <- ld[["by.LG"]][[cr[i]]][which(ld[["by.LG"]][[cr[i]]][["p"]] < 0.001),]
}

f <- function(x) {
  (dim = dim(x)[1])
}
aq = lapply(res$by.LG, f)

aq1= sapply(aq, function(x) x <= 100)
qwe = names(which(aq1 == T))
ld[["by.LG"]] <- within(ld[["by.LG"]], rm(list = qwe))
cr = gtools::mixedsort(names(ld$by.LG)[-1])

# for (i in cr) {
#   assign(paste0("modSPAN01Chr_",i),
#          loess(r2 ~ d, data=ld[["by.LG"]][[i]], span = 0.1))
# }

for (i in cr) {
  assign(paste0("modSPAN01Chr_",i),
         bam(r2 ~ s(d, bs = "cr", k =40), data=ld[["by.LG"]][[i]], span = 0.1))
}

# Savet the files
fd = mixedsort(ls()[grep("modSPAN01Chr_",ls())])
for (i in fd) {
  saveRDS(get(paste0(i)),paste0("Loess_Model_span_01_por_Cromosoma_Base_genotipica_85037_SNPs_Chr_",i,
                                              "_con_220_Genotipos.rds"))
}

# Load the loess curve models

# setwd("/biodatas/biodata10/rstudio/fernando/02.GWAS/01.BD_Genotipicas/LD_Decay/01.Mayo2026/01.Crom_1_10_220Genoti/01.LOESS_Model/")
# 
# nm = list.files()[-1]
# for(iop in 1:length(nm)){
#   assign( gsub(paste(c("Loess_Model_span_01_por_Cromosoma_Base_genotipica_85037_SNPs_Chr_modSPAN01Chr_","_con_220_Genotipos.rds"),
#                      collapse = "|"),"modSPAN01Chr_",
#                nm[iop]),
#           readRDS(nm[iop]))
# }


lk = fd = mixedsort(ls()[grep("modSPAN01Chr_",ls())])
b = 1000 # Number of bases in which the prediction jumps
target = 0.1
tolerance = 0.0025
plotlimit = 5*1e6 # Limit of the X axis to create a graph with lower size

lid = data.frame(Chromosome = numeric(), 
                 LD_decay_bp = numeric() ,
                 LD_decay_kb = numeric() ,
                 r2_threshold = numeric(),
                 r2 = numeric())

for( j in seq_along(fd)) {
  
  asd = as.data.frame(ld[["by.LG"]][[cr[j]]])
  asd <- asd[!is.na(asd$d) & !is.na(asd$r2), ]
  
  fit = get(fd[j])
    dr <- max(asd$d, na.rm = TRUE)
    xg <- seq(min(asd$d, na.rm = TRUE),dr,by = b)
  
    lilo <- predict(fit, newdata = data.frame(d = xg))
    valid <- !is.na(lilo)
    xg <- xg[valid]
    lilo <- lilo[valid]
    ix_down <- which(abs(lilo - target) <= tolerance)
    
    if (length(ix_down) > 0) {
      k <- ix_down[1]
      x0 <- xg[k]
      r2_x0 <- lilo[k]
    } else {
      x0 <- NA
    }
    
    plot_limit = dr*0.1
    
    asd10_plot <- asd %>% mutate(
      d_kb = d / 1000,
      d_Mb = d / 1e6
    ) %>%
      filter(d <= plot_limit+1e06)
    
    curve10_plot <- data.frame(
      d = xg,
      d_kb = xg / 1000,
      r2_loess = lilo
    ) %>%
      filter(!is.na(r2_loess), d <= plot_limit+1e06)
    
    x0_10 <- x0
    x0_10_kb <- x0_10 / 1000
    
    
    mi_plot <- ggplot(asd10_plot, aes(x = d_kb, y = r2)) +
      geom_point(color = alpha("grey70", 0.15), size = 0.03) +
      geom_line(data = curve10_plot,aes(x = d_kb, y = r2_loess),color = "black",
                linewidth = 1.0) +
      geom_hline(yintercept = target,color = "tomato",linetype = "solid",
                 linewidth = 1.5) +
      geom_vline(xintercept = x0_10_kb,color = "green",linewidth = 1.5) +
      annotate("text",x = x0_10_kb,y = target + 0.03,
               label = paste0(round(x0_10_kb, 1), " kb"),
               color = "black",hjust = -0.1,size = 10) +
      labs(x = "Distance (kb)",
           y = bquote("LD (" ~ r^2 ~ ")"),
           title = paste0("Chromosome ", j)) + 
      theme_bw(base_size = 16) +
      theme(panel.grid = element_blank(), 
            axis.title = element_text(face = "bold", size = 45), 
            axis.text = element_text(face = "bold", size = 25),
            plot.title = element_text(face = "bold", size = 60, hjust = 0.5)) +
      scale_x_continuous(expand = c(0,0), limits = c(0, (plot_limit/1000)+1000)) +
      scale_y_continuous(expand = c(0,0), limits = c(0,NA))
    
    ggsave(filename = paste0("LD_Decay_Span_01_Chromosome_",j,"_ST.tiff"), 
           plot = mi_plot, units = "in", width = 35, 
           height = 18, dpi = 300)
    
    rm(mi_plot);rm(asd10_plot)
    lid[j, "Chromosome"] <- cr[j]
    lid[j, "LD_decay_bp"] <- x0
    lid[j, "LD_decay_kb"] <- x0 / 1000
    lid[j, "r2"] <- r2_x0
    lid[j, "r2_threshold"] <- target
  
}

ld_decay_summary <- lid %>%
  summarise(
    mean_LD_decay_kb = mean(LD_decay_kb, na.rm = TRUE),
    sd_LD_decay_kb = sd(LD_decay_kb, na.rm = TRUE),
    min_LD_decay_kb = min(LD_decay_kb, na.rm = TRUE),
    max_LD_decay_kb = max(LD_decay_kb, na.rm = TRUE),
    median_LD_decay_kb = median(LD_decay_kb, na.rm = TRUE)
  )
ld_decay_summary

#'Across all chormosomes

ld = res
ld$all.LG <- ld$all.LG[which(ld$all.LG$p < .0001),]

library(mgcv)
mod = bam(r2 ~ s(d, bs = "cr", k =40), data=ld[["all.LG"]])

xg  <- seq(min(ldd$all.LG$d, na.rm=TRUE), 
           max(ldd$all.LG$d, na.rm=TRUE), 
           length.out = 10e06)
lilo  <- predict(mod, newdata = data.frame(d = xg))
valid <- !is.na(lilo)
xg <- xg[valid]
lilo <- lilo[valid]
ix_down <- which(abs(lilo - target) <= tolerance)

if (length(ix_down) > 0) {
  k <- ix_down[1]
  x0 <- xg[k]
  r2_x0 <- lilo[k]

for (oi in 1:10) {
  print(paste0(gsub(paste(c("modSPAN01","modSPAN01Chr_"), collapse = "|"),"",lk[oi]),
               " con dimensiones de ", dim(ld[["by.LG"]][[cr[oi]]])[1] ))
}

target = 0.1
fit = get(lk[1])
asd = as.data.frame(ld[["by.LG"]][[cr[1]]])

xg  <- seq(min(asd$d, na.rm=TRUE), max(asd$d, na.rm=TRUE), length.out = 5000)
yg  <- predict(fit, newdata = data.frame(d = xg))

ix  <- which(diff(sign(yg - target)) != 0)          # brackets where it crosses
stopifnot(length(ix) > 0)                            # no crossing? then 0.1 is outside the curve range

j   <- ix[1]                                         # choose leftmost crossing
x0  <- uniroot(function(x) predict(fit, data.frame(d = x)) - target,
               lower = xg[j], upper = xg[j+1])$root

tiff(filename = "../LD_Decay_Chr_8_Span_01_1.tiff", 
     units="in", width=20, height=12, res=300)

ggplot(asd, aes(x = d, y = r2)) +
  geom_point(color = scales::alpha("cadetblue", 0.2), size = 0.03) +
  geom_smooth(method = "loess", span = 0.1, se = FALSE, linewidth = 0.6) +
  geom_hline(yintercept = target, color = "tomato", linewidth = 0.41) +
  geom_vline(xintercept = x0, color = "tomato", linetype = 1, linewidth = 0.41) +
  annotate("point", x = x0, y = target, color = "tomato", size = 2) +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

dev.off()


## Create one graph
library(magick)
# setwd("R:/02.GWAS/01.BD_Genotipicas/LD_Decay/01.Mayo2026/01.Crom_1_10_220Genoti/03.LD_Plots/")

archivos <- mixedsort(list.files())[-c(1)]

imagenes <- image_read(archivos)
imagenes <- image_background(imagenes, color = "white", flatten = TRUE)
imagenes <- image_trim(imagenes)
imagenes_uniformes <- image_resize(imagenes, geometry = "1000x")
info <- image_info(imagenes_uniformes)

max_width  <- max(info$width)
max_height <- max(info$height)

imagenes_uniformes <- image_extent(imagenes_uniformes,
  geometry = paste0(max_width, "x", max_height),
  gravity = "center",
  color = "white"
)


rasters <- lapply(1:length(imagenes_uniformes), function(i) {
  as.raster(imagenes_uniformes[i])
})

dpi <- 600
ancho_cm <- 18
alto_cm <- 22  # Ajustable según el aspecto final

width_px <- round((ancho_cm / 2.54) * dpi)
height_px <- round((alto_cm / 2.54) * dpi)

tiff(
  filename = "Figure_S1.tiff",
  width = width_px,
  height = height_px,
  res = dpi,
  compression = "lzw",   # Compresión LZW estándar de R (sin fallos)
  type = "cairo"          # Generador de imagen Cairo súper estable
)
library(grid)
grid.newpage()

pushViewport(plotViewport(c(0.5, 0.5, 0.5, 0.5))) # Margen exterior
pushViewport(viewport(layout = grid.layout(5, 2, widths = unit(c(1, 1), "null"))))
# pushViewport(viewport(layout = grid.layout(4, 3, widths = unit(c(1, 1, 1), "null"))))


idx <- 1
for (r in 1:5) {
  for (c in 1:2) {
    if (idx <= length(rasters)) {
      pushViewport(viewport(layout.pos.row = r, layout.pos.col = c))
      grid.raster(rasters[[idx]], width = unit(0.93, "npc"), height = unit(0.93, "npc"))
      popViewport()
      idx <- idx + 1
    }
  }
}

dev.off() # Cerrar el archivo
