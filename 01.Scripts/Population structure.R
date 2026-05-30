############################################################################################################
##############                 POPULATION STRUCTURE                                           ##############
##############                                                                                ##############
############################################################################################################

library(adegenet)
library(ape)
library(dendextend)

#geno correspond to genotype data only SNPs and individuals
geno <-  "//192.168.153.238/biodata10/rstudio/fernando/02.GWAS/papers_2025_sanidad/base_genotipica_121168_SNPs_519_genotipos_ACGT_lista_gwas.csv"

set.seed(31401)
dt <- (t(geno))

d3 <- df2genind(geno, ploidy=10, sep="", NA.char = NA)


grp <- find.clusters(d3, max.n=20, n.pca=50,
                     scale=F,
                     choose.n.clust = F,
                     criterion = "min")

dapc3 <- dapc(d3, grp$grp, n.pca = 12,  n.da = length(unique(grp$grp)) - 1)

############################################################################################################
##############  The following code allows you to obtain a graph of the subpopulations found,  ##############
##############  as well as showing the individuals of interest.                               ##############
############################################################################################################

dapc3_result <- dapc3$posterior
matrix <- dist(dapc3_result)

clusters <- apply(dapc3_result, 1, which.max)
names(clusters) <- rownames(dapc3_result)

hc <- hclust(matrix, method = "average")
dend <- as.dendrogram(hc)

desired_order <- names(sort(clusters))
dend <- rotate(dend, order = desired_order)

phylo_tree <- as.phylo(dend)

my_colors <- c("blue", "aquamarine4", "plum1", "steelblue1")

clusters_final <- clusters[phylo_tree$tip.label]
cluster_colors <- my_colors[clusters_final]

highlight <- highlight <- scan(text = "
  S1
S10
S100
S101
S102
S103
S104
S105
S106
S107
S108
S109
S11
S110
S111
S112
S113
S114
S115
S116
S117
S118
S119
S12
S120
S121
S122
S123
S124
S125
S126
S127
S128
S129
S13
S130
S131
S132
S133
S134
S135
S136
S137
S138
S139
S14
S140
S141
S142
S143
S144
S145
S146
S147
S148
S149
S15
S150
S151
S152
S153
S154
S155
S156
S157
S158
S159
S16
S160
S161
S162
S163
S164
S165
S166
S167
S168
S169
S17
S170
S171
S172
S173
S174
S175
S176
S177
S178
S179
S18
S180
S181
S182
S183
S184
S185
S186
S187
S188
S189
S19
S190
S191
S192
S193
S194
S195
S196
S197
S198
S199
S2
S20
S200
S201
S202
S203
S204
S205
S206
S207
S208
S209
S21
S210
S211
S212
S214
S215
S216
S217
S218
S219
S22
S220
S221
S23
S24
S25
S26
S27
S28
S29
S3
S30
S31
S32
S33
S34
S35
S36
S37
S38
S39
S4
S40
S41
S42
S43
S44
S45
S46
S47
S48
S49
S5
S50
S51
S52
S53
S54
S55
S56
S57
S58
S59
S6
S60
S61
S62
S63
S64
S65
S66
S67
S68
S69
S7
S70
S71
S72
S73
S74
S75
S76
S77
S78
S79
S8
S80
S81
S82
S83
S84
S85
S86
S87
S88
S89
S9
S90
S91
S92
S93
S94
S95
S96
S97
S98
S99
SCC10
SCC101
SCC103
SCC105
SCC106
SCC113
SCC114
SCC115
SCC116
SCC117
SCC121
SCC123
SCC124
SCC126
SCC128
SCC129
SCC13
SCC130
SCC132
SCC133
SCC136
SCC137
SCC138
SCC139
SCC140
SCC142
SCC143
SCC144
SCC145
SCC146
SCC148
SCC149
SCC15
SCC153
SCC154
SCC155
SCC156
SCC157
SCC159
SCC160
SCC161
SCC162
SCC163
SCC164
SCC165
SCC166
SCC167
SCC168
SCC171
SCC174
SCC175
SCC176
SCC177
SCC18
SCC180
SCC182
SCC183
SCC184
SCC185
SCC186
SCC19
SCC190
SCC2
SCC21
SCC22
SCC24
SCC26
SCC27
SCC275
SCC276
SCC277
SCC278
SCC279
SCC283
SCC284
SCC285
SCC287
SCC288
SCC294
SCC295
SCC296
SCC299
SCC3
SCC305
SCC32
SCC324
SCC326
SCC34
SCC35
SCC37
SCC38
SCC40
SCC41
SCC43
SCC44
SCC45
SCC50
SCC51
SCC54
SCC56
SCC57
SCC58
SCC6
SCC66
SCC67
SCC68
SCC69
SCC70
SCC71
SCC72
SCC74
SCC75
SCC76
SCC77
SCC78
SCC79
SCC8
SCC80
SCC81
SCC84
SCC88
SCC90
SCC91
SCC92
SCC93
SCC94
SCC96
SCC99
  
", what = "character")


tiff(
  filename = "population_structure.tiff",
  width = 4000,
  height = 4000,
  res = 600,         # 🔥 alta resolución (publicación)
  compression = "lzw"
)


plot.phylo(
  phylo_tree,
  type = "fan",
  use.edge.length = FALSE,
  tip.color = cluster_colors,
  edge.width = 0.2,
  edge.color = "grey60",
  cex = 0.18,
  label.offset = 0.05,
  show.tip.label = TRUE
  
)


lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

tip_positions <- which(phylo_tree$tip.label %in% highlight)

x <- lastPP$xx[tip_positions]
y <- lastPP$yy[tip_positions]


factor <- 1.04

points(
  x * factor,
  y * factor,
  pch = 20,
  col = "red",
  cex = 0.25
)


dev.off()

