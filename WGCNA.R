library(DESeq2)
library(WGCNA)
library(GO.db)
library(topGO)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyr)
library(dplyr)
library(ggplot2)

datos <- read.csv("TumoresCaracteristicas.csv", header =TRUE, sep=",")
countData <- read.table(file = 'genes_WGCNA.txt', header = TRUE, row.names = 1)


dim(countData)  # Ver tamaño de la tabla después de filtrar

datos$extraccion <- factor(datos$extraccion)
datos$Tratamiento <- factor(datos$Tratamiento)


dds_wgcna <- DESeqDataSetFromMatrix(
   countData = countData,
   colData = datos,
   design = ~Grupo)
  
dds_wgcna <- dds_wgcna[rowSums(counts(dds_wgcna)) > 10, ]
  
#Transformación de cuentras por Estabilización de la Varianza 
vsd_wgcna <- varianceStabilizingTransformation(dds_wgcna, blind=FALSE)
class (vsd_wgcna)
#Extracción de las cuentas del objeto DSQeq  
vsd_wgcna <- assay(vsd_wgcna)
expression.data <- t(vsd_wgcna)
gsg <-goodSamplesGenes(expression.data)
gsg$allOK

#Remoción de los genes que no pasen el filtrado para realizar la red 
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(expression.data)[!gsg$goodGenes], collapse = ", "))); #Identificar los genes outliers 
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", "))); #Identificar muestras que sean outliers 
  expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Remoción de los genes y muestrasque sean outliers 
}

#Clusterización de las muestras basado en la distancia entre la expresión de
#los genes para ver agrupamiento
sampleTree <- hclust(dist(expression.data), method = "average") 

par(cex = 0.6);
par(mar = c(0,4,2,0))
pdf(file = "./cluster_samples_all_lncmRNA.pdf",
    width = 12 , height = 10)
#Plotting the cluster dendrogram
plot(sampleTree, main = "Agrupamiento de muestras", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

#Remoción de muestras que no se agrupen a las muestras y qeu no pasen el valor de corte
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 260, minSize = 1) 
expression.data <- expression.data[cut.sampleTree==1, ]
#Verificar las muestras que quedaron 
rownames(expression.data)


#Construcción de la red de correlación

datExpr <- expression.data
powers <- c(c(1:10), seq(from = 12, to=30, by=2))
#Analiis de la tipología de la red para selección de la potencia que miniminice la conectividad 
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 4,  networkType = "unsigned", corFnc = bicor)


#save(sft, file = "sft_unsigned_all.RData")
#save.image("mi_sesion_wgcna_unisgned_all_lncmRNA.RData")

#Grafica para ver la conectividad y seleecionar la potencia a utilizar 
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
pdf(file = "./scale_and_mean_connectivity_all_unsigned.pdf",
    width = 12 , height = 10)
par(mfrow = c(1,2));
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# 2.2 Seleccionar el soft power
rm(softPower)
softPower <- if (!is.na(sft$powerEstimate)) {
  sft$powerEstimate
} else {
  which.max(sft$fitIndices[,2] > 0.8)  # Primera potencia donde R² > 0.8
}

softPower

# Construcción de la red de coexpresión génica con WGCNA


net_unmer <- blockwiseModules(datExpr, power = softPower,
                              TOMType = "unsigned", minModuleSize = 50,
                              networkType = "unsigned",
                              maxBlockSize = 12000,
                              mergeCutHeight = FALSE, deepSplit = 2,
                              corType = "bicor",
                              numericLabels = FALSE, pamRespectsDendro = FALSE,
                              saveTOMs = TRUE,
                              saveTOMFileBase = "TOM_wgcna_unmerged_unsigned", 
                              verbose = 4)


# CONVERSIÓN DE etiquetas a colores para graficas 
table(net_unmer$colors)
unmergedColors <- labels2colors(net_unmer$colors)
unmergedColors

unique(unmergedColors)

# Grafica de dendrograma y del color de los modulos sin unir 
pdf(file = "./clusterDendogram_all_lncmRNA.pdf",
    width = 20 , height = 10);
plotDendroAndColors(net_unmer$dendrograms[[1]], unmergedColors[net_unmer$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#save(net_unmer, file = "WGCNA_net_unmerged_unsigned_all_lncmRNA.RData")
#write.table(net_unmer$colors, file = "WGCNA_modules_unmerged_unsigned_all_lncmRNA.txt", sep = "\t", quote = FALSE, col.names = NA)

#save.image("mi_sesion_wgcna_unisgned_all_lncmRNA.RData")

#load("WGCNA_net_unmerged.RData")


moduleLabels_unmer <- net_unmer$colors

# 1. Calcular los eigengenes de los módulos detectados sin fusionar
ME_unmerged <- moduleEigengenes(datExpr, colors = moduleLabels_unmer, softPower = softPower)
MEs <- ME_unmerged$eigengenes
head(MEs)
#Calculo de la de la disimilaridad de los eingengenes 
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")

pdf(file = "./clusterDendogram_Module_unsigned_all_lncmRNA.pdf",width = 20 , height = 10);
par(mar = c(0,4,2,0)) 
par(cex = 0.6);
plot(METree)
abline(h=.25, col = "red") #Marcar correlación del 0.75 entre modulos 
dev.off()

# 2. Fusionar los módulos cuyos eigengenes estén altamente correlacionados
merged_clusters_unmerged <- mergeCloseModules(datExpr,
                                              net_unmer$colors,
                                              cutHeight = 0.25,  
                                              MEs = ME_unmerged$eigengenes,
                                              verbose = 3)

# 3. Obtener las nuevas etiquetas y colores después de la fusión
merged_labels <- merged_clusters_unmerged$colors
merged_colors <- labels2colors(merged_labels)
merged_Es <- merged_clusters_unmerged$newMEs

ME_merged <- moduleEigengenes(datExpr, colors =merged_clusters_unmerged$colors,softPower = softPower)
#save(module_colors_merged, ME_merged, file = "WGCNA_merged_network_unsigned_all_lncmRNA.RData")

#Visualizar el dendrograma con los módulos fusionados

pdf(file = "./clusterDendogram_MERGED_unsigned_all.pdf", width = 20 , height = 10);
plotDendroAndColors(net_unmer$dendrograms[[1]], 
                    cbind(unmergedColors[net_unmer$blockGenes[[1]]], merged_labels[net_unmer$blockGenes[[1]]]),
                    groupLabels = c("Original", "Merged"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

dev.off()


# Guardar los colores de los módulos fusionados
module_colors_merged <- merged_clusters_unmerged$colors
# Guardar los eigengenes de los módulos fusionados
MEs_merged <- merged_clusters_unmerged$newMEs
#save.image("mi_sesion_wgcna_unisgned_all_lncmRNA.RData")


#Identificación de los terminos GO asociados a cada modulo 

# Extraer genes asociados a términos GO
geneID2GO <- AnnotationDbi::as.list(org.Hs.egGO)
probe2GO <- lapply(geneID2GO, names)

GOresults_lncmodules<-data.frame()
rm(enrich)
universe <- colnames(datExpr)

#for (module in unique(module_colors_merged[module_colors_merged != "grey"])) 
  for(module in unique(net_unmer$colors)){
    print(module)
    genes <- universe[module_colors_merged==module]
    # head(print(genes))
    # genes<-universe[net_mergeado$colors=="14"]
    geneList <- factor(as.integer(universe %in% genes))
    
    names(geneList) <- universe
  
    # Convierte los símbolos de genes a Entrez IDs
    gene_symbols <- names(geneList)  # Asegúrate de que geneList tiene nombres
    
    entrez_ids <- mapIds(org.Hs.eg.db, 
                         keys = gene_symbols, 
                         column = "ENTREZID", 
                         keytype = "SYMBOL", 
                         multiVals = "first")
    
    # Filtra genes sin ID asignado
    valid_ids <- !is.na(entrez_ids)
    geneList <- geneList[valid_ids]
    names(geneList) <- entrez_ids[valid_ids]
    #head(print(geneList))
    # topGO analysis
    if(exists("enrich")) {remove(enrich)}
    for(on in c("BP"))
    {
      print(on)
      # Make topGO object
      GOdata_wgcna <- new("topGOdata", 
                          ontology = "BP", 
                          allGenes = geneList, 
                          nodeSize = 5, 
                          annot = annFUN.gene2GO, 
                          gene2GO = probe2GO, 
                          geneSelectionFun = function(x) x == 1)  # Selecciona genes significativos    
      # fisher test
      result <- runTest(GOdata_wgcna, algorithm = "classic", statistic = "fisher")
      results.table <- GenTable(GOdata_wgcna, result, topNodes = 100)
      colnames(results.table) <- c("GO.ID", "Term", "Annotated", "Significant", "Expected", "result1")
      results.table$result1 <- gsub('<','', results.table$result1)
      results.table <- as.data.frame(lapply(results.table, function(x) gsub("<", "", x)))
      results.table$result1 <- as.numeric(results.table$result1)
      # Corrección del p-value por multiples pruebas utilizando Benjamini-Hochberg
      results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")
      results.table$ontology<-on
      keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.05,]
      if(exists("enrich")) enrich<- rbind(enrich, keep)
      if(!exists("enrich")) enrich<- keep
      
    }
    if(dim(enrich)[1]>0)
    {
      enrichME<-enrich
      enrichME$MEs_mer=module
      GOresults_lncmodules<-rbind(GOresults_lncmodules,enrichME)   }
    
    
write.table( GOresults_lncmodules,
file="./consensus_GOresults_mergedmodules_unsigned_onlylncmRNA_all_lncmRNA.txt",
sep="\t", row.names=FALSE)
    
#save.image("mi_sesion_wgcna_unisgnedi_all_lncmRNA.RData")


#Preparación de Dataframe con información del numero de CPG y lncRNA por modulo 
colors_dat <- as.data.frame(merged_colors)
numbers_dat <- as.data.frame(merged_labels)
all_modules <- cbind (numbers_dat, colors_dat)
colnames(all_modules) <- c("module_number", "module_color")


lnc_names <- read.table("./lnc_all.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(lnc_names) <- "names"
PCG_names <- read.table("../gencode/protein_cosing_names.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(PCG_names) <- "names"

all_modules <- rownames_to_column(all_modules, var="names")
lnc_modules <- merge(all_modules, lnc_names, by= "names", all = FALSE)
write.table(lnc_modules , file="./lnc_modules_TNC.txt",sep="\t", col.names = TRUE, row.names=FALSE)

pcg_modules <- merge(all_modules, PCG_names, by= "names", all = FALSE)

all_freq <- all_modules %>% 
  group_by(merged_labels) %>%
  summarise(genes = n())

pcg_freq <- pcg_modules%>%                   
  group_by(merged_labels) %>%                                
  summarise(genes = n())

module_lnc_freq <- lnc_modules %>%
  group_by(merged_labels) %>%     
  summarise(freq_lnc = n())

# Asegurarse de que las columnas tengan el mismo nombre
colnames(module_lnc_freq) <- c("module", "freq_lnc")

lnc_prop_df <- full_join(all_freq, pcg_freq, by = "module", suffix = c("_all", "_pcg")) %>%
  full_join(module_lnc_freq, by = "module")



# Calcular porcentajes de PCG y lncRNA en cada modulo 
lnc_prop_df$perc_lnc <- (lnc_prop_df$freq_lnc / lnc_prop_df$genes_all) * 100
lnc_prop_df$perc_codegene <- (lnc_prop_df$genes_pcg / lnc_prop_df$genes_all) * 100

lnc_prop_df$perc_otros <- 100 - lnc_prop_df$perc_lnc - lnc_prop_df$perc_codegene

lnc_prop_df$new_module <- rownames(lnc_prop_df)


plot_df <- lnc_prop_df %>%
  select(module,genes_all, perc_lnc, perc_codegene) %>%
  pivot_longer(cols = c(perc_lnc, perc_codegene,genes_all),
               names_to = "type",
               values_to = "percentage") %>%
  mutate(type = factor(type, levels = c("genes_all","perc_lnc", "perc_codegene"),
                       labels = c("genes","lncRNA", "PCG")))

plot_df <- plot_df %>%
  mutate(new_module = factor(module, levels = sort(unique(as.numeric(as.character(module_number)))))
  )

numeros <- lnc_prop_df %>%
  select(module, genes_pcg, freq_lnc) %>%
  pivot_longer(cols = c("genes_pcg", "freq_lnc"),
               names_to = "type",
               values_to = "cantidad") %>%
  mutate(
    type = factor(type, levels = c("genes_pcg", "freq_lnc"),
                  labels = c("PCG", "lncRNA")),
    module = factor(module, levels = orden_modulos) 
  )

pdf(file = "./porcentajelncRNA_codgenes_bymodule_unsigned_alllncmRNA.pdf",width = 12 , height = 10)
ggplot(plot_df , aes(x = module, y = percentage, fill = type)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            angle = 90, size = 4, color = "white") +
  scale_fill_manual(values = c("lncRNA" = "#1B9E77", "PCG" = "#7570B3")) +
  labs(title = "Proporción de Biotipos por módulo",
       x = "Módulo",
       y = "Porcentaje (%)",
       fill = "Biotipo") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_blank(),          # oculta los textos
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  # Agrega recuadros de color simulando el color del módulo debajo del gráfico
  annotate("rect",
           xmin = seq_along(levels(plot_df$module)) - 0.35,
           xmax = seq_along(levels(plot_df $module)) + 0.35,
           ymin = -5, ymax = 0,
           fill = module_colors[levels(plot_df$module)],
           color = "black")
dev.off()

write.table(plot_df, file = "porcentajes_modulos_unisgned.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lnc_prop_df, file = "freq_modulos_unsigned_TNC.txt", sep = "\t", row.names = FALSE, quote = FALSE)
lnc_prop_df
#plot_df <- read.table("porcentajes_lnc_cod.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
plot_df
ME_merged <- moduleEigengenes(datExpr,merged_labels)$eigengenes;
relacion_mod<- lnc_prop_df[,1]


nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
ME_merged <- moduleEigengenes(datExpr, colors =merged_clusters_unmerged$colors)$eigengenes


MEs = orderMEs(ME_merged)

#save.image("mi_sesion_wgcna_unisgned_lncmrna.RData")




#Calculo de la correlación modulo-caracteristica 
# Asegúrate de que 'Grupo' sea un factor
datos$Grupo <- as.factor(datos$Grupo)

# Crear variables dummy utilizando model.matrix()
dummy_vars <- model.matrix(~ Grupo - 1, data = datos)

# Combinar las variables dummy con las muestras
datTraits <- cbind(Muestra = datos$Muestra, dummy_vars)
datTraits <- as.data.frame(datTraits)
rownames(datTraits) <- NULL
datTraits <- column_to_rownames(datTraits, var = "Muestra")


# Si datTraits no es un data frame, conviértelo
datTraits <- as.data.frame(lapply(datTraits, as.numeric))


# Ahora puedes acceder a la columna 'Muestra'
rownames(datTraits) <- datos$Muestra
datTraits$Muestra <- NULL



tratamiento_bin <- binarizeCategoricalVariable(
  datos$Tratamiento,
  includePairwise = FALSE,
  includeLevelVsAll = TRUE,
  minCount = 1
)

traits2 <- data.frame(Muestra = datos$Muestra, tratamiento_bin)
rownames(traits2) <-datos$Muestra
traits2$Muestra <- NULL
colnames(traits2) <- c("Control","Yodo","Placebo")
Traits <- merge(traits2,datTraits,by=0)
Traits <- column_to_rownames(Traits, var="Row.names")

moduleTraitCor = cor(MEs, Traits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

par(mar =c(6, 8.5, 3, 3))
png("./ME_caracteristicas_unsigned_lncRNAs.png", width = 20, height = 12, units = "in", res = 300)
par(mar =c(8, 15, 3, 3))
#sizeGrWindow(8,4)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(Traits),
               yLabels = rownames(moduleTraitCor),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.7,
               zlim = c(-1,1),
               main = "Relación módulos-caracteristicas",
               legendLabel = "Correlación",  # Título en la barra de color
               cex.legendLabel = 1.6,
               cex.lab.y = 1.8,
               cex.lab.x = 1.8,         # ← Tamaño del texto en eje Y
               cex.main = 2, )

dev.off()


SubLA <- as.data.frame(Traits$GrupoLA)
names(SubLA) = "LA"
rownames(SubLA) <- rownames(Traits)
SubLA <- SubLA[rownames(datExpr), , drop = FALSE]
LA <- SubLA


SubTNC <- as.data.frame(Traits$GrupoTNC)
names(SubTNC) = "TNC"
rownames(SubTNC) <- rownames(Traits)
SubTNC <- SubTNC[rownames(datExpr), , drop = FALSE]

SubTNCY <- as.data.frame(Traits$GrupoTNC_Y)
names(SubTNCY) = "TNC"
rownames(SubTNCY) <- rownames(Traits)
SubTNCY <- SubTNCY[rownames(datExpr), , drop = FALSE]


SubLAY <- as.data.frame(Traits$GrupoLA_Y)
names(SubLAY) = "LAY"
rownames(SubLAY) <- rownames(Traits)

Pla <- as.data.frame(Traits$Placebo)
names(Pla) = "Pla"
rownames(Pla) <- rownames(Traits)
Pla <- Pla[rownames(datExpr), , drop = FALSE]  # Reordena Pla para que siga el orden de datExpr

Yod <- as.data.frame(Traits$Yodo)
names(Yod) = "Yod"
rownames(Yod) <- rownames(Traits)
Yod <- Yod[rownames(datExpr), , drop = FALSE]

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


geneTraitSignificance = as.data.frame(cor(datExpr,Yod, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))


names(geneTraitSignificance) = paste("GS.", names(Yod), sep="")
names(GSPvalue) = paste("p.GS.", names(Yod), sep="")


# Número del módulo como número, no texto
module <- "3"

# Obtener el nombre del eigengene a partir del número (como string)

# Obtener la columna correspondiente en MEs
column <- which(modNames == module)

# Filtrar genes que pertenecen a este módulo
moduleGenes <- moduleColors == module


png("./modulo3_TNC.png", width = 10, height = 10, units = "in", res = 600)

verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance[moduleGenes, 1],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Significancia de genes relacionados con el grupo triple negativo",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main =2.2, cex.lab = 1.8, cex.axis = 1.2, col = module)
dev.off()

#Analisis de pertenecia de los genes a ese modulo 

# Filtrar solo los genes del módulo 7
moduleGenes <- moduleColors == module

# Subset de geneTraitSignificance para esos genes
gs_mod7 <- geneTraitSignificance[moduleGenes, , drop = FALSE]

# Ordenar por la correlación absoluta con el rasgo Yod
idx <- order(abs(gs_mod7$GS.TNC), decreasing = TRUE)

# Obtener los 10 genes más correlacionados
top10_genes_mod7 <- rownames(gs_mod7)[idx[1:10]]

# Mostrar lista para copiar en tesis
cat(paste(top10_genes_mod7, collapse = ", "))

analisis_modulo("6", SubTNC)

NormTrait <- as.data.frame(Traits$GrupoNORM)
names(NormTrait) = "NORM"
rownames(NormTrait) <- rownames(Traits)
geneTraitSignificance = as.data.frame(cor(datExpr, NormTrait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(NormTrait), sep="")
names(GSPvalue) = paste("p.GS.", names(NormTrait), sep="")

png("./modulo34.png", width = 8, height = 10, units = "in", res = 300)
verboseScatterplot((geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in 17 module"),
                   ylab = "Significancia de genes Normal",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
