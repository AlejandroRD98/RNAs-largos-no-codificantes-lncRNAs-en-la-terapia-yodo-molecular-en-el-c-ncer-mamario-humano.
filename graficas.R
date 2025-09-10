library(ggplot2)
library(ComplexHeatmap)
library(colorRamp2)
library(VennDiagram)
library(tidyr)
library(dplyr)
library(readr)

#Grafica de barras de la cantidad de lncRNAs y PCG por grupo (Figura 7b)

#Extracción de los lncRNAs de la tabla general de cuentas y asignación por grupo 
normal <- general[,grep("MNR", colnames(general), value = TRUE)]
lnc_normal <- normal[rownames(lncRNA),]
lnc_normal<-  lnc_normal[rowSums(lnc_normal) >=10,]


cancer_plac <- general[,grep("Y|MNR", colnames(general), invert = TRUE, value = TRUE)]
cancer_plac
lnc_cancer_plac <- cancer_plac[rownames(lncRNA),]
lnc_cancer_plac <- lnc_cancer_plac[rowSums(lnc_cancer_plac) >=10,]

cancer_yodo <- general[,grep("Y", colnames(general), value = TRUE)]
cancer_yodo
lnc_cancer_yodo <- cancer_yodo[rownames(lncRNA),]
lnc_cancer_yodo <- lnc_cancer_yodo[rowSums(lnc_cancer_yodo) >=10,]

cancer_LA <- general[,c(1,2)]
lnc_LA_plac <- cancer_LA [rownames(lncRNA),]
lnc_LA_plac <- lnc_LA_plac[rowSums(lnc_LA_plac) >=10,]

cancer_LA_yodo <- general[,c(3,4,7)]
lnc_LA_yodo <- cancer_LA_yodo[rownames(lncRNA),]
lnc_LA_yodo <- lnc_LA_yodo[rowSums(lnc_LA_yodo) >=10,]

cancer_LB_yodo <- general[,c(8,9)]
lnc_LB_yodo <- cancer_LB_yodo[rownames(lncRNA),]
lnc_LB_yodo <- lnc_LB_yodo[rowSums(lnc_LB_yodo) >=10,]

cancer_TNC <- general[,c(5,6,12,13,14)]
lnc_TNC_plac <- cancer_TNC [rownames(lncRNA),]
lnc_TNC_plac <- lnc_TNC_plac[rowSums(lnc_TNC_plac) >=10,]

cancer_TNC_yodo <- general[,c(15:17)]
lnc_TNC_yodo <- cancer_TNC_yodo[rownames(lncRNA),]
lnc_TNC_yodo <- lnc_TNC_yodo[rowSums(lnc_TNC_yodo) >=10,]

#Matriz con el numero de lncRNA de todos los grupos 
lnc_muestras = matrix(c(length(rownames(lnc_normal)),
                        length(rownames(lnc_cancer_plac)),
                        length(rownames(lnc_LA_plac)), 
                        length(rownames(lnc_TNC_plac)),
                        length(rownames(lnc_cancer_yodo)),
                        length(rownames(lnc_LA_yodo)), 
                        length(rownames(lnc_LB_yodo)),
                        length(rownames(lnc_TNC_yodo))
)

, 
ncol = 1, byrow = TRUE)

rownames(lnc_muestras) <- c("Normal","Cáncer de mama","Luminal A", "Triple Negativo", "Cancer de mama Yodo", "Luminal A Yodo", "Luminal B Yodo","Triple Negativo Yodo")

#Asignación de de los colores a cada grupo 
colores_grupos <- c(
  "Normal" = "black",
  "Cáncer de mama" = "#E46780",
  "Luminal A" = "#7ED957",
  "Triple Negativo" = "#3F89D0",
  "Cancer de mama Yodo" = "#00E6E6",
  "Luminal A Yodo" = "#D633FF",
  "Luminal B Yodo"= "red",
  "Triple Negativo Yodo" = "#FFD700" )

colores_grupos <- c(
  expression("Normal"),
  expression("Cáncer de mama"),
  expression("Luminal A"),
  expression("Triple Negativo"),
  expression("Cáncer de mama"~I[2]),
  expression("Luminal A"~I[2]),
  expression("Luminal B"~I[2]),
  expression("Triple Negativo"~I[2])
)



# Extracción de  colores en el mismo orden que las filas que aparece en lnc_muestras
colores_barras <- colores_grupos[rownames(lnc_muestras)]

# Anotación con barras coloreadas
row_ha3 <- rowAnnotation(
  "N° de lncRNA expresados\n por grupo" = anno_barplot(
    lnc_muestras, 
    bar_width = 0.95,
    border = FALSE,
    gp = gpar(fill = colores_barras),  # <- aquí se aplican los colores
    axis_param = list(
      labels_rot = 80,
      gp = gpar(fontsize = 11)
    ),
    width = unit(4, "cm")
  )
)

rownames(lnc_muestras) <- gsub("Yodo", "I\u2082", rownames(lnc_muestras))



# Heatmap vacío con etiquetas de fila normales (puedes también colorearlas si quieres)
lnc_barras <- Heatmap(
  matrix(ncol = 0, nrow = nrow(lnc_muestras)),
  row_labels = rownames(lnc_muestras),
  row_names_side = "left",
  right_annotation = row_ha3,
  column_title = "LncRNAs",
  show_row_dend = FALSE,
  show_column_names = FALSE
)

draw(lnc_barras)


#Extracción de los PCG por grupo de la tabla general de conteos 
normal <- general[, grep("MNR", colnames(general), value = TRUE)]
pcg_normal <- normal[rownames(PCG), ]
pcg_normal <- pcg_normal[rowSums(pcg_normal) >= 10, ]


cancer_plac <- general[, grep("Y|MNR", colnames(general), invert = TRUE, value = TRUE)]
pcg_cancer_plac <- cancer_plac[rownames(PCG), ]
pcg_cancer_plac <- pcg_cancer_plac[rowSums(pcg_cancer_plac) >= 10, ]

cancer_yodo <- general[, grep("Y", colnames(general), value = TRUE)]
pcg_cancer_yodo <- cancer_yodo[rownames(PCG), ]
pcg_cancer_yodo <- pcg_cancer_yodo[rowSums(pcg_cancer_yodo) >= 10, ]

cancer_LA <- general[, c(1, 2)]
pcg_LA_plac <- cancer_LA[rownames(PCG), ]
pcg_LA_plac <- pcg_LA_plac[rowSums(pcg_LA_plac) >= 10, ]

cancer_LA_yodo <- general[, c(3, 4, 7)]
pcg_LA_yodo <- cancer_LA_yodo[rownames(PCG), ]
pcg_LA_yodo <- pcg_LA_yodo[rowSums(pcg_LA_yodo) >= 10, ]

cancer_LB_yodo <- general[, c(8, 9)]
pcg_LB_yodo <- cancer_LB_yodo[rownames(PCG), ]
pcg_LB_yodo <- pcg_LB_yodo[rowSums(pcg_LB_yodo) >= 10, ]


cancer_TNC <- general[, c(5, 6, 12, 13, 14)]
pcg_TNC_plac <- cancer_TNC[rownames(PCG), ]
pcg_TNC_plac <- pcg_TNC_plac[rowSums(pcg_TNC_plac) >= 10, ]

cancer_TNC_yodo <- general[, c(15:17)]
pcg_TNC_yodo <- cancer_TNC_yodo[rownames(PCG), ]
pcg_TNC_yodo <- pcg_TNC_yodo[rowSums(pcg_TNC_yodo) >= 10, ]


#Matriz con el número total de PCG de los grupos en el mismo orden que la lnc_muestras 
pcg_muestras = matrix(c(length(rownames(pcg_normal)),
                        length(rownames(pcg_cancer_plac)),
                        length(rownames(pcg_LA_plac)), 
                        length(rownames(pcg_TNC_plac)),
                        length(rownames(pcg_cancer_yodo)),
                        length(rownames(pcg_LA_yodo)),
                        length(rownames(pcg_LB_yodo)),
                        length(rownames(pcg_TNC_yodo))
),
ncol = 1, byrow = TRUE)
rownames(pcg_muestras) <- c("Normal","Cáncer de mama","Luminal A", "Triple Negativo", "Cancer de mama Yodo", "Luminal A Yodo", "Luminal B Yodo","Triple Negativo Yodo")
rownames(pcg_muestras) <- gsub("Yodo", "I\u2082", rownames(pcg_muestras))

row_ha2 <- rowAnnotation(
  "N° de PCG expresados\n por grupo" = anno_barplot(
    pcg_muestras, 
    bar_width = 0.95,                # Barras ocupan el 100% del espacio
    border= FALSE,
    gp = gpar(fill = colores_barras),# Sin bordes
    axis_param = list(
      labels_rot = 80,
      gp = gpar(fontsize = 11)),
    width = unit(4, "cm")           # Ajustar según necesidad
  ))



pcg_barras <- Heatmap(matrix(ncol = 0, nrow = nrow(pcg_muestras)),
             row_labels = c("","","","","","","",""),
             row_names_side = c("left"),
             right_annotation = row_ha2,
             column_title = "PCG",
             show_row_dend = FALSE, show_column_names = FALSE)




png("./exploratorios.png", width = 8, height = 6, units = "in", res = 300)

grid.newpage()

# Definir layout: 1 fila, 2 columnas
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))

#Dibujar el heatmap de los lncRNAs
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(lnc_barras, newpage = FALSE)
popViewport()

#Dibujar el heatmap de los PCG
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(pcg_barras, newpage = FALSE)
popViewport()
dev.off()


#Resumen visual DEGs grupos placebos y I2 (Figura 9b)

#Matrices del número de genes y  lncRNAs sobreexpresados y reprimidos por grupo

m = matrix(c(length(rownames(genes_up_cancer_norm)), length(rownames(genes_down_cancer_norm)),
             length(rownames(genes_up_Norm_LA)), length(rownames(genes_down_Norm_LA)),
             length(rownames(genes_up_Norm_TNC)), length(rownames(genes_down_Norm_TNC)),
             length(rownames(genes_up_yodo_norm)), length(rownames(genes_down_yodo_norm)),
             length(rownames(genes_up_Norm_LA_I)), length(rownames(genes_down_Norm_LA_I)),
             length(rownames(genes_up_Norm_TNC_Y)), length(rownames(genes_down_Norm_TNC_Y))
), 
ncol = 2, byrow = TRUE)


lncRNA_count = matrix(c(length(rownames(lncRNA_up_cancer_norm)), length(rownames(lncRNA_down_cancer_norm)),
                        length(rownames(lncRNA_up_Norm_LA)), length(rownames(lncRNA_down_Norm_LA)),
                        length(rownames(lncRNA_up_Norm_TNC)), length(rownames(lncRNA_down_Norm_TNC)),
                        length(rownames(lncRNA_up_cancer_I)), length(rownames(lncRNA_down_cancer_I)),
                        length(rownames(lncRNA_up_Norm_LA_I)), length(rownames(lncRNA_down_Norm_LA_I)),
                        length(rownames(lncRNA_up_Norm_TNC_Y)), length(rownames(lncRNA_down_Norm_TNC_Y))
), 
ncol = 2, byrow = TRUE)

# Columna de genes totales 
m_total= cbind(m, rowSums(m))

# Asignacion de nombres a las columnas
colnames(m_total) <- c("Up", "Down", "Total")
colnames(m) <- c("Up", "Down")




# Función para contar cuántos lncRNA o genes hay en cada lista
count_intersections <- function(up, down, reference) {
  up_count <- sum(up %in% reference)
  down_count <- sum(down %in% reference)
  return(c(up_count, down_count))
}


# DEG_down
DEG_down <- matrix(c(length(rownames(DEG_cancer[DEG_cancer$log2FoldChange <= -1 & 
                                                  complete.cases(DEG_cancer$padj) & 
                                                  DEG_cancer$padj <= 0.05,])),
                     length(rownames(DEG_LA[DEG_LA$log2FoldChange <= -1 & 
                                              complete.cases(DEG_LA$padj) & 
                                              DEG_LA$padj <= 0.05,])),
                     length(rownames(DEG_TNC[DEG_TNC$log2FoldChange <= -1 & 
                                               complete.cases(DEG_TNC$padj) & 
                                               DEG_TNC$padj <= 0.05,])),
                     length(rownames(DEG_yodo[DEG_yodo$log2FoldChange <= -1 & 
                                                complete.cases(DEG_yodo$padj) & 
                                                DEG_yodo$padj <= 0.05,])),
                     length(rownames(DEG_LAY[DEG_LAY$log2FoldChange <= -1 & 
                                               complete.cases(DEG_LAY$padj) & 
                                               DEG_LAY$padj <= 0.05,])),
                     length(rownames(DEG_TNCY[DEG_TNCY$log2FoldChange <= -1 & 
                                                complete.cases(DEG_TNCY$padj) & 
                                                DEG_TNCY$padj <= 0.05,]))),
                   ncol = 1, byrow = TRUE)

DEG_down

# DEG_up matrix
DEG_up <- matrix(c(length(rownames(DEG_cancer[(DEG_cancer$log2FoldChange >= 1 & 
                                                 complete.cases(DEG_cancer$padj) & 
                                                 DEG_cancer$padj <= 0.05),])),
                   length(rownames(DEG_LA[(DEG_LA$log2FoldChange >= 1 & 
                                             complete.cases(DEG_LA$padj) & 
                                             DEG_LA$padj <= 0.05),])),
                   length(rownames(DEG_TNC[(DEG_TNC$log2FoldChange >= 1 & 
                                              complete.cases(DEG_TNC$padj) & 
                                              DEG_TNC$padj <= 0.05),])),
                   length(rownames(DEG_yodo[(DEG_yodo$log2FoldChange >= 1 & 
                                               complete.cases(DEG_yodo$padj) & 
                                               DEG_yodo$padj <= 0.05),])),
                   length(rownames(DEG_LAY[(DEG_LAY$log2FoldChange >= 1 & 
                                              complete.cases(DEG_LAY$padj) & 
                                              DEG_LAY$padj <= 0.05),])),
                   length(rownames(DEG_TNCY[(DEG_TNCY$log2FoldChange >= 1 & 
                                               complete.cases(DEG_TNCY$padj) & 
                                               DEG_TNCY$padj <= 0.05),]))),
                 ncol = 1, byrow = TRUE)
DEG_up


#Cantidad de genes codificantes a proteinas 
CPG_col <- rbind(
  count_intersections(rownames(genes_up_cancer_norm), rownames(genes_down_cancer_norm), rownames(CPG)),
  count_intersections(rownames(genes_up_Norm_LA), rownames(genes_down_Norm_LA),rownames(CPG)),
  count_intersections(rownames(genes_up_Norm_TNC), rownames(genes_down_Norm_TNC),rownames(CPG)),
  count_intersections(rownames(genes_up_yodo_norm), rownames(genes_down_yodo_norm), rownames(CPG)),
  count_intersections(rownames(genes_up_Norm_LA_I), rownames(genes_down_Norm_LA_I),rownames(CPG)),
  count_intersections(rownames(genes_up_Norm_TNC_Y), rownames(genes_down_Norm_TNC_Y),rownames(CPG))
)


# Matriz general con toda la información de los biotipos 
DEG <- cbind(lncRNA_count, CPG_col,DEG_up,DEG_down)
DEG


# Asignar nombres de columna si lo deseas
colnames(DEG) <- c("lncRNA_Up", "lncRNA_Down", "CPG_Up", "CPG_Down","DEG_up","DEG_down")

DEG_down <- DEG[,c(2,4,6)]
DEG_up <- DEG[,c(1,3,5)]

DEG_up <- as.data.frame(DEG_up) %>% 
  mutate(
    Otros_Up = DEG_up - lncRNA_Up - CPG_Up
  )
DEG_up 
DEG_down <- as.data.frame(DEG_down) %>% 
  mutate(
    Otros_Up = DEG_down - lncRNA_Down - CPG_Down
  )

#Función de calculo de procentaje relativo de cada biotipo del total de DEGs
calcular_porcentajes_relativos_idx <- function(df, col_idxs, total_idx) {
  for (i in col_idxs) {
    nueva_col <- paste0(names(df)[i], "_pct")
    df[[nueva_col]] <- round(df[[i]] / df[[total_idx]] * 100, 2)
  }
  return(df)
}

DEG_up_porcentaje <- calcular_porcentajes_relativos_idx(DEG_up, col_idxs = c(1, 2, 4), total_idx = 3)
DEG_down_porcentaje <- calcular_porcentajes_relativos_idx(DEG_down, col_idxs = c(1, 2, 4), total_idx = 3)

DEG_up <- DEG_up[,c(1,2,4)]
colnames(DEG_up) <- c("lncRNA","CPG","Otros")
DEG_up
DEG_down <- DEG_down[,c(1,2,4)]
colnames(DEG_down) <- c("lncRNA","CPG","Otros")
DEG_down

DEG_up_p <- DEG_up_porcentaje %>% 
  mutate(
    number = paste(DEG_up,lncRNA_Up_pct, sep = " ("),
    number = paste(number,CPG_Up_pct,Otros_Up_pct, sep = "/"),
    number = paste(number,"",sep=")")
  )

DEG_down_p <- DEG_down_porcentaje %>% 
  mutate(
    number = paste(DEG_down,lncRNA_Down_pct, sep = " ("),
    number = paste(number,CPG_Down_pct,Otros_Up_pct, sep = "/"),
    number = paste(number,"",sep=")")
  )

DEG <- as.matrix(DEG)

# Conversion de la matriz DEG a numérica
col_biotypes = c("lncRNA" = "#1B9E77", "PCG" = "#7570B3", "Otros" = "#D95F02")

# Matriz de proporciones
biotype_mat <- as.matrix(DEG[, c("lncRNA", "PCG", "Otros")])


row_Biotype_ha1 <- rowAnnotation("N° DEG (UP%)" = row_anno_barplot(DEG_up_porcentaje[,5:7], axis = TRUE, 
                                                                   gp = gpar(fill = col_biotypes)), # color
                                 width = unit(1, "cm"))
row_Biotype_ha2 <- rowAnnotation("N° DEG (Down%)" = row_anno_barplot(DEG_down_porcentaje[,5:7], axis = TRUE, 
                                                                     gp = gpar(fill = col_biotypes)), # color
                                 width = unit(1, "cm"))

# Crear la anotación de los biotipos
biotype_anno <-  rowAnnotation(
  "Biotype composition" = row_anno_barplot(
    biotype_mat,               # Usamos biotype_mat para las proporciones
    gp = gpar(fill = col_biotypes), # Asignamos los colores a cada biotipo
    bar_width = 1,              # Ancho de las barras
    axis = TRUE                 # Mostrar el eje
  ),
  width = unit(4, "cm")         # Ancho de la barra lateral
)       # Ancho de la barra lateral

type_anno <- rowAnnotation("Tratamiento" = metadata$Group, col = list("Tratamiento" = col_type))
type_anno2 <- rowAnnotation("Grupo" = metadata$Typo, col = list("Grupo" = col_Typo))






group_ht <-  Heatmap(
  matrix(ncol = 0, nrow = nrow(m)),                         # Usamos una matriz dummy para mostrar las anotaciones
  name = "log2FC",                   # Nombre para la leyenda (aunque no haya datos reales)
  show_heatmap_legend = FALSE,       # No mostrar leyenda de valores de Heatmap
  row_split = metadata$Group, 
  column_title = " ",# Dividir por grupos (asegúrate de que metadata$Group existe)
  left_annotation = row_Biotype_ha1,
  col = col_biotypes                # Definir los colores para la leyenda
  # Añadir la anotación lateral con la barra apilada
)

lgd_list = list(Legend(labels = c("lncRNA", "PCG", "Otros"), title = "Biotipo", 
                       legend_gp = gpar(fill = col_biotypes), nrow = 1) # , title_position = "leftcenter")
)

m = t(apply(m, 1, function(x) x/sum(x)))
m <- as.data.frame(m)
colnames(m) <- c("UP","DOWN")
m <- cbind(m_total,m)

m <- m %>%
  mutate(
    Total = round(Total),
    Up = round(UP * 100, 2),      # Convertir a porcentaje
    Down = round(DOWN * 100, 2),
    number = paste(Total, Up, sep = " ("),
    number = paste(number,Down, sep = "/"),  # Formato requerido
    number = paste0(number, ")")  # Agregar el símbolo %
  )

#Anotaciones
col_type <- c("Placebo" = "#66A61E", "Yodo" = "yellow")

col_expr <- c("UP" = "firebrick3", "DOWN" = "dodgerblue3")

col_Typo <- c("Cáncer de mama" = "#E46780","Luminal A"="#7ED957","Triple Negativo"="#3F89D0")



metadata <- data.frame(
  Group = c("Placebo", "Placebo", "Placebo", "Yodo","Yodo","Yodo"),
  Typo = c("Cáncer de mama", "Luminal A", "Triple Negativo","Cáncer de mama", "Luminal A", "Triple Negativo")
)

type_anno <- rowAnnotation("Tratamiento" = metadata$Group, col = list("Tratamiento" = col_type))
type_anno2 <- rowAnnotation("Grupo" = metadata$Typo, col = list("Grupo" = col_Typo))


biotype_anno <- rowAnnotation(
  "Biotype" = anno_simple(
    col = list(Biotype = col_biotypes),   # Asigna los colores a la anotación
    which = "row"  # Asegura que se asigne a las filas
  )
)


# Definir anotaciones
group_labels <- metadata$Group

# Barra para DEG UP/DOWN
row_ha2 <- rowAnnotation("N°DEG(UP/DOWN%)" = row_anno_barplot(m[, c("Up", "Down")],
                                                              axis = TRUE,gp = gpar(fill = col_expr)),width = unit(4, "cm"))



# Anotación combinada a la izquierda
left_anno <- rowAnnotation(
  "Tratamiento" = metadata$Group,
  "Grupo" = metadata$Typo,
  col = list(
    "Tratamiento" = col_type,
    "Grupo" = col_Typo
  ),
  annotation_name_side = "top"
)
resumen_barras <- Heatmap(matrix(ncol = 0, nrow = nrow(m)), 
             row_labels = m[,6],
             right_annotation = row_ha2,
             left_annotation = left_anno,
             row_split = metadata$Group,
             show_row_dend = FALSE, 
             show_column_names = FALSE)


heatmap_with_barplot <- resumen_barras  +  
  rowAnnotation(rn = anno_text(m$number, gp = gpar(fontsize = 14), 
                               location = unit(0, "npc"), just = "left")) +  row_Biotype_ha1 +
  rowAnnotation(DOWN = anno_text(DEG_up_p$number, gp = gpar(fontsize = 14),
                                 location = unit(0, "npc"), just = "left"))+
  row_Biotype_ha2 +
  # DOWN genes (number)
  rowAnnotation(DOWN = anno_text(DEG_down_p$number, gp = gpar(fontsize = 14),
                                 location = unit(0, "npc"), just = "left"))


png("./DEG_summary.png", width = 10, height = 8, units = "in", res = 300)

draw(heatmap_with_barplot, 
     heatmap_legend_list = lgd_list,
     heatmap_legend_side = 'top',
     annotation_legend_side = 'top')

dev.off()



# Diagramas de Venn (Figuras 10a y c)

#Funcion para crear los diagramas de Venn 
display_venn <- function(x, title = "Diagrama de Venn", tamaño = 2, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL,
                              main = title, 
                              main.cex = tamaño,  # Tamaño del título
                              main.fontface = "bold", # Negrita
                              main.pos = c(0.5,1.2),
                              print.mode = "raw",
                              ...)
  grid.draw(venn_object)
}





# Lista con el nombre de los lncRNAs con expresión diferencial en los grupos placebos
todos_lnc =list(
  Cancer=rownames(lncRNA_norm_dif),
  LA=rownames(lncRNA_Norm_LA_dif),
  TNC=rownames(lncRNA_Norm_TNC_dif)
)



png("./DEG_cancer.png", width = 9, height = 8, units = "in", res = 300)
display_venn(
  todos_lnc,
  title = "Largos no codificantes sobreexpresados y reprimidos \n en el cáncer de mama y por subtipo molecular",
  category.names = c("Cáncer de mama","Luminal A ", "Triple Negativo"),
  lwd = 2,
  lty = 'blank',
  fill = c( "#E46780","#7ED957", "#3F89D0"),
  cex = 2,
  fontface = "italic",
  cat.cex = 2,
  cat.pos = c(-30, 30, 0), 
  cat.dist = c(0.05, 0.05, 0.03)
)
dev.off()

#Lista de lncRNAs expresados diferencialmente en los grupos I2
todos_yodo =list(
  Cancer=rownames(lncRNA_yodo_dif),
  LA=rownames(lncRNA_Norm_LA_I_dif),
  TNC=rownames(lncRNA_Norm_TNC_Y_dif)
)

category.names = c(
  "Cáncer de mama\n I₂",
  "Luminal A\n I₂",
  "Triple Negativo\n I₂"
)


png("./DEG_cancer_YODO.png", width = 9, height =8 , units = "in", res = 300)
display_venn(
  todos_yodo,
  title = "Largos no codificantes sobreexpresados y reprimidos \n en el cáncer de mama y por subtipo molecular\n inducidos por el I₂",
  category.names = c(
    "Cáncer de mama\n I₂",
    "Luminal A\n I₂",
    "Triple Negativo\n I₂"
  ),
  lwd = 2,
  lty = 'blank',
  fill = c("#00E6E6", "#D633FF", "#FFD700"),
  cex = 2,
  fontface = "italic",
  cat.cex = 2
)

dev.off()

#Diagramas UpserR (Figuras 10b y d)

#Preparación de los lncRNAs sobreexpresados y reprimidos por grupo 

###Cancer de mama Placebo###
a <- rownames(lncRNA_up_cancer_norm)
b<- rownames(lncRNA_down_cancer_norm)

###Luminal A Placebo###
c <- rownames(lncRNA_up_Norm_LA) 
d<- rownames(lncRNA_down_Norm_LA)

###TN Placebo###
e <- rownames(lncRNA_up_Norm_TNC)
f<- rownames(lncRNA_down_Norm_TNC)

###Cancer de mama I2###
g <- rownames(lncRNA_up_yodo_norm)
h<- rownames(lncRNA_down_yodo_norm)

###Luminal A I2###
i <- rownames(lncRNA_up_Norm_LA_I) 
j<- rownames(lncRNA_down_Norm_LA_I) 

###TN I2###
k <- rownames(lncRNA_up_Norm_TNC_Y) 
l<- rownames(lncRNA_down_Norm_TNC_Y) 


todos_pacebos= list(
  Sobreexpresado_cancer = a,
  Reprimidos_cancer = b,
  Sobreexpresado_LA = c,
  Reprimidos_LA = d,
  Sobreexpresado_TN = e,
  Reprimidos_TN = f
)

todos_yodo= list(
  Sobreexpresado_yodo = g,
  Reprimidos_yodo = h,
  Sobreexpresado_LA_yodo = i,
  Reprimidos_LA_yodo = j,
  Sobreexpresado_TN_yodo = k,
  Reprimidos_TN_yodo = l
)

cancer_LA =list(
  Sobreexpresado_cancer = a,
  Reprimidos_cancer = b,
  Sobreexpresado_LA = c,
  Reprimidos_LA = d
)

cancer_TN =list(
  Sobreexpresado_cancer = a,
  Reprimidos_cancer = b,
  Sobreexpresado_TN = e,
  Reprimidos_TN = f
)


LA_TN =list(
  Sobreexpresado_LA = c,
  Reprimido_LA = d,
  Sobreexpresado_TN = e,
  Reprimido_TN = f
)

LA_TNC_Yodo =list(
  Sobreexpresado_LA_Y = i,
  Reprimidos_LA_Y = j,
  Sobreexpresado_TN_Y = k,
  Reprimidos_TN_Y = l
)

library(UpSetR)
todos_binary <- fromList(todos_placebos)
todos_yodo_binary <- fromList(todos_yodo)
cancer_LA_binary <- fromList(cancer_LA)
cancer_TN_binary <- fromList(cancer_TN)
LA_TN_binary <- fromList(LA_TN)
LA_TN_yodo_binary <- fromList(LA_TN_Yodo)
plac_yodo_binary <- fromList(plac_yodo)
detach("package:UpSetR", unload = TRUE)
library(ComplexUpset)



names(LA_TN_Yodo) <- gsub("_Y$", "~I[2]", names(LA_TN_Yodo))
plac_yodo
names(plac_yodo) <- gsub("_yodo$", "-~I[2]", names(plac_yodo))

#Diagrama UpsetR con los lncRNAs sobreexpresados y reprimdos en común entre los subtipos LA y TN 

png("./UpsetR_LA_TNC.png", width = 8, height = 6, units = "in", res = 300)
upset(
  LA_TNC_binary,
  sort_intersections = FALSE,
  intersections = list(
    "Sobreexpresado_TN",
    "Reprimido_TN",
    "Sobreexpresado_LA",
    "Reprimido_LA",
    c("Sobreexpresado_LA", "Sobreexpresado_TN"),
    c("Reprimido_LA", "Reprimido_TN"),
    c("Sobreexpresado_TN", "Reprimido_LA"),
    c("Reprimido_TN", "Sobreexpresado_LA")
  ),
  height_ratio = 1,
  intersect = colnames(LA_TNC_binary),
  base_annotations = list(
    'Intersection size' = (
      intersection_size(
        text = list(size = 8)
      ) + ylab('LncRNAs compartidos')
    )
  ),
  themes = upset_modify_themes(
    list(
      'Intersection size' = theme(
        axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 12, face = 'italic')
      ),
      'intersections_matrix' = theme(
        text = element_text(size = 15)
      ),
      'overall_sizes' = theme(
        axis.text.x = element_text(angle = 90, size = 10)
      )
    )
  ),
  set_sizes = FALSE,
  queries = list(
    upset_query(
      intersect = c("Sobreexpresado_TN", "Sobreexpresado_LA"),
      color = 'red',
      fill = 'red'
    ),
    upset_query(
      intersect = c("Reprimido_TN", "Reprimido_LA"),
      color = 'blue',
      fill = 'blue'
    )
  ),
  wrap = TRUE
) +
  ggtitle('Por subtipo molecular') +
  theme(
    plot.title = element_text(size = 20, vjust = 1), # Tamaño y centrado del título
    axis.text.x = element_text(size = 70),      # Tamaño de los números en el eje X
    axis.text.y = element_text(size = 30),      # Tamaño de los números en el eje Y
    axis.title.y = element_text(size = 100)  
  )

dev.off()


#Diagrama UpsetR con los lncRNAs sobreexpresados y reprimdos en común entre los subtipos LA y TN trtados con I2

png("./UpsetR_LA_TN_Y.png", width = 8, height = 6, units = "in", res = 300)
ComplexUpset::upset(
  LA_TNC_yodo_binary,
  intersect = c(
    "Sobreexpresado_TN_Y",
    "Reprimidos_TN_Y",
    "Sobreexpresado_LA_Y",
    "Reprimidos_LA_Y"
  ),
  sort_intersections = FALSE,
  intersections = list(
    "Sobreexpresado_TN_Y",
    "Reprimidos_TN_Y",
    "Sobreexpresado_LA_Y",
    "Reprimidos_LA_Y",
    c("Sobreexpresado_LA_Y", "Sobreexpresado_TN_Y"),
    c("Reprimidos_LA_Y", "Reprimidos_TN_Y")
  ),
  height_ratio = 1,
  base_annotations = list(
    'Intersection size' = (
      intersection_size(text = list(size = 8)) +
        ylab('LncRNAs compartidos')
    )
  ),
  themes = upset_modify_themes(
    list(
      'Intersection size' = theme(
        axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 12, face = 'italic')
      ),
      'intersections_matrix' = theme(text = element_text(size = 15)),
      'overall_sizes' = theme(axis.text.x = element_text(angle = 90))
    )
  ),
  queries = list(
    upset_query(
      intersect = c("Sobreexpresado_TN_Y", "Sobreexpresado_LA_Y"),
      color = 'red',
      fill = 'red'
    ),
    upset_query(
      intersect = c("Reprimidos_TN_Y", "Reprimidos_LA_Y"),
      color = 'blue',
      fill = 'blue'
    )
  ),
  wrap = TRUE,
  set_sizes = FALSE
) +
  ggtitle('Por subtipo molecular') +
  theme(
    plot.title = element_text(size = 20, vjust = 1),
    axis.text.x = element_text(size = 70),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 100)
  )

dev.off()


#Construcción de los data frame para heatmaps (Figura12a)

res_norm_shrunken <- as.data.frame(res_norm_shrunken)
res_yodo_shrunken <- as.data.frame(res_yodo_shrunken)
res_Norm_LA_shrunken <- as.data.frame(res_Norm_LA_shrunken)
res_Norm_LA_I_shrunken <- as.data.frame(res_Norm_LA_I_shrunken)
res_Norm_TNC_shrunken <- as.data.frame(res_Norm_TN_shrunken)
res_Norm_TNC_Y_shrunken <- as.data.frame(res_Norm_TN_Y_shrunken)

#Dataframe grupo cáncer de mama 
LFC_general <- merge(res_norm_shrunken,res_yodo_shrunken, by=0)
LFC_general <- column_to_rownames(LFC_general,var="Row.names")
LFC_general <- LFC_general[,c(2,7)]
colnames(LFC_general) <- c("Placebo","Yodo")

write.table(LFC_general , file = "./LFC_general.txt",quote = FALSE)

#Dataframe grupo subtipo LA
LFC_LA <- merge(res_Norm_LA_shrunken,res_Norm_LA_I_shrunken, by=0)
LFC_LA <- column_to_rownames(LFC_LA,var="Row.names")
LFC_LA <- LFC_LA[,c(2,7)]
colnames(LFC_LA) <- c("Placebo","Yodo")

write.table(LFC_LA , file = "./LFC_LA.txt",quote = FALSE)

#Dataframe grupo subtipo TN
LFC_TNC <- merge(res_Norm_TNC_shrunken,res_Norm_TNC_Y_shrunken, by=0)
LFC_TNC<- column_to_rownames(LFC_TNC,var="Row.names")
LFC_TNC <- LFC_TNC[,c(2,7)]
colnames(LFC_TNC) <- c("Placebo","Yodo")
write.table(LFC_TNC, file = "./LFC_TNC.txt",quote = FALSE)

#Función para crear el heatmap: Dataframe con el valor de LFC del grupo placebo y tratmiento, dataframe de los genes a 
#visualizar con los nombres como nombre de fila , tamaño del titulo, nombre del titulo, tamaño, y row names

crear_heatmap_lncRNA <- function(datos, lncRNA_up_cancer_norm, column_title = "Titulo", titulo, tamaño = 5, rows= TRUE) {
  # Unión del dataframe con el valor de LFC con los genes con expresión aberrante 
  lncRNA_up_ambos <- merge(datos, lncRNA_up_cancer_norm, by = 0, sort = FALSE)
  row.names(lncRNA_up_ambos) <- lncRNA_up_ambos$Row.names
  lncRNA_up_ambos <- lncRNA_up_ambos[, 2:3]
  colnames(lncRNA_up_ambos) <- c("Placebo","Yodo")
  
  col_fun2 <- colorRamp2(
    c(-9, -5, -1, -0.5, 0, 0.5, 1, 5, 9),c("darkblue", "blue", "lightblue", "lightcyan", "lightyellow", "yellow", "orange", "red", "darkred")
  )
  
  
  # Crear el heatmap
  heatmap <- Heatmap(
    as.matrix(dataframe),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = rows,
    row_names_gp = gpar(fontsize = tamaño),   # Tamaño de las etiquetas de las filas
    column_names_gp = gpar(fontsize = 12),   # Tamaño de las etiquetas de las columnas
    heatmap_legend_param = list(
      title = "Log2FC",                      # Título de la leyenda
      at = c(-9, -5, -1, -0.5, 0, 0.5, 1, 5, 9),       # Puntos específicos en la escala
      labels = c("-9", "-5", "-1", "-0.5", "0", "0.5", "1", "5", "9")
    ),
    name = "lncRNAs",                          # Nombre del heatmap
    col = col_fun2,                          # Escala de colores
    column_title = column_title,             # Título del heatmap
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  return(heatmap)
}

#Creación de los Heatmaps para lncRNAs desregulados en los placebos (Figura 12a y b )
png("./lnc_cancer_DEG.png", width = 8, height = 6, units = "in", res = 300)
crear_heatmap_lncRNA(LFC_general,lncRNA_norm_dif, 
                     column_title = "Cáncer de mama", 5,tamaño=10,FALSE)
dev.off()

#Diagramas de Venn de los lncRNAs modulados por el I2 (Figura 13a)

#Comparaciones entre los lncRNAs expresados diferencialmente 
#en los grupos placebos con los que no cambiaron el I2 y
#tienen significancia estadistica 

#lncRNAs que disminuyen con el I2 del grupo placebo 
lncRNA_up_cancer_down_yodo <- merge(lncRNA_up_cancer_norm,lncRNA_nochange_yodo_norm, by=0)
lncRNA_up_cancer_down_yodo <- column_to_rownames(lncRNA_up_cancer_down_yodo,var="Row.names")


#lncRNAs que aumentan con el I2 del grupo placebo 
lncRNA_down_cancer_up_yodo <- merge(lncRNA_down_cancer_norm, lncRNA_nochange_yodo_norm, by=0)
lncRNA_down_cancer_up_yodo  <- column_to_rownames(lncRNA_down_cancer_up_yodo ,var="Row.names")

#lncRNAs que disminuyen con el I2 del grupo  LA 
lncRNA_up_LA_down_yodo <- merge(lncRNA_up_Norm_LA,lncRNA_nochange_yodo_LA, by=0)
lncRNA_up_LA_down_yodo <- column_to_rownames(lncRNA_up_LA_down_yodo ,var="Row.names")

#lncRNAs que aumentan con el I2 del grupo  LA
lncRNA_down_LA_up_yodo <- merge(lncRNA_down_Norm_LA,lncRNA_nochange_yodo_LA, by=0)
lncRNA_down_LA_up_yodo  <- column_to_rownames(lncRNA_down_LA_up_yodo ,var="Row.names")

#lncRNAs que disminuyen con el I2 del grupo  TN
ncRNA_up_TNC_down_yodo <- merge(lncRNA_up_Norm_TNC,lncRNA_nochange_yodo_TNC, by=0)
lncRNA_up_TNC_down_yodo <- column_to_rownames(lncRNA_up_TNC_down_yodo ,var="Row.names")


#lncRNAs que aumentan con el I2 del grupo  TN
lncRNA_down_TNC_up_yodo <- merge(lncRNA_down_Norm_TNC,lncRNA_nochange_yodo_TNC, by=0)
lncRNA_down_TNC_up_yodo  <- column_to_rownames(lncRNA_down_TNC_up_yodo ,var="Row.names")

#Union de los modulados en un solo dataframe 
lnc_cancer_mod <- rbind(lncRNA_up_cancer_down_yodo, lncRNA_down_cancer_up_yodo)
lnc_cancer_mod

lnc_LA_mod <- rbind(lncRNA_up_LA_down_yodo, lncRNA_down_LA_up_yodo)
lnc_LA_mod

lnc_TNC_mod <- rbind(lncRNA_up_TNC_down_yodo, lncRNA_down_TNC_up_yodo)
lnc_TNC_mod

#Unión de los lncRNAs modulados en una sola lista 
lnc_modulados =list(
  cancer=rownames(lnc_cancer_mod),
  LA=rownames(lnc_LA_mod),
  TN=rownames(lnc_TNC_mod)
)

png("./lncRNAs_moduladosI2.png", width = 8, height = 6, units = "in", res = 300)
display_venn(
  lnc_modulados,
  title = expression(bold("Largos no codificantes modulados por el" ~ I[2])),
  category.names = c("Cáncer de mama" , "Luminal A " , "Triple Negativo"),
  lwd = 2,
  lty = 'blank',
  fill = c("#E46780", "#7ED957","#3F89D0"),
  cex = 2,
  fontface = "italic",
  cat.cex = 1.5
) 
dev.off()

#Diagramas de Venn de los lncRNAs sobreexpresados por el I2 (Figura 13b)

#Comparaciones entre los lncRNAs que no estaban expresados diferencialmente 
#en los grupos placebos con los que cambiaron en los grupos I2 


#LncRNas grupo cáncer de mama I2
lncRNA_nuevos_yodo <- setdiff(row.names(lncRNA_yodo_dif),row.names(lncRNA_norm_dif))

lncRNA_nuevos_yodo <- setdiff(lncRNA_nuevos_yodo,rownames(lncRNA_nochange_yodo_norm))
lncRNA_nuevos_yodo <- as.data.frame(lncRNA_nuevos_yodo)
colnames(lncRNA_nuevos_yodo) <-"Row.names"

lncRNA_nuevos_yodo <- column_to_rownames(lncRNA_nuevos_yodo, var = "Row.names")

lnc_nuevos_yodo <- merge(LFC_general,lncRNA_nuevos_yodo, by=0)

#Eliminación de aquellos que que no tengan una diferfencia de 1 en el LFC
#limpiando aquellos que si tenian una desregulación en el LFC, pero no la significancia estadistica del grupo placebo 
lnc_nuevos_yodo <- lnc_nuevos_yodo %>%  
  mutate(Diff = abs(Yodo - Placebo)) %>%  
  filter(Diff >= 0.5, between(Placebo, -0.49, 0.49))


lnc_nuevos_yodo <- column_to_rownames(lnc_nuevos_yodo, var = "Row.names")
lnc_nuevos_yodo <- lnc_nuevos_yodo[,-3]
colnames(lnc_nuevos_yodo) <- c("Placebo","I~[2]")


#LncRNas grupo LA I2
lncRNA_nuevos_yodo_LA <- setdiff(row.names(lncRNA_Norm_LA_I_dif),row.names(lncRNA_Norm_LA_dif))
lncRNA_nuevos_yodo_LA <- setdiff(lncRNA_nuevos_yodo_LA,rownames(lncRNA_nochange_yodo_LA))

lncRNA_nuevos_yodo_LA <- as.data.frame(lncRNA_nuevos_yodo_LA)
colnames(lncRNA_nuevos_yodo_LA) <-"Row.names"

lncRNA_nuevos_yodo_LA <- column_to_rownames(lncRNA_nuevos_yodo_LA, var = "Row.names") 

lncRNA_nuevos_yodo_LA <- merge(LFC_LA,lncRNA_nuevos_yodo_LA, by=0)

lncRNA_nuevos_yodo_LA <- lncRNA_nuevos_yodo_LA  %>%  
  mutate(Diff = abs(Yodo - Placebo)) %>%  
  filter(Diff >= 0.5, between(Placebo, -0.49, 0.49))

#LncRNas grupo TN I2
lncRNA_nuevos_yodo_TNC <- setdiff(row.names(lncRNA_Norm_TNC_Y_dif),row.names(lncRNA_Norm_TNC_dif))
lncRNA_nuevos_yodo_TNC <- setdiff(lncRNA_nuevos_yodo_TNC,rownames(lncRNA_nochange_yodo_TNC))

lncRNA_nuevos_yodo_TNC <- as.data.frame(lncRNA_nuevos_yodo_TNC)
colnames(lncRNA_nuevos_yodo_TNC) <-"Row.names"

lncRNA_nuevos_yodo_TNC <- column_to_rownames(lncRNA_nuevos_yodo_TNC, var = "Row.names")

lncRNA_nuevos_yodo_TNC <- merge(LFC_TNC,lncRNA_nuevos_yodo_TNC, by=0)

lncRNA_nuevos_yodo_TNC <- lncRNA_nuevos_yodo_TNC  %>%  
  mutate(Diff = abs(Yodo - Placebo)) %>%  
  filter(Diff >= 0.5, between(Placebo, -0.49, 0.49))


lnc_exp_yodo = list(
  cancer_nuevos = rownames(lnc_yodo_nuevos),
  TNC_nuevos =rownames(lncRNA_nuevos_yodo_TNC),
  LA_nuevos =rownames(lncRNA_nuevos_yodo_LA )
)


png("/lncRNA_yodo_new.png", width = 9, height = 8, units = "in", res = 300)
display_venn(
  lnc_exp_yodo,
  title =   "Largos no codificantes expresados diferencialmente \n debido al I",
  category.names = c(
    "Cáncer de mama I₂",
    "Triple Negativo I₂",
    "Luminal A I₂"
  ),
  lwd = 4,
  lty = 'blank',
  fill = c("#00E6E6","#FFD700","#D633FF" ),
  cex = 1.5,
  fontface = "italic",
  cat.cex = 1.3) 
dev.off()