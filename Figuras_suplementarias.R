library(ggVennDiagram)
library(ggplot2)

#Diagramas de Veen con los genes expresados diferencialmente en comun entre placebos 
#y I2 (Figura Suplementaria 2)

display_venn <- function(x, title = "Diagrama de Venn", tamaño = 2, ...){
  library(VennDiagram)
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


LA_LAY =list(
  LA=rownames(Norm_LA_dif),
  LA_Y=rownames(Norm_LA_I_dif)
)


png("C:/Users/Alex Ruiz/Desktop/Maestria/seqNueva/0502/graficas/DEG_LA_LAY.png", width = 9, height = 8, units = "in", res = 300)
display_venn(
  LA_LAY,
  title = "",
  category.names = c("Luminal A","Luminal A I₂"),
  lwd = 4,
  lty = 'blank',
  fill = c("#7ED957","#D633FF" ),
  cex = 1.5,
  fontface = "italic",
  cat.cex = 1.3) # Mueve las etiquetas para que no se superpongan
# Aumenta la distancia de las etiquetas respecto al diagrama
dev.off()


TNC_TNCY =list(
  TNC=rownames(Norm_TNC_dif),
  TNC_Y=rownames(Norm_TNC_Y_dif)
)

png("C:/Users/Alex Ruiz/Desktop/Maestria/seqNueva/0502/graficas/DEG_TNC_TNCY.png", width = 9, height = 8, units = "in", res = 300)
display_venn(
  TNC_TNCY,
  title = "",
  category.names = c("Triple negativo","Triple negativo\n I₂"),
  lwd = 4,
  lty = 'blank',
  fill = c("#3F89D0","#FFD700" ),
  cex = 1.5,
  fontface = "italic",
  cat.cex = 1.3) # Mueve las etiquetas para que no se superpongan
# Aumenta la distancia de las etiquetas respecto al diagrama
dev.off()

can_canY =list(
  can=rownames(norm_dif),
  can_Y=rownames(yodo_dif)
)


png("C:/Users/Alex Ruiz/Desktop/Maestria/seqNueva/0502/graficas/DEG_can_canY.png", width = 9, height = 8, units = "in", res = 300)
display_venn(
  can_canY,
  title = "Genes expresados diferencialmente entre \n grupos placebo y tratados con I₂",
  category.names = c("Cáncer de mama","Cáncer de mama \nI₂"),
  lwd = 4,
  lty = 'blank',
  fill = c("#E46780","#00E6E6" ),
  cex = 1.5,
  fontface = "italic",
  cat.cex = 1.3) # Mueve las etiquetas para que no se superpongan
# Aumenta la distancia de las etiquetas respecto al diagrama
dev.off()

#Función para crear las tablas con los lncRNAs
#expresados diferencialmente por subtipo molecular
#coloreando de rojo los sobreexpresados y de azul los reprimidos 
#(Figura suplementaria 3,6 y 7)

crear_tabla_coloreada <- function(df, num_columnas, titulo, rojos, azules) {
  # Calcular el número total de celdas necesarias
  total_celdas <- num_columnas * ceiling(length(df) / num_columnas)
  
  # Rellenar con espacios vacíos para completar la tabla
  nombres_filas <- c(df, rep("", total_celdas - length(df)))
  
  # Crear la matriz de nombres
  matriz_nombres <- matrix(nombres_filas, ncol = num_columnas, byrow = TRUE)
  colnames(matriz_nombres) <- NULL  # Evitar encabezados repetidos
  
  # Crear vector de colores
  colores <- rep("black", length(nombres_filas))  # Negro por defecto
  colores[matriz_nombres %in% azules] <- "blue"   # Azul para genes reducidos
  colores[matriz_nombres %in% rojos] <- "red"     # Rojo para genes sobreexpresados
  colores[matriz_nombres == ""] <- "black"        # Negro para celdas vacías
  
  # Definir el tema de la tabla
  tema_personalizado <- ttheme_default(
    core = list(
      fg_params = list(fontsize = 10, col = colores),  # Colores dinámicos
      bg_params = list(fill = c("white", "#56B4E9"), alpha = 0.5)  # Alternar colores de fondo
    ),
    colhead = list(
      fg_params = list(fontsize = 12, fontface = "bold", col = "white"),  # Encabezado en negritas y blanco
      bg_params = list(fill = "black", alpha = 0.8)  # Fondo del encabezado
    )
  )
  
  # Crear la tabla con el tema personalizado
  tabla <- tableGrob(matriz_nombres, theme = tema_personalizado)
  
  # Crear el encabezado general
  encabezado <- textGrob(titulo, gp = gpar(fontsize = 14, fontface = "bold", col = "black"))
  
  # Agregar el encabezado a la tabla
  tabla <- gtable_add_rows(tabla, heights = unit(1.5, "lines"), pos = 0)
  tabla <- gtable_add_grob(tabla, encabezado, t = 1, l = 1, r = num_columnas)
  
  # Dibujar la tabla en la pantalla
  grid.arrange(tabla)
}

#crear una lista de los lncRNA expresados diferencialmente en todos los grupos
#Para despues comprar y extraer los comúnes y especificos.
lnc_todos =list(
  LA=rownames(lncRNA_Norm_LA_dif),
  Cancer=rownames(lncRNA_norm_dif),
  TN=rownames(lncRNA_Norm_TNC_dif)
)

#Calcular los compartidos 
intersecciones <- calculate.overlap(todos)
intersecciones 
#Extraer los especificos del cancer de mama y subtipo especifico
#cambiar dependiendo el contenido de intersecciones_todos 
lncRNAs_comunes_cancer  <- intersecciones$a5 
lncRNAs_la <- intersecciones$a1 
lncRNAs_tnc <- intersecciones$a7


#Generar la tabla y colorear dependiendo el valoer de su LFC
crear_tabla_coloreada(lncRNAs_comunes_cancer, 
                      num_columnas = 4, 
                      titulo = "DElncRNAs cáncer de mama", 
                      rojos = rownames(lncRNA_norm_dif[lncRNA_norm_dif$log2FoldChange>0,]), 
                      azules = rownames(lncRNA_down_cancer))


crear_tabla_coloreada(lncRNAs_la, 
                     num_columnas = 10,
                     titulo = "DElncRNAs luminal A", 
                     rojos = rownames(lncRNA_Norm_LA_dif[lncRNA_Norm_LA_dif$log2FoldChange>0,]), 
                     azules = rownames(lncRNA_Norm_LA_dif[lncRNA_Norm_LA_dif$log2FoldChange<0,]))

png("./DElncRNAS_TNC.png", width = 20, height = 20, units = "in", res = 300)  # Inicia la imagen
crear_tabla_coloreada(lncRNAs_tnc, 
                      num_columnas = 12,
                      titulo = "DElncRNAs triple negativo", 
                      rojos = rownames(lncRNA_Norm_TNC_dif[lncRNA_Norm_TNC_dif$log2FoldChange>0,]), 
                      azules = rownames(lncRNA_Norm_TNC_dif[lncRNA_Norm_TNC_dif$log2FoldChange<0,]))
dev.off()


# extracción de las cuentas normalizadas 
norm_counts <- counts(ddsDE,normalized=TRUE)


#Extracción de las cuentas de los lncRNAs comunes en el cáncer de mama 
lncRNAs_comunes_cancer_counts <- norm_counts[lncRNAs_comunes_cancer,]

# Conversion del dataframe de expresión con las cuentas a formato largo para 
#graficar 

lncRNAs_comunes_cancer_long <- as.data.frame(lncRNAs_comunes_cancer_counts)  %>%
  rownames_to_column(var = "Gene") %>%   # Convertir rownames a columna explícita
  pivot_longer(cols = -Gene,             # Todas las columnas excepto 'Gene'
               names_to = "sample", 
               values_to = "Expresión")

# Unir con la información de grupos la información de la metadata   
lncRNAs_comunes_cancer_long  <- lncRNAs_comunes_cancer_long  %>%
  left_join(datos[, c("sample", "Grupo","Tratamiento")], by = "sample")

#Agregar información de la expresión (sobreexpresado = 1)
lncRNAs_comunes_cancer_long <- lncRNAs_comunes_cancer_long %>%
  mutate(Expresión = ifelse(Expresión == 0, 1, Expresión)) %>%
  filter(Grupo %in% c("NORM", "LA", "TN","LA-I₂","TN-I₂"))
         
# Convertir 'Grupo' a factor con orden deseado para graficar 
lncRNAs_comunes_cancer_long$Grupo <- factor(lncRNAs_comunes_cancer_long$Grupo, 
                              levels = c("NORM", "LA", "TN","LA-I₂","TN-I₂"))

#Agrefar información de LFC 
lncRNAs_comunes_cancer_long <- lncRNAs_comunes_cancer_long %>%
  left_join(lncRNA_norm_dif[, c("Gene", "log2FoldChange")], by = "Gene") %>%
  mutate(Direccion = ifelse(log2FoldChange > 0, "Up", "Down"))  # Agregar dirección

#Establecer colores de los grupos 
colores_grupo <- c(
  "NORM" = "black",
  "LA" = "#7ED957",
  "TN" = "#3F89D0",
  "LA-I₂" = "#D633FF",
  "TN-I₂" = "#FFD700"
  
)

#orden de los lncRNAs para que sea de forma ascendente los sobreexpresados 
orden_genes <- lncRNAs_comunes_cancer_long  %>%
  filter(Direccion == "Up", Grupo == "NORM") %>%
  group_by(Gene) %>%
  summarise(prom_exp = mean(log2(Expresión), na.rm = TRUE)) %>%
  arrange(prom_exp) %>%
  pull(Gene)

#Generar la grafica de los lncRNAs sobreexpresados en el cáncer de mama 
up_comun <- ggplot(lncRNAs_comunes_cancer_long [lncRNAs_comunes_cancer_long $Direccion == "Up",], aes(x = Gene,
                                                            y = log2(Expresión),
                                                            fill = Grupo)) +
  geom_boxplot(width = 0.9, outlier.shape = NA) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  theme_minimal(base_size = 10) +
  labs(title = "LncRNAs sobreexpresados en el cáncer de mama",
       x = "Genes",
       y = expression(log[2]("Expresión Normalizada"))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = colores_grupo)


#orden de los lncRNAs para que sea de forma descendente de los reprimidos 
orden_genes_down <- lncRNAs_comunes_cancer_long  %>%
  filter(Direccion == "Down", Grupo == "NORM") %>%
  group_by(Gene) %>%
  summarise(prom_exp = mean(log2(Expresión), na.rm = TRUE)) %>%
  arrange(desc(prom_exp)) %>%
  pull(Gene)



lncRNAs_comunes_cancer_long$Gene <- factor(lncRNAs_comunes_cancer_long $Gene, levels = orden_genes_down)

# Paso 3: graficar con genes ordenados
down_comun <- ggplot(df_long[df_long$Direccion == "Down",], aes(x = Gene,
                                                                y = log2(Expresión),
                                                                fill = Grupo)) +
  geom_boxplot(width = 0.9, outlier.shape = NA) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  theme_minimal(base_size = 10) +
  labs(title = "LncRNAs reprimidos en el cáncer de mama",
       x = "Genes",
       y = expression(log[2]("Expresión Normalizada"))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank()) +
  scale_fill_manual(values = colores_grupo)



#Conbinar ambas graficas en una sola 
png("./lncRNA_cancer_exp.png", width = 16, height = 18, units = "in", res = 300)
up_comun / down_comun
dev.off()


#Grafico UpsetR con los genes sobreexpresados y reprimidos de todos
#los grupos placebos y yodo (Figura suplementaria 5 )

plac_yodo = list(
  Sobreexpresado_cancer = a,
  Reprimido_cancer = b,
  Sobreexpresado_LA = c,
  Reprimido_LA = d,
  Sobreexpresado_TN = e,
  Reprimido_TN = f,
  Sobreexpresado_yodo = g,
  Reprimido_yodo = h,
  Sobreexpresado_LA_yodo = i,
  Reprimido_LA_yodo = j,
  Sobreexpresado_TN_yodo = k,
  Reprimido_TN_yodo = l)


library(UpSetR)
plac_yodo_binary <- fromList(plac_yodo)
detach("package:UpSetR", unload = TRUE)
library(ComplexUpset)

upset(
  plac_yodo_binary,
  colnames(plac_yodo_binary),width_ratio=1,
  set_sizes = FALSE,
  queries = list(
    upset_query(
      intersect = c("Sobreexpresado_TN-~I[2]", "Sobreexpresado_LA_-~I[2]","Sobreexpresado_-~I[2]","Sobreexpresado_cancer","Sobreexpresado_TN",
                    "Sobreexpresado_LA"),
      color = 'red',
      only_components=c('intersections_matrix')
    ),
    upset_query(
      intersect = c("Reprimido_TN-~I[2]", "Reprimido_LA-~I[2]", "Reprimido_-~I[2]",
                    "Reprimidos_cancer", "Reprimidos_TN", "Reprimidos_LA"
      ),
      color = 'blue',
      only_components=c('intersections_matrix')
      
    )
  ))+
  labs(
    title = "LncRNAs sobreexpresados y reprimidos",
    x = "Grupo",                   # Cambia el nombre del eje X
    y = "Tamaño de Intersección"   # Cambia el nombre del eje Y
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 1), # Tamaño, centrado y ajuste vertical del título
    axis.text = element_text(size = 10),      # Aumenta el tamaño de los números en los ejes
    axis.title = element_text(size = 16),    # Aumenta el tamaño de los títulos de los ejes     # Aumenta el tamaño del título (si lo hay)
    axis.text.y = element_text(size = 12)   # Tamaño de los números en el eje Y
  )

