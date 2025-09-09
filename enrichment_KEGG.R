# Cargar librerías necesarias
library(clusterProfiler)  
library(org.Hs.eg.db)     
library(ggplot2)          


gene_symbols <- rownames(genes)
# Convertir símbolos a ENTREZ IDs
genes_entrez <- bitr(gene_symbols, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Enriquecimiento KEGG
kegg_result <- enrichKEGG(gene = genes_entrez$ENTREZID,
                          organism = 'hsa',
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

#Grafica barplot de los rutas 
barplot(kegg_result, showCategory = 20) +
  ggtitle(bquote(bold("Rutas enriquecidas"))) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

