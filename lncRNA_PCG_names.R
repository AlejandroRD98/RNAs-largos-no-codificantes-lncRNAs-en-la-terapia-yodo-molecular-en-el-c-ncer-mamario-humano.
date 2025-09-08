
# Lectura del archivo con los nombres de lncRNs extraidos de genode
lncRNA <- readLines("./lnc_names.txt")
lncRNA <- as.data.frame(lncRNA)
colnames(lncRNA) <- "lncRNA"
row.names(lncRNA) <- NULL
lncRNA <- column_to_rownames(lncRNA, var="lncRNA")

# Lectura del archivo con los nombre de los PCG extraidos de gencode
PCG <- readLines("./protein_cosing_names.txt")
PCG <- as.data.frame(PCG)
row.names(PCG) <- NULL
PCG <- column_to_rownames(PCG, var="PCG")






