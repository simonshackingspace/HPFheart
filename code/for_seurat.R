library(Matrix)

sc_data = read.csv("/home/zz2565/data/thesis_data/diff_seqs/sc_data.csv", 
                   stringsAsFactors=F,
                   header=TRUE, 
                   sep=",",
                   row.names=1)
proteins = scan("/home/zz2565/data/thesis_data/protein_encoding_genes.txt", character())
sc_data = sc_data[, colnames(sc_data) %in% proteins]

sparse.sc_data = Matrix(as.matrix(sc_data), sparse=T)
    
writeMM(obj = sparse.sc_data, 
        file="/home/zz2565/data/thesis_data/sc_heart_seurat_matrix/matrix.mtx")
write(x = rownames(sc_data),
      file = "/home/zz2565/data/thesis_data/sc_heart_seurat_matrix/features.tsv")
write(x = colnames(sc_data),
      file = "/home/zz2565/data/thesis_data/sc_heart_seurat_matrix/barcodes.tsv")
    


