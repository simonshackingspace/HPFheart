library(ggplot2)
library(umap)

cell_scores <- read.csv("/home/zz2565/Documents/thesis_code/thesis_output/scHPF_best_score/cell_score.txt", sep='\t', header=F)
cell_types <- read.csv("/home/zz2565/data/thesis_data/diff_seqs/celltypes.csv", header=F)$V2

cell_scores.umap <- umap::umap(cell_scores)
df <- data.frame(x = cell_scores.umap$layout[,1],
                 y = cell_scores.umap$layout[,2],
                 celltype = cell_types)
                 
ggplot(df, aes(x, y, colour = celltype)) + geom_point()


