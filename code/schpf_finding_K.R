library(pheatmap)
library(grid)
library(dplyr)

args <- commandArgs(trailingOnly=T)
cell_scores_file <- args[1]
K <- args[2]
meta_data_file <- args[3]

cell_scores <- read.csv(cell_scores_file,
                         sep='\t', header=F)
meta_data = read.csv(meta_data_file, 
                      header=F,
                      stringsAsFactors=T,
                      sep=",",
                      row.names=1)
cell_scores$celltype = meta_data[,1]
cell_factor <- aggregate(.~celltype, cell_scores, mean)
rownames(cell_factor) <- cell_factor$celltype
cell_factor$celltype=NULL
colnames(cell_factor) <- c(1:length(colnames(cell_factor)))

draw_colnames_rite <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.25 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(7,"bigpts"), vjust = 0.5, hjust = 1, rot = 0, gp = gpar(...))
  return(res)}
 
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_rite",
                 ns=asNamespace("pheatmap"))

original_wd = getwd()
setwd("thesis_imgs/tmp_imgs")
pheatmap(as.matrix(cell_factor), 
           cellwidth=15, 
           cellheight=15,
           main=paste("K", K, sep="="),
           filename=paste(paste("K", K, sep="_"), "png", sep="."))
setwd(original_wd)
 
