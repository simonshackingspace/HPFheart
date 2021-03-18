library(pheatmap)
library(grid)
library(readxl)
library(dplyr)

genes <- read.csv("/home/zz2565/data/thesis_data/diff_seqs/scHPF_pretrain/genes.txt", sep='\t', header=F)
gene_ids <- genes[,2]
gene_scores <- read.csv("thesis_output/scHPF_best_score/gene_score.txt",
                        sep='\t', header=F)
colnames(gene_scores) <- c(1:dim(gene_scores)[2])
rownames(gene_scores) <- gene_ids

# marker genes
marker_genes <- c("VWF", "ACTN2", "TNNI3", "CASP3", "ACTA2", "CKAP4","OGN", "ALDH1A2", "ITLN1")
marker_genes_scores <- subset(gene_scores, rownames(gene_scores)%in%marker_genes)

# CHD genes
# mess_genes <- c("GDF1", "MYH6", "FLT4")
CHD_genes_file <- read_xlsx("/home/zz2565/data/thesis_data/41588_2017_BFng3970_MOESM3_ESM.xlsx", sheet=2)
CHD_genes <- dplyr::pull(CHD_genes_file, 1)[2:254]
# new_genes <- c("AKAP12", "ANK3", "CLUH", "CTNNB         1", "KDM5A", "KMT2C", "MINK1", "MYRF", "PRRC2B","RYR3", "U2SURP", "WHSC1")
# CHD_genes <- c(CHD_genes, new_genes)
CHD_genes_scores <- subset(gene_scores, rownames(gene_scores)%in%CHD_genes)

# de novo mutation
try <- read_xlsx("/home/zz2565/data/thesis_data/41588_2017_BFng3970_MOESM3_ESM.xlsx", sheet=10)
try1 <- dplyr::pull(try, 7)[-1]
try2 <- dplyr::pull(try, 8)[-1]
try3 <- dplyr::pull(try, 14)[-1]
try_d = data.frame(gene=try1, class=try2, pli=as.numeric(try3))
try_d <- na.omit(try_d)
try_d <- try_d[try_d$class!="syn"&try_d$pli>0.899999,]
try_d_scores <- subset(gene_scores, rownames(gene_scores)%in%unique(try_d$gene))

# subtype risk genes
subtg <- c("CHD7", "NOTCH1","FLT4", "KMT2D", "RBFOX2", "GATA6")
subtg_scores <- subset(gene_scores, rownames(gene_scores)%in%unique(subtg))

# write column names right
draw_colnames_rite <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.25 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(7,"bigpts"), vjust = 0.5, hjust = 1, rot = 0, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_rite",
                  ns=asNamespace("pheatmap"))

# marker gene enrichment
pheatmap(as.matrix(marker_genes_scores),
         scale="row",
         cellwidth=5,
         cellheight=5,
         fontsize = 4,
         filename="thesis_imgs/enrichment_imgs/1_marker_gene_enrichment.png")

# CHD risk gene enrichment
pheatmap((as.matrix(CHD_genes_scores)),
         scale="row",
         cellwidth=15,
         cellheight=2,
         fontsize_col=5,
         fontsize_row=2,
         fontsize=5,
         filename="thesis_imgs/enrichment_imgs/2_CHD_gene_enrichment.pdf")

  # de novo mutation enrichment
pheatmap(as.matrix(try_d_scores),
         scale="row",
         cellwidth=10,
         cellheight=2,
         fontsize=2,
         filename="thesis_imgs/enrichment_imgs/3_de_novo_gene_enrichment.pdf")
            
# subtg scores
pheatmap(as.matrix(subtg_scores),
         cluster_cols=F,
         cluster_rows=F,
         scale="row",
         cellwidth=25,
         cellheight=25,
         fontsize = 20,
         filename="thesis_imgs/enrichment_imgs/4_subtg_enrichment.png")                       
                        