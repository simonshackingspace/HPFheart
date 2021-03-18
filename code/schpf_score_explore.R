library(readxl)
library(dplyr)
require(ROSE)

# get data matrix 
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

# revise mutation type
mut.redef <- function(df) {
  # df[df[,"Variant_Class"]=="misD","Variant_Class"] = "mis"
  df[(df[,"Variant_Class"]=="startloss") | (df[,"Variant_Class"]=="stoploss"),"Variant_Class"] = "non"
  return(df)
}

# get num of de novo mut
genes.dnv <- function(vec, df, lab) {
  if (lab=="ctrl") {
    rdf <- df[rownames(df) %in% vec, ]
    rdf$label <- lab
    return(rdf)
  }
  genes.table <- table(vec)
  genes_multi_dnv <- names(which(genes.table>1))
  rdf <- df[rownames(df) %in% genes_multi_dnv, ]
  rdf$label <- lab
  return(rdf)
}

# schpf genes and scores
genes <- read.table("/home/zz2565/data/thesis_data/diff_seqs/scHPF_pretrain/genes.txt")[,2]
genes.scores <- read.csv("thesis_output/scHPF_best_score/gene_score.txt", sep='\t', header=F)
rownames(genes.scores) <- genes

# dnv info
denovo.ctrl <- read_xlsx("/home/zz2565/data/thesis_data/41588_2017_BFng3970_MOESM3_ESM.xlsx", sheet=10)
denovo.case <- read_xlsx("/home/zz2565/data/thesis_data/41588_2017_BFng3970_MOESM3_ESM.xlsx", sheet=9)
denovo.ctrl <- mut.redef(as.data.frame(header.true(denovo.ctrl), stringAsFactors=T))
denovo.case <- mut.redef(as.data.frame(header.true(denovo.case), stringAsFactors=T))
# lof dnv
lof_muts <- c("non", "frameshift", "splice")
denovo.ctrl_lof <- denovo.ctrl[denovo.ctrl$Variant_Class %in% lof_muts, "Gene"]
# go for multiple denovo.ctrl_lof <- 
denovo.case_lof <- denovo.case[denovo.case$Variant_Class %in% lof_muts, "Gene"]

# for classification
genes.case <- genes.dnv(denovo.case_lof, genes.scores, "case")
genes.ctrl <- genes.dnv(denovo.ctrl_lof, genes.scores, "ctrl")
genes.intersect <- intersect(rownames(genes.ctrl), rownames(genes.case))
genes.all <- rbind(genes.case, genes.ctrl)
genes.all <- genes.all[!(rownames(genes.all) %in% genes.intersect),]
genes.all$label <- as.factor(genes.all$label)

smp_siz <- floor(0.75*nrow(genes.all))
set.seed(123)   # set seed to ensure you always have same random numbers generated
train_ind <- sample(seq_len(nrow(genes.all)),size = smp_siz)  # Randomly identifies therows equal to sample size ( defined in previous instruction) from  all the rows of Smarket dataset and stores the row number in train_ind
train <- genes.all[train_ind,] #creates the training dataset with row numbers stored in train_ind
test <- genes.all[-train_ind,]  # creates the test dataset excluding the row numbers mentioned in train_ind

train_balanced_both <- ovun.sample(label ~ ., data = train, method = "both", p=0.5, N=1000)$data     

model <- glm(formula=label~., family="binomial", data=train_balanced_both)
summary(model)
res <- predict(model, newdata=test, type="response")

pred_true <- names(which(res > 0.5))
real_true <- rownames(test[test$label=="case", ])
