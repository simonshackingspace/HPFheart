library(readxl)
library(plyr)
library(dplyr)
library(denovolyzeR)

##################################################### functions #######################################################
# get data matrix 
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

# revise mutation type
mut.redef <- function(df) {
  df[df$Variant_Class=="misD","Variant_Class"] = "mis"
  df[(df$Variant_Class=="startloss") | (df$Variant_Class=="stoploss"),"Variant_Class"] = "non"
  return(df)
}


# get num of de novo mut
mut.num <- function(df1, df2) {
  num_c <- c()
  for (i in 1:ncol(df1)){
    num_c <- c(num_c, sum(df1[,i] %in% df2$Gene))
  }
  return(num_c)
}

# get ready for denovolyeR
mut.dnv <- function(vec, df) {
  rdf <- df[df$Gene %in% vec, c("Gene", "Variant_Class")]
  colnames(rdf) <- c("gene", "class")
  return(rdf)
}

# case control study for diff factors
mut.ccs <- function(df1, df2, num) {
  vec.ccs <- vector(mode = "list", length = ncol(df1))
  for (i in 1:ncol(df1)) {
    dnv <- mut.dnv(df1[,i], df2)
    vec.ccs[[i]] <- denovolyzeByClass(genes=dnv$gene, 
                      classes=dnv$class,
                      nsamples=num, 
                      includeGenes=df1[,i])
  }
  return(vec.ccs)
}

# change subtypes names
st.rnm <- function(dt) {
  st <- dt$`Cardiac Category`
  st <- mapvalues(st, from=c("CTD (TGA)", "other", "other (AVC)", "Other (AVC)", "Other"), 
                  to=c("D-TGA", "OTHER", "OTHER", "OTHER", "OTHER"))
  return(st)
}

# case control for diff subtypes
st.ccs <- function(dt, stn){
  vec.ccs <- vector(mode = "list", length = nrow(stn))
  for (item in unique(dt$`Cardiac Category`)){
    dti <- dt[dt$`Cardiac Category`==item,]
    dnv_p <- denovolyzeByGene(genes=dti$Gene,
                            classes=dti$Variant_Class,
                            nsamples=stn[item,])
    # head(dnv_p[order(dnv_p$lof_pValue),])
    dnv_un <- dnv_p[dnv_p$lof_pValue < 0.05 | dnv_p$prot_pValue < 0.05,]
    vec.ccs[[item]] <- dnv_un[order(dnv_un$lof_pValue),]
  }
  return(vec.ccs)
}
#################################################### functions ########################################################

#################################################### analysis #########################################################
# laod data
top.genes <- read.csv("thesis_output/scHPF_best_score/ranked_genes.txt", sep='\t', header=F)[1:500,]
denovo.ctrl <- read_xlsx("/home/zz2565/data/thesis_data/41588_2017_BFng3970_MOESM3_ESM.xlsx", sheet=10)
denovo.case <- read_xlsx("/home/zz2565/data/thesis_data/41588_2017_BFng3970_MOESM3_ESM.xlsx", sheet=9)
denovo.ctrl <- mut.redef(as.data.frame(header.true(denovo.ctrl), stringAsFactors=T))
denovo.case <- mut.redef(as.data.frame(header.true(denovo.case), stringAsFactors=T))
# denovo.ctrl.hr <- denovo.ctrl[denovo.ctrl$`pLI Score` >= 0.9 & denovo.ctrl$Variant_Class != "syn", ]
# denovo.case.hr <- denovo.case[denovo.case$`pLI score` >= 0.9 & denovo.case$Variant_Class != "syn", ]

# count number of mutation
# num_mut_ctrl <- mut.num(top.genes, denovo.ctrl)
# num_mut_case <- mut.num(top.genes, denovo.case)
# num_mut_ctrl_hr <- mut.num(top.genes, denovo.ctrl.hr)
# num_mut_case_hr <- mut.num(top.genes, denovo.case.hr)

# general case control study
gccs.case <- denovolyzeByClass(genes=denovo.case$Gene, 
                  classes=denovo.case$Variant_Class,
                  nsamples=2645)
write.csv(gccs.case, "thesis_output/gccs.case.csv", row.names=F)
gccs.ctrl <- denovolyzeByClass(genes=denovo.ctrl$Gene, 
                               classes=denovo.ctrl$Variant_Class,
                               nsamples=1789)
write.csv(gccs.ctrl, "thesis_output/gccs.ctrl.csv", row.names=F)

# case control study for each factor
ccs.ctrl <- mut.ccs(top.genes, denovo.ctrl, 1789)
ccs.case <- mut.ccs(top.genes, denovo.case, 2645)
for (i in 1:length(ccs.ctrl)){
  sfx <- paste(i, "csv", sep=".")
  write.csv(ccs.ctrl[[i]], paste("thesis_output/ccs.ctrl", sfx, sep="_"), row.names=F)
  write.csv(ccs.case[[i]], paste("thesis_output/ccs.case", sfx, sep="_"), row.names=F)
}

  # case control with subtypes
denovo.case$`Cardiac Category` <- st.rnm(denovo.case)
subtype_nums <- data.frame(No=c(872, 251, 272, 797, 679),
                           row.names=c("CTD", "D-TGA", "HTX", "LVO", "OTHER"))
ccs.subt <- st.ccs(denovo.case, subtype_nums)

# find important genes for factor 1
dnv <- mut.dnv(top.genes[,4], denovo.case)
dnv_p <- denovolyzeByGene(genes=dnv$gene,
                          classes=dnv$class,
                          nsamples=2645)
head(dnv_p[order(dnv_p$lof_pValue),])
dnv_un <- dnv_p[dnv_p$lof_pValue < 0.05 | dnv_p$prot_pValue < 0.05,]
dnv_un[order(dnv_un$lof_pValue),]

write.csv(dnv_un$gene, "unexpected.txt", row.names=FALSE)
write.csv(unique(dnv$gene), "dnv.csv", row.names=F)
