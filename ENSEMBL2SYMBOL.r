library(clusterProfiler)
library(org.Hs.eg.db)

remove_duplicates <- function(df, column_name) {
  df$sum <- rowSums(df[, -which(names(df) == column_name)])
  df <- df[order(df[[column_name]], -df$sum), ]
  df <- df[!duplicated(df[[column_name]]), ]
  df$sum <- NULL
  return(df)
}

disease = 'PAAD'
data_type = 'tpm_unstranded'#unstranded;tpm_unstranded


input = paste0("/home/liuchenglong/Documents/TCGA/3.intergra/",disease,"_tpm_unstranded.txt")
genetype = 'ensembl_gene_id'
data = read.csv(input,sep = '\t',header=TRUE,check.names = FALSE)
result = bitr(data$gene_id,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
colnames(result) = c('gene_id','symbol')
data1 = as.data.frame(merge(data,result,by = 'gene_id'))
data1 = data1[data1$symbol !='',]
data1 = data1[,-1]
data1 = remove_duplicates(data1,'symbol')
rownames(data1) = data1$symbol
data1 = data1[,-ncol(data1)]
rm(data)
colnames(data1) <- substr(colnames(data1), 1, 12)
write.table(data1,paste0('/home/liuchenglong/Documents/TCGA/4.ENSEMBL2SYMBOL/',disease,'_tpm_unstranded.txt'),sep = '\t',row.names = TRUE)