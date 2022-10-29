#printing dataset
data_file = read.csv("C:/Users/Dell !/OneDrive/Documents/MSc sem 3/cancer genomics/GSE214779_Raw_gene_counts_matrix_METASTASIS.csv",sep = ",",header = T,row.names = 1)
print(data_file)

#create count per matrix
cpmatrix=data_file 
for (i in 1:ncol(data_file)) {
  cpmatrix[,i]=(data_file[,i]/sum(data_file[,i]))*1000000
}
cpmatrix[is.na(cpmatrix)]=0

#log
logcpm = log2(cpmatrix+1)
summary(logcpm)

#zscore
library(matrixStats)
z_score = (logcpm-rowMeans(logcpm))/rowSds(as.matrix(logcpm))[row(logcpm)]
print(z_score)

#variance
vargenes = apply(z_score,1,var)
print(vargenes)

vargenes = sort(vargenes,decreasing = T)
top50 = vargenes[1:50]
pmat = z_score[names(top50),]
m=as.matrix(pmat)

#heatmap
library(ComplexHeatmap)
heatmap(pmat)
