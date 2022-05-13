#bootstrap resampling#
library(Matrix)
library(Seurat)
library(dplyr)
library(readxl)
library(SeuratWrappers)
library(batchelor)

#starting with "uglia_only" object as generated in the clustering_diffex_pipeline#
uglia_only <- subset(uglia_only, idents = 13, invert = TRUE)

###bootstrap clustering - repeat clustering on randomly selected subsets of cells to assess cluster robustness###
###this step can take a long time, so adjust the number of iterations as needed###
subset_fraction=0.75
num_iterations=100
#clustering pipeline parameters
num_var_genes = 4500
pc_dimensions = 40
resolution = 0.5

###set up sampling matrix###
sampmat=matrix(0,nrow=length(Idents(uglia_only)),ncol=num_iterations)
rowfrac=round(subset_fraction*ncol(sampmat))
allinds=rep(0,rowfrac*length(Idents(uglia_only)))
for (ii in 1:nrow(sampmat)) {
  startind=((ii-1)*rowfrac)+1
  endind=(ii*rowfrac)
  set.seed(ii)
  allinds[startind:endind]=sample(1:num_iterations,rowfrac)
}
rowinds=rep(1:nrow(sampmat),each=rowfrac)
sampmat[cbind(rowinds,allinds)]=1
rownames(sampmat)=names(Idents(uglia_only))

clmat=list()
for (ii in 1:num_iterations) {
  keep.cells=colnames(uglia_only)[which(sampmat[,ii]==1)]
  subset=subset(uglia_only,cells = keep.cells)
  subset <- FindNeighbors(subset, reduction = "mnn", dims = 1:pc_dimensions)
  subset <- FindClusters(subset, resolution = resolution)
  clmat[[ii]]=Idents(subset) ###may need to change this
  print(ii)
}
outmat2=matrix(NA,nrow=nrow(sampmat),ncol=ncol(sampmat))
rownames(outmat2)=names(Idents(uglia_only))
for (ii in 1:num_iterations) {
  outmat2[names(clmat[[ii]]),ii]=clmat[[ii]]
}

###assess cluster robustness using bootstrap clustering results###
topcl=as.numeric(as.character(Idents(uglia_only)))
allcl=unique(topcl)
clval=allcl[1]
coclustmat=matrix(0,nrow=length(allcl),ncol=length(allcl))
rownames(coclustmat)=allcl[order(allcl)]
colnames(coclustmat)=rownames(coclustmat)

####generate "adjacency scores" from averages of paired cell classifications####
for (ii in 1:12) {
  keepcells1=which(topcl==rownames(coclustmat)[ii])
  for (jj in ii:nrow(coclustmat)) {
    if (ii!=jj) {
      keepcells2 = which(topcl == rownames(coclustmat)[jj])
      tmpmat = outmat2[keepcells1,]
      overall_mean = 0
      final_pairs = 0
      for(i in 1:(nrow(tmpmat)-1)){
        overall_mean = overall_mean + mean(apply(outmat2[keepcells2,], 1, FUN = function(x){
          sum(tmpmat[i,] == x, na.rm = T)/sum(!is.na(tmpmat[i,]) & !is.na(x))}))
        final_pairs = final_pairs + 1
        print(final_pairs)
      }
      coclustmat[ii,jj] = overall_mean/final_pairs
    } else {
      tmpmat = outmat2[keepcells1,]
      overall_mean = 0
      final_pairs = 0
      for(i in 1:(nrow(tmpmat)-2)){
        comparison_cell = tmpmat[i,]
        overall_mean = overall_mean + mean(apply(tmpmat[-c(1:i),], 1, FUN = function(x){
          sum(comparison_cell == x, na.rm = T)/sum(!is.na(tmpmat[i,]) & !is.na(x))}))
        final_pairs = final_pairs + 1
        print(final_pairs)
      }
      coclustmat[ii,jj] = overall_mean/final_pairs
    }
    print(c(ii,jj,coclustmat[ii,jj]))
  }
}