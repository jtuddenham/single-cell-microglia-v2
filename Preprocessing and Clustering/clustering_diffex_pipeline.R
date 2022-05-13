#loading in required packages#
library(Matrix)
library(Seurat)
library(dplyr)
library(readxl)
library(SeuratWrappers)
library(batchelor)
library(DropletUtils)
source("clustering_diffex_functions.R")

##initial processing and merging of data##
allfiles = list.files("../raw_data/")
sampdata=read_excel("../intermediate_data/Table S1 - Overview of demographics, hashing strategy, and quality control, related to STAR Methods..xlsx", sheet = 1)
featuredata = read_excel("../intermediate_data/Table S1 - Overview of demographics, hashing strategy, and quality control, related to STAR Methods..xlsx", sheet = 3)

#set up sample lists for aggregation#
allfiles_hash = allfiles[allfiles %in% featuredata$`Sequencing ID`]
allfiles_nohash = allfiles[!(allfiles %in% allfiles_hash)]

#make full data objects#
dat_barcodes <- aggfilter_nonhash(allfiles_nohash, base_directory = "../raw_data/", umimin = 500, umimax = 10000, RP_filt = T, MT_filt = T)
alldat = dat_barcodes[[1]]; barcodes = dat_barcodes[[2]]
hash_counts <- aggfilter_hash(allfiles_hash, base_directory = "../raw_data/", umimin = 500, umimax = 10000, 
                              HTO_filter = 0, feature_sheet = "../intermediate_data/Table S1 - Overview of demographics, hashing strategy, and quality control, related to STAR Methods..xlsx",
                              RP_filt = T, MT_filt = T)
allhash = hash_counts[[1]]; counts = hash_counts[[2]]

#generate the combined seurat object with added metadata and finalized hash assignments#
uglia.seurat <- makeseurat_metadata(alldat = alldat, allhash = allhash, sampdata = sampdata, featuredata = featuredata)

#inspection of general range for number of PCs by way of elbow plot#
uglia.seurat<-NormalizeData(uglia.seurat)
uglia.seurat <-FindVariableFeatures(uglia.seurat, nfeatures = 6000)
uglia.seurat <- ScaleData(uglia.seurat)
uglia.seurat <- RunPCA(uglia.seurat, features = VariableFeatures(object = uglia.seurat))
ElbowPlot(uglia.seurat, ndims = 50)

#run top-level clustering - this step takes a long time and a large amount of memory#
var_genes = 4500
PCs = 40
rez = c(0.35, 0.5, 0.75)
full.combined = sct_mnn(uglia.seurat, regress = c("nCount_RNA", "tech"), num_var_genes = var_genes, resolution = rez, pc_dimensions = PCs, auto_merge = T, conserve.memory = T)

#using top-level cluster assignments, perform subclustering of microglia#
assigned_idents = read.csv("../intermediate_data/4.5k40_batch_SCTv3_0.5_idents.csv")
names = assigned_idents[,1]
assigned_idents = data.frame(assigned_idents[,2])
row.names(assigned_idents) = names
uglia.seurat = AddMetaData(object = uglia.seurat, metadata = assigned_idents, col.name = "assigned_idents")
uglia.seurat = SetIdent(object = uglia.seurat, value = uglia.seurat$assigned_idents)

##get microglia
uglia_only = subset(uglia.seurat, idents = c(0:8, 10:13, 15))
###now downsample the v3 count matrix###
v3_only = subset(uglia_only, cells = names(uglia_only$tech)[uglia_only$tech == "v3"])
v2_only = subset(uglia_only, cells = names(uglia_only$tech)[uglia_only$tech == "v2"])
count_mat = v3_only[["RNA"]]@counts
count_mat = downsampleMatrix(count_mat, prop = 0.5)
v3_only = CreateSeuratObject(count_mat, meta.data = v3_only@meta.data)
v2_only = CreateSeuratObject(v2_only[["RNA"]]@counts, meta.data = v2_only@meta.data)
uglia_only = merge(v3_only, v2_only)
uglia_only = NormalizeData(uglia_only)
uglia_list = SplitObject(uglia_only, split.by = "orig.ident")
rm(uglia_only, v2_only, v3_only, count_mat); gc()

#perform microglial subclustering - this step takes a long time and a large amount of memory#
num_var_genes = 4500
pc_dimensions = 40
resolution = c(0.35, 0.5, 0.7)
uglia_list = lapply(X=uglia_list, FUN = function(x){SCTransform(x, do.correct.umi = T, return.only.var.genes = T,
                                                                variable.features.n = num_var_genes, verbose = T, do.center = T, conserve.memory = T) })
uglia_only <- RunFastMNN(object.list = uglia_list, features = num_var_genes, auto.merge = T)
uglia_only <- RunUMAP(uglia_only, reduction = "mnn", dims = 1:pc_dimensions)
uglia_only <- FindNeighbors(uglia_only, reduction = "mnn", dims = 1:pc_dimensions)
uglia_only <- FindClusters(uglia_only, resolution = resolution)

#using pre-calculated IDs#
uglia_idents = read.csv("../intermediate_data/SCT_batch_downsamp_split_0.7_idents.csv")
names = uglia_idents[,1]
uglia_idents = data.frame(uglia_idents[,2])
row.names(uglia_idents) = names
uglia_only = AddMetaData(object = uglia_only, metadata = uglia_idents, col.name = "uglia_idents")
uglia_only$uglia_idents = factor(uglia_only$uglia_idents, levels = c(0:max(unique(uglia_idents))))
uglia_only = SetIdent(object = uglia_only, value = uglia_only$uglia_idents)
new.ids = 1:length(unique(uglia_only$uglia_idents))
names(new.ids) <- 0:(length(unique(uglia_only$uglia_idents))-1)
uglia_only <- RenameIdents(uglia_only, new.ids)

##generate the DE genes for microglia using MAST-based differential expression testing##
out_stub = "../intermediate_data/pairwise_diffex/downsampled/split/"
all_pairwise_diffex(uglia_only, outdir = paste(out_stub, "RNA_MAST/RNA_MAST_", sep =""),
                    assay = "RNA", slot = "data", test = "MAST")
microglial_pops = c(1:12)
markers_up = pairwise_genes(DEG_file = paste(out_stub, "RNA_MAST/RNA_MAST_alldegenes.csv", sep =""), allpops = microglial_pops, upreg = T)
write.csv(markers_up, paste(out_stub, "RNA_MAST/microglia_pairwise_up.csv", sep = ""))
markers_down = pairwise_genes(DEG_file = paste(out_stub, "RNA_MAST/RNA_MAST_alldegenes.csv", sep =""), allpops = microglial_pops, upreg = F)
write.csv(markers_down, paste(out_stub, "RNA_MAST/microglia_pairwise_down.csv", sep = ""))
