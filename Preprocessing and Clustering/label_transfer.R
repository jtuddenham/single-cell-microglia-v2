library(caret)
library(e1071)
library(dplyr)
library(Matrix)
library(Seurat)	
library(ROSE)
library(tidyverse)
library(xgboost)
library(DropletUtils)

#label transfer pipeline#
source(clustering_diffex_functions)
#using the "uglia_only" object with pre-assigned IDs from the clustering/diffex pipeline#
DefaultAssay(uglia_only) = "RNA"
uglia_only$uglia_idents = Idents(uglia_only)
uglia_only[["SCT"]] = NULL
uglia_only <- UpdateSeuratObject(uglia_only)

##############
#mNN workflow#
##############

##we will use the DEGs from our pairwise diffex approach as the genes to train on between distinct clusters##
allDEG = read.csv("RNA_MAST_alldegenes.csv")
allcl = 1:12
feature_list = list()
for(i in allcl[-12]){
  for(j in (i+1):12){
    DEG_list = allDEG %>% filter(cluster == i & vs ==j) ##select the i vs. j pairs of DEG
    feature_list[[paste(i,j,sep = "v")]] = DEG_list$gene[abs(DEG_list$avg_logFC) > 0.25] ###and retain those with avg logFC b/w cluster >0.25 (arbitrary)
  }
}

#parameters for mNN merging
num_var_genes = length(unique(allDEG$gene))
pc_dimensions = 40
resolution = c(0.35, 0.5, 0.7)

#representative example- hasselmann et al. xenograft data - other datasets follow a similar approach#
iPSC_seurat = load_dataset("xeno")

#no need to downsample for this dataset- for v3 datasets, will need to downsample to equalize counts between original data, query dataset#
iPSC_seurat = CreateSeuratObject(counts = final_mat, min.cells = 3, names.delim = "-", names.field = 1, project = "uglia-complete")
iPSC_seurat$tech = "v2"

#merge datasets#
train_merge= merge(uglia_only, iPSC_seurat) %>% NormalizeData()
rm(iPSC_seurat); gc()
#use predetermined ordering from original mNN integration + add query data one at a time#
merge_order = c(unlist(uglia_only@tools$RunFastMNN@metadata$merge.info[111,1]), unlist(uglia_only@tools$RunFastMNN@metadata$merge.info[111,2]), 
                setdiff(unique(train_merge$orig.ident), unique(uglia_only$orig.ident)))
train_merge = mNN_pipeline(seurat = train_merge, merge_order = merge_order, DEG_genes = unique(allDEG$gene), num_var_genes = num_var_genes, pc_dimensions = pc_dimensions, resolution = resolution)

##############
#SVM workflow#
##############

uglia_only = train_merge

#separate out the query data#
uglia_only$SCT_snn_res.0.7[is.na(uglia_only$SCT_snn_res.0.7)] = 14
uglia_only$SCT_snn_res.0.7 = factor(uglia_only$SCT_snn_res.0.7, levels = c(0:14))
uglia_only = SetIdent(object = uglia_only, value = uglia_only$SCT_snn_res.0.7)
new.ids = 1:length(unique(Idents(uglia_only)))
names(new.ids) <- 0:(length(unique(uglia_only$SCT_snn_res.0.7))-1)
uglia_only <- RenameIdents(uglia_only, new.ids)
xeno_uglia = subset(uglia_only, idents = 15, invert = F)
###and now set up the old reference data###
uglia_only = subset(uglia_only, idents = c(1:12), invert = F)
new.ids = new.ids[1:12]

##randomly downsample data to 0.2x size while maintaining identical cluster proportions to reduce compute time##
clust_sizes = table(Idents(uglia_only))
prop_retained = floor(clust_sizes * 0.2)
retain_set = c()
validation_set = c()
holdout_set = c() ###just the leftovers###
metadata=as.character(Idents(uglia_only)) ###vector of cluster identities
whole_set = 1:length(metadata)

##instead, sample 50/50 from the smallest minority classes to improve accuracy on these##
prop_retained[10:12] = floor(clust_sizes[10:12]/2)

#set seed for reproducibility-create a training and validation set#
set.seed(length(metadata))
for(i in 1:length(prop_retained)){
  sampling_set = whole_set[metadata == i]
  retain_samp = sample(sampling_set, prop_retained[[i]], replace = F)
  validation_samp = sample(sampling_set[!(sampling_set %in% retain_samp)], prop_retained[[i]], replace = F)
  retain_set = c(retain_set,retain_samp)
  validation_set = c(validation_set, validation_samp)
}
# generate the validation and training matrices #
DefaultAssay(uglia_only) = "mnn.reconstructed"
validation_only = uglia_only[,validation_set]
train_only = uglia_only[,retain_set]
tpm_meta = as.numeric(Idents(train_only)) ###vector of cluster identities
vpm_meta = as.numeric(Idents(validation_only)) ###vector of cluster identities

#reformat the metadata based on high-level groupings#
overall_levels = c("1,6,7", "2,4", "3", "5", "8", "9", "10", "11", "12")
tpm_meta[tpm_meta %in% c(1,6:7)] = "1,6,7"
tpm_meta[tpm_meta %in% c(2,4)] = "2,4"
tpm_meta = factor(tpm_meta, levels = overall_levels)
vpm_meta[vpm_meta %in% c(1,6:7)] = "1,6,7"
vpm_meta[vpm_meta %in% c(2,4)] = "2,4"
vpm_meta = factor(vpm_meta, levels = overall_levels)

#as with mNN, we will use the DEGs from our pairwise diffex approach as the genes to train on between distinct clusters #
allDEG = read.csv("RNA_MAST_alldegenes.csv")
feature_list = list()
final_DEG_list = c()
for(i in 1:(length(overall_levels)-1)){
  for(j in (i+1):length(overall_levels)){
    DEG_list = allDEG %>% filter(cluster %in% strsplit(overall_levels[i], ",")[[1]] & vs %in% strsplit(overall_levels[j], ",")[[1]]) ##select the i vs. j pairs of DEG
    feature_list[[paste(overall_levels[i],overall_levels[j],sep = "v")]] = DEG_list$gene[abs(DEG_list$avg_logFC) > 0.25] ###and retain those with avg logFC b/w cluster >0.25 (arbitrary)
    final_DEG_list = union(final_DEG_list, feature_list[[paste(overall_levels[i],overall_levels[j],sep = "v")]])
  }
}

#for PCA#
#subset to only genes of interest#
VariableFeatures(validation_only) = final_DEG_list
validation_only = ScaleData(object = validation_only)
VariableFeatures(train_only) = final_DEG_list
train_only = ScaleData(object = train_only)
DefaultAssay(xeno_uglia) = "mnn.reconstructed"
VariableFeatures(xeno_uglia) = final_DEG_list
xeno_uglia = ScaleData(object = xeno_uglia)

#generate scaled data for input to PCA
tpm = t(train_only@assays$mnn.reconstructed@scale.data)
vpm = t(validation_only@assays$mnn.reconstructed@scale.data)
ipm = t(xeno_uglia@assays$mnn.reconstructed@scale.data)
ipm = ipm[,colnames(tpm)]

###set up some variables to specify how we construct the NN###
outcome_all = list()
##now loop through and examine each pair - this loop takes a long time and is memory intensive##
for(i in 1:(length(overall_levels)-1)){
  print(overall_levels[i])
  for(j in(i+1):length(overall_levels)){
    print(overall_levels[j])
    #version 1- PCA
    #T
    temp_tpm = tpm[tpm_meta %in% c(overall_levels[i], overall_levels[j]),
                   colnames(tpm) %in% feature_list[[paste(overall_levels[i],overall_levels[j],sep = "v")]]] ##reduce size of matrix and metadata to only those in this comparison
    temp_meta = tpm_meta[tpm_meta %in% c(overall_levels[i], overall_levels[j])]
    temp_meta = factor(temp_meta, levels = c(overall_levels[i], overall_levels[j]))
    #V- classifying the ENTIRE matrix each time to prepare for ensemble classification
    temp_vpm = vpm[,colnames(vpm) %in% feature_list[[paste(overall_levels[i],overall_levels[j],sep = "v")]]] ##reduce size of matrix and metadata to only those in this comparison
    temp_metaV = vpm_meta
    #G- classifying entire matrix each time for ensemble classification on the query dataset
    temp_ipm = ipm[,colnames(ipm) %in% feature_list[[paste(overall_levels[i],overall_levels[j],sep = "v")]]]
    
    #to help reduce class imbalance, perform over/under resampling#
    temp_df = data.frame(temp_tpm); temp_df$class = temp_meta
    rebalanced_data = ovun.sample(class~., temp_df, method = "both", p = 0.5)
    #perform SVM training with caret- 10-fold CV with PCA in each iteration#
    train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
    svm_rad <- train(x = as.matrix(rebalanced_data$data[,-dim(rebalanced_data$data)[2]], sparse= F), 
                     y = as.factor(rebalanced_data$data[,dim(rebalanced_data$data)[2]]),
                     method = "svmRadial", 
                     trControl = train_control,  
                     preProcess = c("pca"),
                     tunelength = 10)
    
    #now classify holdout, ipm data#
    test_vpm = as.matrix(temp_vpm[temp_metaV %in% c(overall_levels[i], overall_levels[j]),])
    test_meta = factor(temp_metaV[temp_metaV %in% c(overall_levels[i], overall_levels[j])], levels = c(overall_levels[i], overall_levels[j]))
    test_pred_rad <- predict(svm_rad, newdata = as.matrix(data.frame(test_vpm)))
    print(confusionMatrix(test_pred_rad, test_meta))
    all_test_pred <- predict(svm_rad, newdata = as.matrix(data.frame(temp_vpm)))
    ipm_pred <- predict(svm_rad, newdata = as.matrix(data.frame(temp_ipm)))
    
    outcome_all[[paste(i,j,sep = "v")]] = list(model = svm_rad, test_pred = test_pred_rad, test_meta = temp_metaV, all_test = all_test_pred, 
                                               all_test_meta = temp_metaV, ipm_pred = ipm_pred)
  }
}
save(outcome_all, file = "../intermediate_data/2022_SVM_pw_xeno_PCAcomponents.rda")

########################
####XGBoost workflow####
########################

## we will use the DEGs from our pairwise diffex approach as the genes to train on between distinct clusters ##
allDEG = read.csv("RNA_MAST_alldegenes.csv")
allcl = 1:12
feature_list = list()
for(i in allcl[-12]){
  for(j in (i+1):12){
    DEG_list = allDEG %>% filter(cluster == i & vs ==j) ##select the i vs. j pairs of DEG
    feature_list[[paste(i,j,sep = "v")]] = DEG_list$gene[abs(DEG_list$avg_logFC) > 0.25] ###and retain those with avg logFC b/w cluster >0.25 (arbitrary)
  }
}

#re-using the same "uglia_only" and "xeno_uglia" objects separated during the SVM workflow#
##randomly downsample data to 0.3x size while maintaining identical cluster proportions to reduce compute time##
clust_sizes = table(Idents(uglia_only))
prop_retained = floor(clust_sizes / 3)
retain_set = c()
validation_set = c()
holdout_set = c() ###just the leftovers###
metadata=as.character(Idents(uglia_only)) ###vector of cluster identities
whole_set = 1:length(metadata)

###optional-instead, fully drop the last 3 groupings###
start = 8
prop_retained[start:12] = 0

###set seed for reproducibility - generate train/test/holdout###
set.seed(length(metadata))
for(z in 1:length(prop_retained)){
  if(z %in% c(start:12)){
    next
  }
  sampling_set = whole_set[metadata == z]
  retain_samp = sample(sampling_set, prop_retained[[z]], replace = F)
  validation_samp = sample(sampling_set[!(sampling_set %in% retain_samp)], prop_retained[[z]], replace = F)
  holdout_samp = sampling_set[!((sampling_set %in% retain_samp) | (sampling_set %in% validation_samp))] ###just for reference, here are the leftovers###
  print(isTRUE(sum(length(holdout_samp), length(retain_samp), length(validation_samp)) == clust_sizes[[z]])) #just check to make sure you're accounting for all cells
  retain_set = c(retain_set,retain_samp)
  validation_set = c(validation_set, validation_samp)
  holdout_set = c(holdout_set,holdout_samp)
}
sample_sets = list(retain_set, validation_set, holdout_set)

#re-prime the query object#
DefaultAssay(xeno_uglia) = "mnn.reconstructed"
VariableFeatures(xeno_uglia) = unique(allDEG$gene)
xeno_uglia = ScaleData(object = xeno_uglia)

#train for each of the 3 partitions#
for(j in 1:3){
  print(j)
  train_set = sample_sets[[j]]
  val_set = unlist(sample_sets[-j])
  
  ##set up the train and validation objects## 
  validation_only = uglia_only[,val_set]
  train_only = uglia_only[,train_set]
  
  ###initialize some relevant variables and clean up
  tpm_meta = as.numeric(Idents(train_only)) ###vector of cluster identities
  vpm_meta = as.numeric(Idents(validation_only)) ###vector of cluster identities
  
  #reformat the metadata based on high-level groupings#
  overall_levels = c("167", "24", "3", "5", "8", "9", "10", "11", "12")
  tpm_meta[tpm_meta %in% c(1,6:7)] = "167"
  tpm_meta[tpm_meta %in% c(2,4)] = "24"
  tpm_meta = factor(tpm_meta, levels = overall_levels)
  vpm_meta[vpm_meta %in% c(1,6:7)] = "167"
  vpm_meta[vpm_meta %in% c(2,4)] = "24"
  vpm_meta = factor(vpm_meta, levels = overall_levels)
  
  #flat xgboost classifier#
  train_meta = as.integer(tpm_meta) - 1 ###xgboost requires the first class to be 0, and needs factors to be converted to integers
  test_meta = as.integer(vpm_meta) - 1
  
  #subset to only genes of interest#
  VariableFeatures(validation_only) = unique(allDEG$gene)
  validation_only = ScaleData(object = validation_only)
  VariableFeatures(train_only) = unique(allDEG$gene)
  train_only = ScaleData(object = train_only)
  
  ##generate PCA as input to the XGBoost model
  train_tpm = t(train_only@assays$mnn.reconstructed@scale.data)
  test_vpm = t(validation_only@assays$mnn.reconstructed@scale.data)
  
  #run PCA and assemble datasets for training#
  num_pcs = 50
  train_pca <- prcomp(x = train_tpm, scale. = F, center = F, rank. = num_pcs)
  pc_mat_train = train_pca$x[,1:num_pcs]
  pc_mat_test = test_vpm %*% train_pca$rotation[,1:num_pcs]
  
  ###load in the xeno data###
  ipm = t(xeno_uglia@assays$mnn.reconstructed@scale.data)
  ipm = ipm[,colnames(train_tpm)]
  projected_ipm = ipm %*% train_pca$rotation[,1:num_pcs]
  
  ##set up XGboost matrices###
  xgb.train = xgb.DMatrix(data=pc_mat_train,label=train_meta)
  xgb.test = xgb.DMatrix(data=pc_mat_test,label=test_meta)
  xgb.ibm = xgb.DMatrix(data = projected_ipm)
  num_class = length(levels(tpm_meta))
  rm(train_only, validation_only); gc()
  
  ###model optimization with caret##
  xgb_trcontrol = trainControl(
    method = "cv",
    number = 5,  
    allowParallel = TRUE,
    verboseIter = FALSE,
    returnData = FALSE
  )
  ##choose nrounds##
  params <- list(booster = "gbtree", 
                 objective = "multi:softprob", eval_metric="mlogloss",
                 eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample=0.8, colsample_bytree=0.8, num_class = num_class)
  xgbcv <- xgb.cv( params = params, data = xgb.train, nrounds = 300, nfold = 5, showsd = T, stratified = T, print.every.n = 10, early.stop.round = 20, maximize = F, watchlist=list(val1=xgb.train,val2=xgb.test))
  
  xgbGrid <- expand.grid(nrounds = xgbcv$best_iteration,  
                         max_depth = c(3,5,7,9),
                         colsample_bytree = 0.8,
                         eta = 0.3,
                         gamma = 0.1,
                         min_child_weight = seq(5,11,2),
                         subsample = 0.8
  )
  set.seed(0)
  xgb_model = train(
    xgb.train, train_meta,  
    trControl = xgb_trcontrol,
    tuneGrid = xgbGrid,
    method = "xgbTree"
  )
  
  ##now retune new params - gamma this time##
  xgbGrid <- expand.grid(nrounds = xgbcv$best_iteration,
                         max_depth = xgb_model$bestTune$max_depth,
                         colsample_bytree = 0.8,
                         eta = 0.3,
                         gamma = seq(0.4, 1.2, 0.2),
                         min_child_weight = xgb_model$bestTune$min_child_weight,
                         subsample = 0.8
  )
  set.seed(0)
  xgb_model = train(
    xgb.train, train_meta,  
    trControl = xgb_trcontrol,
    tuneGrid = xgbGrid,
    method = "xgbTree"
  )
  ##recalculate number of rounds##
  params <- list(booster = "gbtree", 
                 objective = "multi:softprob", eval_metric="mlogloss",
                 eta=0.3, gamma=xgb_model$bestTune$gamma, max_depth=xgb_model$bestTune$max_depth, 
                 min_child_weight=xgb_model$bestTune$min_child_weight, subsample=0.8, colsample_bytree=0.8, num_class = num_class)
  xgbcv <- xgb.cv(params = params, data = xgb.train,
                  nrounds = 300, nfold = 5, showsd = T, stratified = T,
                  print_every_n = 10, early_stopping_rounds = 20, maximize = F, #scale_pos_weight = sqrt(table(temp_meta)[[1]]/table(temp_meta)[[2]]),
                  watchlist=list(val1=xgb.train,val2=xgb.test))
  ##now retune new params - last 2  this time##
  xgbGrid <- expand.grid(nrounds = xgbcv$best_iteration,
                         max_depth = xgb_model$bestTune$max_depth,
                         colsample_bytree = seq(0.7, 1.2, 0.1),
                         eta = 0.3,
                         gamma = xgb_model$bestTune$gamma,
                         min_child_weight = xgb_model$bestTune$min_child_weight,
                         subsample = seq(0.7, 1.2, 0.1)
  )
  set.seed(0)
  xgb_model = train(
    xgb.train, train_meta,  
    trControl = xgb_trcontrol,
    tuneGrid = xgbGrid,
    method = "xgbTree"
  )
  
  params <- list(booster = "gbtree", objective = "multi:softprob", eval_metric="mlogloss",
                 eta=0.3, gamma=xgb_model$bestTune$gamma, max_depth=xgb_model$bestTune$max_depth, 
                 min_child_weight=xgb_model$bestTune$min_child_weight, subsample=xgb_model$bestTune$subsample, 
                 colsample_bytree=xgb_model$bestTune$colsample_bytree, num_class = num_class)
  xgb_final = xgb.train(params = params, data = xgb.train, nrounds = xgb_model$bestTune$nrounds, print_every_n = 10, 
                        early_stopping_rounds = 20, maximize = F, watchlist=list(val1=xgb.train,val2=xgb.test))
  
  #predict test data#
  xgbpred <- predict(xgb_final, xgb.test, reshape = T)
  colnames(xgbpred) = overall_levels
  predicted_cluster = colnames(xgbpred)[max.col(xgbpred,ties.method="first")]
  names(predicted_cluster) = names(vpm_meta)

  ####classify and analyze the query ###
  ipm_pred <- predict(xgb_final, xgb.ibm, reshape = T)
  colnames(ipm_pred) = overall_levels
  i_predicted = colnames(ipm_pred)[max.col(ipm_pred,ties.method="first")]
  names(i_predicted) = rownames(projected_ipm)
  igb_pred_val =  apply(X = ipm_pred, MARGIN = 1, FUN = max)
  names(igb_pred_val) = rownames(projected_ipm)
  i_predicted = factor(i_predicted, levels = overall_levels)
  print(table(i_predicted))

  ###save the outputs
  xgb_outs = list(final_model = xgb_final, val_pred = xgbpred, query_pred = xgbpred, vpm_meta = vpm_meta,
                  val_set = val_set, train_set = train_set, PCA = train_pca, sample_sets = sample_sets, ipm_pred = ipm_pred
  )
  save(xgb_outs, file = paste0("../intermediate_data/mNN_XGB_xeno_", j, "_subset_PCA.rda"))
}

#####################
#merging predictions#
#####################
load(paste0("2020_Analysis_Outs/classification_final/", i, "/", i, "_mNNPCA_SVM.rda"))
load(paste0("2020_Analysis_Outs/classification_final/", i,  "/mNN_XGB_", i, "_1_subset_PCA.rda"))

##final SVM predictions for query dataset##
ipm_mat = matrix(ncol = length(unique(outcome_all$`1v2`$test_meta)), nrow = length(outcome_all[[1]]$ipm_pred), data = 0)
colnames(ipm_mat) = unique(outcome_all$`1v2`$test_meta)
rownames(ipm_mat) = rownames(outcome_all[[1]]$ipm_pred)
for(z in 1:length(outcome_all)){
  acc = 1 
  model_labels=sparseMatrix(i=1:length(outcome_all[[z]]$ipm_pred),j=sapply(as.character(outcome_all[[z]]$ipm_pred), FUN = function(x){switch(x, "1,6,7" = 1, "2,4" = 2, "3" = 3, "5"=4, "8"=5, "9"=6, "10"=7, "11"=8, "12"=9)}),
                            dims = c(length(outcome_all[[z]]$ipm_pred), length(unique(outcome_all$`1v2`$test_meta))), x = 1*acc) ###could scale the X by the efficacy of the model if we are going by raw votes
  ipm_mat = ipm_mat + model_labels
}
svm_ipm_pred = factor(colnames(ipm_mat)[max.col(ipm_mat,ties.method="first")], levels = overall_levels)

#final XGBoost predictions for query dataset#
i_predicted = factor(colnames(xgb_outs$ipm_pred)[max.col(xgb_outs$ipm_pred,ties.method="first")], labels = overall_levels)
names(i_predicted) = colnames(iPSC_seurat)
igb_pred_val =  apply(X = xgb_outs$ipm_pred, MARGIN = 1, FUN = max) #provides the XGBoost confidence for the query

#now, merge predictions from the two.
full_ipred = as.integer(i_predicted)
svm_ipred = as.integer(svm_ipm_pred)
full_ipred[cells_101112] = svm_ipred[cells_101112]

#remove the low-confidence cells classified via XGBoost# 
full_ipred = full_ipred[igb_pred_val > confidence_thresh | cells_101112]
full_ipred = factor(full_ipred, labels = overall_levels)
names(full_ipred) = colnames(iPSC_seurat)
#final predictions#
print(table(full_ipred))