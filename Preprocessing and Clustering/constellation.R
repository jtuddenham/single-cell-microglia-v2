#constellation diagram#
library(keras)
library(dplyr)
library(Matrix)
library(Seurat)	

#NN post-hoc cluster stability evaluation#
#starting with "uglia_only" object as generated in the clustering_diffex_pipeline using the pre-identified labels##
uglia_only <- subset(uglia_only, idents = 13, invert = TRUE)

#dataset is too big for NN training using all features- so instead, choose the top 10k variable genes#
DefaultAssay(uglia_only) = "RNA"
uglia_only <- NormalizeData(uglia_only, assay = "RNA", normalization.method = "LogNormalize")
uglia_only = FindVariableFeatures(uglia_only, assay = "RNA", nfeatures = 10000)
var_features = VariableFeatures(uglia_only)
#filter to retain only the most variable genes#
tpm = t(uglia_only@assays$RNA@data) ##normalized count data for each cell
tpm = tpm[,colnames(tpm) %in% var_features]

metadata=as.character(Idents(uglia_only)) ###vector of cluster identities
allcl=sort(as.numeric(unique(metadata))) ###unique cluster IDs
rm(uglia_only);gc()	

iters = 50 #iterations of classification
hidden_dim = 500 #hidden layer dimension
hidden_layer = 1 #number of hidden layers
for(i in new.ids){
  if(i > 0){
    other_ids = new.ids[new.ids > i]
    print(i)
    for(j in other_ids){
      if((i > 0 & j > i)){
        print(j)
        temp_tpm = tpm[metadata %in% c(i, j),]
        temp_meta = metadata[metadata %in% c(i, j)]
        iter_results = data.frame(Matrix(ncol = iters, nrow = dim(temp_tpm)))
        for(z in 1:iters){
          print(z)
          k = 4 ###fold split
          ###set seed here for reproducibility###
          set.seed(z*i*j)
          ##generate k partitions for training: in practice, we do this n times##
          k_partitions = list()
          selected_locs = c()
          sampling_set = 1:length(temp_meta)
          for(x in 1:k){
            if(x == 1){
              trainset = sample(sampling_set, length(temp_meta)/k, replace = F)
            } else {
              trainset = sample(sampling_set[-selected_locs], length(temp_meta)/k, replace = F)
            }
            selected_locs = c(selected_locs, trainset)
            if(x ==4){
              trainset = c(trainset, setdiff(sampling_set, selected_locs))
              selected_locs = c(selected_locs, setdiff(sampling_set, selected_locs))
            }
            
            k_partitions[[x]] = trainset
          }
          temp_meta = factor(temp_meta, levels = c(i,j))
          ###"one-hot" encoded output for training###
          y_train=sparseMatrix(i=1:length(temp_meta),j=as.numeric(temp_meta), x = 1)
          classification_comparison=c()
          classification_accuracy = c()
          all_old_classifications = data.frame(Matrix(ncol = length(unique(temp_meta)), nrow = 0))
          weights = list("0" = 1, "1" = table(temp_meta)[[1]]/table(temp_meta)[[2]])
          for(y in 1:k){
            ###input/output layers###
            model_1_input = layer_input(shape=ncol(temp_tpm)) ###Input layer: one node for each gene 
            model_1_output = model_1_input %>% 
              layer_dense(units=hidden_dim,activation="tanh") %>% ###activation function for input layer to hidden layer
              layer_dense(units=length(unique(temp_meta)),activation = "sigmoid")   ###activation function for hidden layer to output layer
            ###assemble and train model###
            model_1 = keras_model(model_1_input,model_1_output)  ###build the model
            model_1 %>% compile(  ###specify the training & error functions for the model
              loss = "mean_squared_error",
              optimizer = optimizer_rmsprop(),
              metrics = c("accuracy")
            )
            history <- model_1 %>% fit(  ###train the model - this will output a plot of accuracy as the model trains
              temp_tpm[-k_partitions[[y]],], y_train[-k_partitions[[y]],],
              epochs = 10, batch_size = 100,
              validation_split = 0.2,
              class_weight = weights
            )
            predvec=(model_1 %>% predict(temp_tpm[k_partitions[[y]],]))  ###run predictions for the test set, look only at the values not included in the training set
            row.names(predvec) = row.names(temp_tpm[k_partitions[[y]],])
            colnames(predvec) = levels(temp_meta)
            all_old_classifications = rbind(all_old_classifications, predvec)
          }
          predicted_cluster = colnames(all_old_classifications)[max.col(all_old_classifications,ties.method="first")]
          names(predicted_cluster) = rownames(all_old_classifications)
          predicted_cluster = predicted_cluster[order(names(predicted_cluster))]
          iter_results[,z] = predicted_cluster
          rownames(iter_results) = names(predicted_cluster)
        }
        colnames(iter_results) = 1:iters
        write.csv(iter_results, file = paste("constellation/constellationNN_", i, "_", j, ".csv", sep = ""))
      }
    }
  }
}