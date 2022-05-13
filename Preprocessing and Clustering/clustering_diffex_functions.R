###############################
#         jft2113             #
#         05/12/22            #
#         R4.0,Py3.7          #
###############################

######################################
#clustering/diffex pipeline functions#
######################################

######################################
##functions for reading in all data ##
###usage: pass character vector of file names (10X sample IDs), the directory in which data is stored, UMI cutoffs, and choose whether to filter RP genes and MT % ###
aggfilter_nonhash <- function(files, base_directory = "../raw_data/", umimin = 500, umimax = 10000, RP_filt = T, MT_filt = T){
  allmat=list()
  allgen=c()
  allcells=c()
  all_barcodes=c()
  #generate individual matrices per dataset
  for (i in 1:length(files)) {
    tempname=paste(base_directory, files[i], "/", sep ="")
    tempmat=Read10X_h5(file=paste(tempname, "filtered_feature_bc_matrix.h5", sep =""), use.names = T,unique.features = T)
    if(typeof(tempmat) == "list"){
      tempmat = tempmat[[1]]
    }
    barcodes = colnames(tempmat)
    tempcol = files[i]
    colnames(tempmat)=paste(tempcol,"-",1:ncol(tempmat),sep="")
    names(barcodes) = colnames(tempmat)
    ###mitochondrial filtering: choose the more $stringent$ of median or IQR filtering###
    if(MT_filt == T){
      mtpct=colSums(tempmat[grep("^MT-",rownames(tempmat)),])/colSums(tempmat)*100
      med_cutoff = median(mtpct) + 2*mad(mtpct)
      IQR_cutoff = as.numeric(quantile(mtpct)[4] + IQR(mtpct)*1.5)
      mitofilter = min(med_cutoff, IQR_cutoff)
      mitofilter = max(mitofilter, 10)
      print(paste(files[i], " mito: med-", med_cutoff, ", IQR-", IQR_cutoff, sep = ""))
      tempmat=tempmat[,which(mtpct<mitofilter)]  ###remove mito-heavy cells
      barcodes=barcodes[which(mtpct<mitofilter)]
    }
    ##filter genes unlikely to be biologically relevant out, such as ribosomal/mitochondrial genes, then filter on UMI count##
    if(RP_filt == T){
      tempmat=tempmat[grep("^MT-|^RP[0-9]|^BC[0-9]|^RPL|^RPS|^MTRNR|-PS",rownames(tempmat),invert=T),] ###remove all RP, MT, and pseudogenes
    }else if(MT_filt == T){
      tempmat=tempmat[grep("^MT-|^BC[0-9]|^MTRNR|-PS",rownames(tempmat),invert=T),] ###remove MT + PS
    } else {
      tempmat=tempmat[grep("^BC[0-9]|-PS",rownames(tempmat),invert=T),] ###remove MT + PS
    }
    umicount=colSums(tempmat)
    keepcols=which(umicount>umimin)
    tempmat=tempmat[,keepcols]
    barcodes = barcodes[keepcols]
    umicount=colSums(tempmat)
    keepcols=which(umicount<umimax)
    tempmat=tempmat[,keepcols]
    barcodes = barcodes[keepcols]
    allmat[[i]]=Matrix(tempmat,sparse=T)
    allgen=c(allgen,rownames(tempmat))
    allcells=c(allcells,colnames(tempmat))
    all_barcodes = c(all_barcodes, barcodes)
  }
  print("merging")
  #merge all matrices
  allgen=unique(allgen)
  alldat=Matrix(0,nrow=length(allgen),ncol=length(allcells),sparse=T)
  rownames(alldat)=allgen
  colnames(alldat)=allcells
  for (ii in 1:length(allmat)) {
    alldat[rownames(allmat[[ii]]),colnames(allmat[[ii]])]=allmat[[ii]]
  }
  return(list(alldat, all_barcodes))
}
###usage: pass character vector of file names (10X sample IDs), the directory in which data is stored, UMI cutoffs, HTO min, and choose whether to filter RP genes and MT % ###
aggfilter_hash <- function(files, base_directory = "../raw_data/", outdir = "../intermediate_data/demuxmix", umimin = 500, umimax = 10000, 
                           HTO_filter = 0, feature_sheet = "../intermediate_data/Table S1 - Overview of demographics, hashing strategy, and quality control, related to STAR Methods..xlsx", 
                           RP_filt = F, MT_filt = T){
  #retrieve information about hashtags found in each sample and check against CellRanger output#
  featuredata = read_excel(feature_sheet, sheet = 3)
  H <- new.env(hash = T)
  exists_hash <- Vectorize(exists, vectorize.args = "x")
  for(i in 1:length(featuredata$`Sequencing ID`)){
    if(exists_hash(featuredata$`Sequencing ID`[[i]], H)){
      H[[featuredata$`Sequencing ID`[[i]]]]= c(H[[featuredata$`Sequencing ID`[[i]]]], paste("Hashtag", colnames(featuredata)[4:10][!is.na(featuredata[i, 4:10])], sep = ""))
    } else{
      H[[featuredata$`Sequencing ID`[[i]]]]= paste("Hashtag", colnames(featuredata)[4:10][!is.na(featuredata[i, 4:10])], sep = "")
    }
  }
  #create seurat objects for each separate sample#
  seuratmat = list()
  hashtag_frame = data.frame(Matrix(nrow = 0, ncol = 27))
  for (i in 1:length(files)) {
    tempname=paste(base_directory, files[i], "/", sep ="")
    mmm=Read10X_h5(paste(tempname, "filtered_feature_bc_matrix.h5", sep = ""), use.names = T, unique.features = T)
    hashmat=mmm[[2]]
    cellmat=mmm[[1]]
    joint.bcs <- intersect(colnames(cellmat), colnames(hashmat))
    cellmat = cellmat[,joint.bcs]
    hashmat = hashmat[,joint.bcs]
    ###mito filtering: choose the more $stringent$ of median or IQR filtering###
    if(MT_filt == T){
      mtpct=colSums(cellmat[grep("^MT-",rownames(cellmat)),])/colSums(cellmat)*100
      med_cutoff = median(mtpct) + 2*mad(mtpct)
      IQR_cutoff = as.numeric(quantile(mtpct)[4] + IQR(mtpct)*1.5)
      mitofilter = min(med_cutoff, IQR_cutoff)
      mitofilter = max(mitofilter, 10)
      print(paste(files[i], " mito: med-", med_cutoff, ", IQR-", IQR_cutoff, sep = ""))
      cellmat=cellmat[,which(mtpct<mitofilter)]  ###remove mito-heavy cells
      hashmat=hashmat[,which(mtpct<mitofilter)]
    }
    ###filter genes unlikely to be biologically relevant out, such as ribosomal/mitochondrial genes, then filter on UMI count##
    if(RP_filt == T){
      cellmat=cellmat[grep("^MT-|^RP[0-9]|^BC[0-9]|^RPL|^RPS|^MTRNR|-PS",rownames(cellmat),invert=T),] ###to remove all RP, MT, and pseudogenes
    }else if (MT_filt == T){
      cellmat=cellmat[grep("^MT-|^BC[0-9]|^MTRNR|-PS",rownames(cellmat),invert=T),] ###to remove MT + PS
    } else{
      cellmat=cellmat[grep("^BC[0-9]|-PS",rownames(cellmat),invert=T),] ###to remove MT + PS
    }
    #filter including those cells with no hashtag reads#
    keepcells=which(colSums(cellmat) > umimin & colSums(cellmat) < umimax & colSums(hashmat) >= HTO_filter)
    #retrieve the hashtag data and add it to a seurat object#
    hashtag_frame = rbind(hashtag_frame, t(rowSums(hashmat[,keepcells])))
    rownames(hashtag_frame)[dim(hashtag_frame)[1]] = files[i]
    keephash = which(row.names(hashmat) %in% H[[files[i]]])
    dat_seurat = CreateSeuratObject(counts=cellmat[,keepcells])
    hashmat = hashmat[keephash,keepcells]
    dat_seurat[["HTO"]]=CreateAssayObject(counts=hashmat)
    tempcol=dat_seurat$orig.ident
    barcodes = names(tempcol)
    names(tempcol)=paste(files[i],"-",1:length(tempcol),sep="")
    dat_seurat = RenameCells(dat_seurat, new.names = names(tempcol))
    tempcol = rep(files[i], length(tempcol))
    dat_seurat$orig.ident = tempcol
    names(barcodes) = names(tempcol)
    dat_seurat = AddMetaData(dat_seurat, metadata = barcodes, col.name = "barcodes")
    
    #run hashtag deconvolution with demuxmix#
    pdf(paste(outdir, dat_seurat$orig.ident[[1]], "_", umimin, "_demuxMM.pdf", sep = ""), onefile = T, paper = "a4r", width = 9, height = 12)
    for(j in dat_seurat@assays$HTO@counts@Dimnames[[1]]){
      dat_seurat <- demuxMM(object = dat_seurat, tag = j, errorProb = 0.000001, quantile = 0.99)
      plotMM(dat_seurat, tag = j)
    }
    dev.off()
    
    ###assign labels via demuxmix###
    dat_seurat = assignLabelsX(object = dat_seurat, tags = dat_seurat@assays$HTO@counts@Dimnames[[1]])
    #run HTODemux and MULTI demux for alternate use#
    dat_seurat = NormalizeData(dat_seurat, assay = "HTO", normalization.method = "CLR", margin = 1)
    dat_seurat = HTODemux(dat_seurat, assay = "HTO", kfunc = "clara", nsamples = 100, positive.quantile = 0.999)
    dat_seurat = MULTIseqDemux(dat_seurat, assay = "HTO", autoThresh = T, maxiter = 10, qrange = seq(from = 0.05, to = 0.95, by = 0.05))
    seuratmat[[i]]=dat_seurat
  }
  return(list(seuratmat, hashtag_frame))
}

############################################################
##prototype of demuxmix used for deconvolution of hashtags##
#       see https://github.com/cu-ctcn/demuxmix            #
mmNegBinom <- function(c, tol=0.001, maxIter=100, rm.outlier=TRUE) {
  # Initial estimation using k-means. Cluster with smaller mean is id=1
  pi <- c(0.5, 0.5)
  ind <- kmeans(log(c+1), centers = 2)$cluster == 1
  if (median(c[ind]) > median(c[!ind])) {
    ind <- !ind
  }
  if (rm.outlier) {
    outlierTh <- mean(c[!ind]) + 3 * sd(c[!ind])
    outlierInd <- c > outlierTh
    if (sum(outlierInd) > 0) {
      message(paste("Removed", sum(outlierInd), "outlier cells before fitting mixture model. Threshold:", outlierTh, " reads."))
    }
    c <- c[!outlierInd]
    ind <- ind[!outlierInd]
  }
  mu <- c(mean(c[ind]), mean(c[!ind]))
  var <- c(var(c[ind]), var(c[!ind]))
  size <- c(mu[1]^2/ (var[1] - mu[1]), mu[2]^2/ (var[2] - mu[2]))
  ## runs into problems when mu >var, so set size to Inf if negative to make this a poisson in that edge case ##
  if(mu[1] == var[1]){
    if(mu[1] == 0){
      size[is.na(size)] = 0.0000001
      mu[1] = 0.0001
      var[1] = 0.0002
    }
    message("identical (0) mu and var @ input to EM")
  }
  if(sum(size < 0) > 0){
    message("negative size param @ input to EM")
  }
  size[size<0] = Inf
  
  iter <- 1
  logLik <- 1
  deltaLogLik <- Inf
  logData <- data.frame()
  while (iter < maxIter & abs(deltaLogLik/logLik) > tol) {
    
    # E-step
    f <- matrix(nrow=length(c), ncol=2)
    f[, 1] <- dnbinom(c, size=abs(size[1]), mu=mu[1])
    f[, 2] <- dnbinom(c, size=size[2], mu=mu[2])
    pi.f <- t(pi * t(f))
    class <- pi.f / apply(pi.f, 1, sum)
    pi <- apply(class, 2, sum) / nrow(class)
    
    # M-step
    w <- apply(class, 2, function(x) {x / sum(x)})
    mu[1] <- sum(w[, 1] * c)
    var[1] <- sum(w[, 1] * (c - mu[1])^2)
    size[1] <- mu[1]^2/ (var[1] - mu[1])
    mu[2] <- sum(w[, 2] * c)
    var[2] <- sum(w[, 2] * (c - mu[2])^2)
    size[2] <- mu[2]^2/ (var[2] - mu[2])
    
    iter <- iter + 1
    oldLogLik <- logLik
    logLik <- sum(log(apply(pi.f, 1, sum)))
    deltaLogLik <- logLik - oldLogLik
    
    logData <- rbind(logData,
                     data.frame(Iteration=iter - 1,
                                LogLikelihood=logLik,
                                DeltaLogLikelihood=deltaLogLik,
                                mu1=mu[1],
                                mu2=mu[2],
                                size1=size[1],
                                size2=size[2]))
  }
  
  results <- data.frame(c=c,
                        prob1=class[, 1],
                        prob2=class[, 2],
                        classification=apply(class, 1, which.max))
  ## runs into problems when mu >var, so set size to Inf if negative to make this a poisson in that edge case ##
  if(sum(size < 0) > 0){
    message("negative size param @ exit from EM")
  }
  size[size<0] = Inf
  
  return(list(results=results,
              logData=logData,
              converged=iter<maxIter,
              mu=mu,
              size=size,
              var=var,
              pi=pi))
}
# general setting for errorProb 1e-06 #
# general setting for quantile 0.99 #
demuxMM <- function(object, tag="Hashtag7", errorProb, quantile, seed=168) {
  if ((missing(errorProb) & !missing(quantile)) | (!missing(errorProb) & missing(quantile))) {
    stop("Either errorProb or quantile have to be defined.")
  }
  set.seed(seed=seed)
  counts <- GetAssayData(object=object, assay="HTO", slot="counts")
  c <- counts[tag, ]
  r <- mmNegBinom(c)
  
  if (!missing(quantile)) {
    csfQ <- quantile
    pbmcQ <- 1 - quantile
    csfT <- qnbinom(csfQ, mu=r$mu[1], size=r$size[1])
    pbmcT <- qnbinom(pbmcQ, mu=r$mu[2], size=r$size[2])
    csfErrorP <- pnbinom(csfT, mu=r$mu[2], size=r$size[2], lower.tail=TRUE)
    pbmcErrorP <- pnbinom(pbmcT, mu=r$mu[1], size=r$size[1], lower.tail=FALSE)
  }
  if (!missing(errorProb)) {
    csfErrorP <- errorProb
    pbmcErrorP <- errorProb
    csfT <- qnbinom(errorProb, mu=r$mu[2], size=r$size[2], lower.tail=TRUE)
    pbmcT <- qnbinom(errorProb, mu=r$mu[1], size=r$size[1], lower.tail=FALSE)
    csfQ <- pnbinom(csfT, mu=r$mu[1], size=r$size[1], lower.tail=TRUE)
    pbmcQ <- pnbinom(pbmcT, mu=r$mu[2], size=r$size[2], lower.tail=FALSE)
  }
  
  if (csfT >= pbmcT) {
    stop("Distribution overlap too much. Choose larger error probabilities.")
  }
  
  classification <- data.frame(tag=rep("uncertain", ncol(object)), outlier=FALSE, stringsAsFactors = FALSE)
  colnames(classification) <- c(tag, paste(tag, ".", "outlier", sep=""))
  rownames(classification) <- colnames(object)
  classification[c <= csfT, 1] <- "negative"
  classification[c >= pbmcT, 1] <- "positive"
  classification[, 2] <- ! rownames(classification) %in% rownames(r$results)
  object <- AddMetaData(object=object, metadata=classification)
  
  info <- data.frame(tag=tag, component=c("1", "2"), mu=r$mu, var=r$var, size=r$size, pi=r$pi,
                     threshold=c(csfT, pbmcT), quant=c(csfQ, pbmcQ), errorP=c(csfErrorP, pbmcErrorP),
                     stringsAsFactors=FALSE)
  if ("HtoDemux" %in% names(object@misc)) {
    object@misc$HtoDemux <- rbind(object@misc$HtoDemux[object@misc$HtoDemux$tag != tag, ], info)
  } else {
    object@misc$HtoDemux <- info
  }
  
  return(object)
}
#assign labels to object based on demuxmix results
assignLabelsX <- function(object, tags=c()) {
  
  stopifnot(length(tags) > 1)
  classification <- object@meta.data[, tags[1], drop=FALSE]
  colnames(classification) <- "Label"
  classification$Label <- rep(NA, nrow(classification))
  
  tag_X = list()
  ###make your tag lists###
  for(i in 1:length(tags)){
    tag_X[[i]] = object@meta.data[, tags[i], drop=TRUE]
  }
  ###assign all positives first###
  for(i in 1:length(tags)){
    classification$Label[tag_X[[i]] == "positive"] <- tags[i]
  }
  ###assign doublets next###
  for(i in 1:length(tag_X)){
    locs_X = which(tag_X[[i]] == "positive")
    other_locs = c()
    for(j in tag_X[-i]){
      other_locs = union(other_locs, which(j == "positive"))
    }
    classification$Label[intersect(locs_X, other_locs)] = "doublet"
  }
  ###now assign uncertain/negative###
  locs_uncertain = which(tag_X[[1]] %in% c("negative", "uncertain"))
  locs_neg = which(tag_X[[1]] == "negative")
  for(j in tag_X[-1]){
    locs_uncertain = intersect(locs_uncertain, which(j %in% c("negative", "uncertain")))
    locs_neg = intersect(locs_neg, which(j == "negative"))
  }
  classification$Label[locs_uncertain] = "uncertain"
  classification$Label[locs_neg] = "negative"
  
  ###assign labels###
  object <- AddMetaData(object=object, metadata=classification)
  return(object)
}
#for plotting QC graphs from demuxmix#
plotMM <- function(object, tag="Hashtag7", bins=100, xlim, ylim) {
  counts <- GetAssayData(object=object, assay="HTO", slot="counts")
  c <- counts[tag, !object@meta.data[, paste(tag, ".outlier", sep="")]]
  
  if (missing(xlim)) {
    xlim <- NULL
  } else {
    bins=bins * (max(c) / (xlim[2]-xlim[1]))
  }
  if (missing(ylim)) {
    ylim <- NULL
  }
  
  xRange <- seq(max(min(c), xlim[1]), min(max(c), xlim[2]))
  d <- object@misc$HtoDemux
  dens1 <- data.frame(x=xRange,
                      y=d[d$tag == tag & d$component == "1", "pi"] * dnbinom(xRange, mu=d[d$tag == tag & d$component == "1", "mu"], size=d[d$tag == tag & d$component == "1", "size"]),
                      Component="1", stringsAsFactors=FALSE)
  dens2 <- data.frame(x=xRange,
                      y=d[d$tag == tag & d$component == "2", "pi"] * dnbinom(xRange, mu=d[d$tag == tag & d$component == "2", "mu"], size=d[d$tag == tag & d$component == "2", "size"]),
                      Component="2", stringsAsFactors=FALSE)
  densM <- data.frame(x=xRange,
                      y=dens1$y + dens2$y,
                      Component="M", stringsAsFactors=FALSE)
  dens <- rbind(dens1, dens2, densM)
  
  if (is.null(xlim) & is.null(ylim)) {
    hist(c, breaks=bins, probability=TRUE, main="", xlab=tag)
  } else if (!is.null(xlim) & is.null(ylim)) {
    hist(c, breaks=bins, probability=TRUE, xlim=xlim, main="", xlab=tag)
  } else if (is.null(xlim) & !is.null(ylim)) {
    hist(c, breaks=bins, probability=TRUE, ylim=ylim, main="", xlab=tag)
  } else {
    hist(c, breaks=bins, probability=TRUE, ylim=ylim, xlim=xlim, main="", xlab=tag)
  }
  lines(x=dens$x[dens$Component == "M"], y=dens$y[dens$Component == "M"], lwd=2)
  lines(x=dens$x[dens$Component == "1"], y=dens$y[dens$Component == "1"], col="dodgerblue", lwd=2)
  lines(x=dens$x[dens$Component == "2"], y=dens$y[dens$Component == "2"], col="red", lwd=2)
  abline(v=object@misc$HtoDemux$threshold[object@misc$HtoDemux$tag == tag &  object@misc$HtoDemux$component == "1"], col="dodgerblue", lwd=2, lty=2)
  abline(v=object@misc$HtoDemux$threshold[object@misc$HtoDemux$tag == tag &  object@misc$HtoDemux$component == "2"], col="red", lwd=2, lty=2)
}
#for integrating assignments from distinct deconvolution approaches#
###select final classifications by manually choosing sample/hashtag pairings where some should be drawn from HTOdemux/MULTI and others from demuxMM### 
integrate_seurat_assignments <- function(sample = "PM077", seurat, demuxMM_tags = c("Hashtag7", "Hashtag8", "Hashtag9"), alt_demux = "HTO_demux", other_tags = c("Hashtag6", "Hashtag10"))
{
  sample_locs = seurat$orig.ident == sample
  hashtag_final = data.frame(seurat$orig.ident[sample_locs], seurat$Label[sample_locs], seurat$Label[sample_locs],stringsAsFactors = F)
  if(alt_demux == "HTO_demux"){
    HTO_outputs = data.frame(seurat$orig.ident[sample_locs], seurat$HTO_classification[sample_locs])
  } 
  if(alt_demux == "MULTI"){
    HTO_outputs = data.frame(seurat$orig.ident[sample_locs], seurat$MULTI_classification[sample_locs])
  }
  colnames(HTO_outputs) = c("orig.ident", "HTO_classification")
  ###aggregate the demuxMM assignments for easier reference
  for(i in demuxMM_tags){
    hashtag_final = cbind(hashtag_final, seurat@meta.data[sample_locs, i])
  }
  colnames(hashtag_final) = c("orig.ident", "label", "final_id", demuxMM_tags)
  
  ###identify singlets/doublets for the Hashtags where we plan to use the other approach###
  seurat_doublet = list()
  all_doublet = c()
  for(i in 1:length(other_tags)){
    locs = grep(other_tags[[i]], HTO_outputs$HTO_classification)
    for(j in demuxMM_tags){
      demux_pos = intersect(locs, which(hashtag_final[,j] %in% "positive"))
      all_doublet = union(all_doublet, demux_pos)
    }
    hashtag_final$final_id[locs] = HTO_outputs$HTO_classification[locs]
    seurat_doublet[[i]] = intersect(locs, grep("_", HTO_outputs$HTO_classification))
  }
  hashtag_final$final_id[all_doublet] = "doublet"
  
  ###sort the list of seurat doublets to identify only those that come from the other approach only, since we have already assigned other demuxmix doublets
  true_seurat_doublets = c()
  for(i in 1:length(seurat_doublet)){
    for(j in seurat_doublet[-i]){
      true_seurat_doublets = c(true_seurat_doublets, intersect(seurat_doublet[[i]], j))
    }
  }
  ##assign all seurat-seurat barcode doublets
  hashtag_final$final_id[true_seurat_doublets] = "doublet"
  ##finally, assign the final labels here based on the combo of demuxmix and the other calling method
  for(i in grep("_", hashtag_final$final_id)){
    if(hashtag_final$label[i] %in% c(other_tags, "uncertain", "negative", "doublet")){
      hashtag_final$final_id[i] = strsplit(hashtag_final$final_id[i], "_")[[1]][strsplit(hashtag_final$final_id[i], "_")[[1]] %in% other_tags][1]
    } else{
      hashtag_final$final_id[i] = hashtag_final$label[i]
    }
  } 
  seurat$Label[sample_locs] = hashtag_final$final_id
  return(seurat)
}

###########################################################
#functions for building  seurat object and adding metadata#
#build the seurat object, add metadata, and fix hashes#
makeseurat_metadata <- function(alldat, allhash,sampdata, featuredata){
  hash = merge(allhash[[1]], allhash[-1])
  ##Integrate different hashtags depending on the data in question##
  hash = SetIdent(object = hash, value = hash$Label)
  hash = integrate_seurat_assignments(sample = "PM064", seurat = hash, other_tags = c("Hashtag7", "Hashtag8"), alt_demux = "MULTI",
                                      demuxMM_tags = c("Hashtag9", "Hashtag10"))
  hash = integrate_seurat_assignments(sample = "PM066", seurat = hash, other_tags = c("Hashtag9"), alt_demux = "MULTI",
                                      demuxMM_tags = c("Hashtag7", "Hashtag8"))
  hash = integrate_seurat_assignments(sample = "PM078", seurat = hash, other_tags = c("Hashtag9", "Hashtag10"), alt_demux = "HTO_demux",
                                      demuxMM_tags = c("Hashtag8"))
  hash = integrate_seurat_assignments(sample = "PM080", seurat = hash, other_tags = c("Hashtag6", "Hashtag9"), alt_demux = "MULTI",
                                      demuxMM_tags = c("Hashtag10"))
  hash = integrate_seurat_assignments(sample = "PM082", seurat = hash, other_tags = c("Hashtag8", "Hashtag9"), alt_demux = "MULTI",
                                      demuxMM_tags = c("Hashtag10"))
  hash = integrate_seurat_assignments(sample = "PM083", seurat = hash, other_tags = c("Hashtag8", "Hashtag9"), alt_demux = "MULTI",
                                      demuxMM_tags = c("Hashtag10"))
  hash = integrate_seurat_assignments(sample = "PM084", seurat = hash, other_tags = c("Hashtag8", "Hashtag9"), alt_demux = "MULTI",
                                      demuxMM_tags = c("Hashtag10"))
  hash = SetIdent(object = hash, value = hash$Label)
  hash.subset <- subset(hash, idents = c("negative", "uncertain", "doublet"), invert = TRUE)
  ####create Seurat object and add metadata####
  uglia.seurat = CreateSeuratObject(counts = alldat, min.cells = 3, names.delim = "-", names.field = 1, project = "uglia-complete")
  if(!is.null(hash.subset)){
    uglia.seurat = merge(uglia.seurat, hash.subset)
  }
  ##add metadata##
  meta_params = c("Diagnosis (neurology)", "Brain region", "10x Chemistry Version", "Sex", "Age", "Donor ID")
  names(meta_params) = c("diagnosis", "region", "tech", "sex", "age", "donor_ID")
  for(i in 1:length(meta_params)){
    meta_col = rep(NA, length(Idents(uglia.seurat)))
    for(sample in unique(uglia.seurat$orig.ident)){
      samploc = grep(sample, sampdata$"Sequencing ID")
      col = as.character(sampdata[samploc, meta_params[i]])
      finalcol = grep(sample, uglia.seurat$orig.ident)
      meta_col[finalcol] = col
    }
    uglia.seurat= AddMetaData(object = uglia.seurat, metadata = meta_col, col.name = names(meta_params)[i])
  }

  ###correct the incorrect metadata
  uglia.seurat$donor_ID[uglia.seurat$orig.ident == "PM077"] = "ALS/FTD2"
  uglia.seurat$diagnosis[uglia.seurat$orig.ident == "PM077"] = "ALS/FTD"
  uglia.seurat$tech[uglia.seurat$orig.ident == "PM077"] = "v3"
  uglia.seurat$region[uglia.seurat$orig.ident == "PM077"] = "BA9+BA4+SC"
  uglia.seurat$sex[uglia.seurat$orig.ident == "PM077"] = "F"
  uglia.seurat$age[uglia.seurat$orig.ident == "PM077"] = 60
  ###
  uglia.seurat$donor_ID[uglia.seurat$orig.ident == "PM086"] = "LOAD35"
  uglia.seurat$diagnosis[uglia.seurat$orig.ident == "PM086"] = "LOAD"
  uglia.seurat$tech[uglia.seurat$orig.ident == "PM086"] = "v3"
  uglia.seurat$region[uglia.seurat$orig.ident == "PM086"] = "BA9+BA20+H"
  uglia.seurat$sex[uglia.seurat$orig.ident == "PM086"] = "F"
  uglia.seurat$age[uglia.seurat$orig.ident == "PM086"] = 78
  uglia.seurat$Label[is.na(uglia.seurat$Label)] = "None"
  
  #fix the two mixed samples
  uglia.seurat$diagnosis[uglia.seurat$orig.ident == "PM077" & (uglia.seurat$Label == "Hashtag6" | uglia.seurat$Label == "Hashtag7")] = "DLBD-PD"
  uglia.seurat$region[uglia.seurat$orig.ident == "PM077" & (uglia.seurat$Label == "Hashtag6" | uglia.seurat$Label == "Hashtag7")] = "BA9+SN"
  uglia.seurat$sex[uglia.seurat$orig.ident == "PM077" & (uglia.seurat$Label == "Hashtag6" | uglia.seurat$Label == "Hashtag7")] = "F"
  uglia.seurat$age[uglia.seurat$orig.ident == "PM077" & (uglia.seurat$Label == "Hashtag6" | uglia.seurat$Label == "Hashtag7")] = 82
  uglia.seurat$donor_ID[uglia.seurat$orig.ident == "PM077" & (uglia.seurat$Label == "Hashtag6" | uglia.seurat$Label == "Hashtag7")] = "DLBD-PD5"
  uglia.seurat$diagnosis[uglia.seurat$orig.ident == "PM086" &
                           (uglia.seurat$Label == "Hashtag5" | uglia.seurat$Label == "Hashtag6" | uglia.seurat$Label == "Hashtag7")] = "HD"
  uglia.seurat$region[uglia.seurat$orig.ident == "PM086" &
                        (uglia.seurat$Label == "Hashtag5" | uglia.seurat$Label == "Hashtag6" | uglia.seurat$Label == "Hashtag7")] = "BA9+BA20+H"
  uglia.seurat$donor_ID[uglia.seurat$orig.ident == "PM086" &
                            (uglia.seurat$Label == "Hashtag5" | uglia.seurat$Label == "Hashtag6" | uglia.seurat$Label == "Hashtag7")] = "HD1"
  ###reassign regions for hashtag deconvolution###
  uglia.seurat = SetIdent(object = uglia.seurat, value = uglia.seurat$Label)
  featuredata = featuredata[featuredata$`Sequencing ID` %in% unique(uglia.seurat$orig.ident),]
  for(i in 1:length(featuredata$`Sequencing ID`)){
    hashes = unique(uglia.seurat$Label[uglia.seurat$orig.ident == featuredata$`Sequencing ID`[i]])
    value = substr(hashes, 8, 9)
    for(j in 1:length(hashes)){
      uglia.seurat$region[uglia.seurat$orig.ident == featuredata$`Sequencing ID`[i] &
                            uglia.seurat$Label == hashes[j]] = as.character(featuredata[i,value[j]])
    }
  }
  donor_region = paste(uglia.seurat$donor_ID, uglia.seurat$region, sep = "_")
  uglia.seurat = AddMetaData(object = uglia.seurat, metadata = donor_region, col.name = "donor_region")
  rm(alldat, allhash, hash, hash.subset); gc()
  return(uglia.seurat)
}

###SCT + mNN
sct_mnn <- function(seurat, regress = c("nCount_RNA", "tech"), num_var_genes = 3000, resolution = c(0.25, 0.5, 1.00), 
                    pc_dimensions = 20, topnum = 5, auto_merge = F, conserve.memory = F){
  samp.combined <- SCTransform(seurat, do.correct.umi = T, vars.to.regress = regress, do.scale = F, do.center = T, conserve.memory = conserve.memory,
                               return.only.var.genes = T, variable.features.n = num_var_genes, seed.use = 1448145, verbose = T) ### could try this with clipping at residual variance cutoff
  samp.combined <- RunFastMNN(object.list = SplitObject(samp.combined, split.by = "orig.ident"), assay = "SCT", features = num_var_genes, auto.merge = auto_merge)
  samp.combined <- RunUMAP(samp.combined, reduction = "mnn", dims = 1:pc_dimensions, seed.use = 42L)
  samp.combined <- FindNeighbors(samp.combined, reduction = "mnn", dims = 1:pc_dimensions)
  samp.combined <- FindClusters(samp.combined, resolution = resolution)
  return(samp.combined)
}

###################################################################################################
###generate pairwise gene lists for each cluster vs. other clusters using MAST and Seurat diffex###
#perform all pairwise diffex testing#
all_pairwise_diffex <- function(seurat, outdir = "2020_Analysis_Outs/Pairwise_gene_lists/ALS_+RP500/",
                                assay = "SCT", slot = "data", test = "wilcox", fc = 0.25)
{
  #choose all clusters in active identity
  clusters = sort(as.numeric(as.character(unique(seurat@active.ident))), decreasing = F)
  alldegenes = data.frame()
  ###take the pairwise genes for each cluster vs. all other clusters based on SCT normalized data
  for(i in 1:length(clusters)){
    other_clust = clusters[-i]
    for(j in other_clust){
      pwise_markers = FindMarkers(seurat, ident.1 = clusters[i], ident.2 = j, test.use = test,
                                  min.pct = 0.05, only.pos = F, assay = assay, slot = slot, logfc.threshold = fc)    
      pwise_markers["cluster"] = clusters[i]
      pwise_markers["vs"] = j
      pwise_markers$gene = row.names(pwise_markers)
      write.csv(pwise_markers, paste(outdir, clusters[i], "vs", j, "markers.csv", sep = ""))
      alldegenes = rbind(alldegenes, pwise_markers[pwise_markers$p_val_adj < 0.05,]) ### pre-filter for genes <0.05
    }
  }
  write.csv(alldegenes, paste(outdir, "alldegenes.csv", sep = ""))
}
#retain final lists of pairwise genes#
#1.pre-filter DEG list to include only genes with adjusted p-values<0.05###
#2.find all genes upregulated in one population versus any other populations with fold-change>2 and fraction of expression>0.05, but not downregulated in said population versus anything else, regardless of fraction or fold-change
#usage: pass a threshold, fraction of expression, a choice of populations to compare (i.e. 0-4 or 0, 2, and 8), and choose either upreg or downreg
pairwise_genes <- function(DEG_file = "2020_Analysis_Outs/Pairwise_gene_lists/ALS_+RP500/alldegenes.csv", 
                           logFC_threshold=log(2), fraction_expression=0.05, allpops=0:8, upreg = T){
  alldegenes = read.csv(DEG_file)
  fullmat=c()   ###output matrix
  if(upreg == F & logFC_threshold > 0){
    logFC_threshold = -logFC_threshold
  }
  for (popval in allpops) {
    if(upreg == T){
      subgenes=setdiff(alldegenes$gene[alldegenes$cluster==popval & alldegenes$vs %in% allpops & alldegenes$avg_logFC>logFC_threshold & alldegenes$pct.1>fraction_expression],
                       alldegenes$gene[alldegenes$cluster==popval & alldegenes$vs %in% allpops & alldegenes$avg_logFC<0])
    } else {
      subgenes=setdiff(alldegenes$gene[alldegenes$cluster==popval & alldegenes$vs %in% allpops & alldegenes$avg_logFC<logFC_threshold & alldegenes$pct.1>fraction_expression],
                       alldegenes$gene[alldegenes$cluster==popval & alldegenes$vs %in% allpops & alldegenes$avg_logFC>0])
    }
    newmat=c()
    for (jj in subgenes) {
      if(upreg == T){
        submat=alldegenes[alldegenes$gene==jj & alldegenes$avg_logFC> logFC_threshold & alldegenes$cluster==popval & alldegenes$vs %in% allpops,]
      } else {
        submat=alldegenes[alldegenes$gene==jj & alldegenes$avg_logFC< logFC_threshold & alldegenes$cluster==popval & alldegenes$vs %in% allpops,]
      }
      numpops=length(submat$vs)
      listpop=paste(submat$vs,collapse=";")
      submat$p_val = submat$p_val + 1e-200 ###slight perturbation to make values non-0
      if(upreg == T){
        tempmat=data.frame(gene=jj,up_type=popval,down_types=listpop,num_down_types=numpops,sum_log10pvals=sum(log10(submat$p_val)),mean_log10pval=mean(log10(submat$p_val)),min_log10pval=min(log10(submat$p_val)),sum_logFC=sum(submat$avg_logFC),mean_logFC=mean(submat$avg_logFC),min_logFC=min(submat$avg_logFC),stringsAsFactors=F)
      } else{
        tempmat=data.frame(gene=jj,down_type=popval,up_types=listpop,num_up_types=numpops,sum_log10pvals=sum(log10(submat$p_val)),mean_log10pval=mean(log10(submat$p_val)),min_log10pval=min(log10(submat$p_val)),sum_logFC=sum(submat$avg_logFC),mean_logFC=mean(submat$avg_logFC),min_logFC=min(submat$avg_logFC),stringsAsFactors=F)
      }
      newmat=rbind(newmat,tempmat)
      if(upreg == T){
        newmat=newmat[order(-newmat$num_down_types),]
      }
      else{
        newmat=newmat[order(-newmat$num_up_types),]
      }
    }
    fullmat=rbind(fullmat,newmat)
  }
  return(fullmat) 
}

###################################
#label transfer pipeline functions#
#upstream mNN integration#
mNN_pipeline <- function(seurat, merge_order, DEG_genes, num_var_genes = 4500, pc_dimensions = 40, resolution = c(0.35, 0.5, 0.75)){
  VariableFeatures(seurat) = unique(DEG_genes)
  uglia_list = SplitObject(seurat, split.by = "orig.ident")
  train_merge <- RunFastMNN(object.list = uglia_list, features = num_var_genes, auto.merge = F, merge.order = match(merge_order, names(uglia_list)))
  train_merge <- RunUMAP(train_merge, reduction = "mnn", dims = 1:pc_dimensions)
  train_merge <- FindNeighbors(train_merge, reduction = "mnn", dims = 1:pc_dimensions)
  train_merge <- FindClusters(train_merge, resolution = resolution)
  return(train_merge)
}
#load in a query dataset#
load_dataset <- function(dataset){
  if(dataset == "xeno"){
    ###load in the Blurton-Jones iPSC data-Hassellmann###
    data_stub = "../other_data/Hasselmann_2019_in_vivo/"
    mtx_files = list.files(data_stub, pattern = "\\.mtx$")
    barcode_files = list.files(data_stub, pattern = "\\_barcodes")
    gene_files = list.files(data_stub, pattern = "\\_genes")
    ###aggregate and process to add back in gene/sample names
    final_mat = Matrix()
    for(i in 1:length(mtx_files)){
      tmp_mat = readMM(file = paste(data_stub, mtx_files[[i]], sep= ""))
      tmp_BC = read.table(file = paste(data_stub, barcode_files[[i]], sep = ""), sep = "\t")
      tmp_gene = read.table(file = paste(data_stub, gene_files[[i]], sep = ""), sep = "\t")
      tmp_gene$gene = sapply(strsplit(tmp_gene$V2, "_"),"[", 2)
      ###add gene and BC to matrix###
      row.names(tmp_mat) = tmp_gene$gene; 
      colnames(tmp_mat) = paste(strsplit(mtx_files[[i]],"_")[[1]][2], 1:dim(tmp_mat)[[2]], sep = "-")
      if(i == 1){
        final_mat = tmp_mat
      }  else{
        tmp_mat = tmp_mat[rownames(final_mat,),]
        ###make sure all the rows/genes are the same###
        if(any(!(row.names(tmp_mat) == row.names(final_mat)))){
          stop("misaligned/missing genes")
        }
        final_mat = cbind(final_mat, tmp_mat)
      }
    }
    seurat = CreateSeuratObject(counts = final_mat, min.cells = 3, names.delim = "-", names.field = 1, project = "uglia-complete")
    seurat <- NormalizeData(seurat)
  }else{
    print("unknown_input")
    stop()
  }
  return(seurat)
}