#' evalaute performance
synthetic_performance <- function (test_spots_metadata_mtrx, spot_composition_mtrx) 
{
    if (!is.matrix(test_spots_metadata_mtrx)) 
        stop("ERROR: test_spots_metadata_mtrx must be a matrix object!")
    if (!is.matrix(spot_composition_mtrx)) 
        stop("ERROR: syn_spots_ls must be the list obtained from the function syn_spot_comb_topic_fun().")
    colnames(spot_composition_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", 
        ".", x = colnames(spot_composition_mtrx), perl = TRUE)
    colnames(test_spots_metadata_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", 
        ".", x = colnames(test_spots_metadata_mtrx), perl = TRUE)
    suppressMessages(require(philentropy))
    true_jsd_mtrx <- matrix(nrow = nrow(test_spots_metadata_mtrx), 
        ncol = 1)
    tp <- 0
    tn <- 0
    fp <- 0
    fn <- 0
    for (i in seq_len(nrow(test_spots_metadata_mtrx))) {
        x <- rbind(test_spots_metadata_mtrx[i, ], spot_composition_mtrx[i, 
            ])
        if (sum(spot_composition_mtrx[i, ]) > 0) {
            true_jsd_mtrx[i, 1] <- suppressMessages(JSD(x = x, 
                unit = "log2", est.prob = "empirical"))
        }
        else {
            true_jsd_mtrx[i, 1] <- 1
        }
        for (index in colnames(test_spots_metadata_mtrx)) {
            if (x[1, index] > 0 & x[2, index] > 0) {
                tp <- tp + 1
            }
            else if (x[1, index] == 0 & x[2, index] == 0) {
                tn <- tn + 1
            }
            else if (x[1, index] > 0 & x[2, index] == 0) {
                fn <- fn + 1
            }
            else if (x[1, index] == 0 & x[2, index] > 0) {
                fp <- fp + 1
            }
        }
        rm(index)
    }
    rm(i)
    accuracy <- round((tp + tn)/(tp + tn + fp + fn), 2)
    sensitivity <- round(tp/(tp + fn), 2)
    specificity <- round(tn/(tn + fp), 2)
    precision <- round(tp/(tp + fp), 2)
    recall <- round(tp/(tp + fn), 2)
    F1 <- round(2 * ((precision * recall)/(precision + recall)), 
        2)
    quants_jsd <- round(quantile(matrixStats::rowMins(true_jsd_mtrx, 
        na.rm = TRUE), c(0.25, 0.5, 0.75)), 5)
    cat(sprintf("The following summary statistics are obtained:\n              Accuracy: %s,\n              Sensitivity: %s,\n              Specificity: %s,\n              precision: %s,\n              recall: %s,\n              F1 score: %s,\n              JSD quantiles: %s[%s-%s]", 
        accuracy, sensitivity, specificity, precision, recall, 
        F1, quants_jsd[[2]], quants_jsd[[1]], quants_jsd[[3]]), 
        sep = "\n")
    cat("raw statistics are returned in the list - TP, TN, FP, FN, JSD quantiles", 
        sep = "\n")
    return(list(JSD = true_jsd_mtrx))
}

#' normalize function
normalize_data <- function(count.list){
    norm.list <- vector('list')
    var.features <- vector('list')
    for ( i in 1:length(count.list)){
        norm.list[[i]] <- as.matrix(Seurat:::NormalizeData.default(count.list[[i]]))
        hvf.info <- Seurat:::FindVariableFeatures.default(count.list[[i]],selection.method='vst')
        hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 0), ]
        hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
        var.features[[i]] <- head(rownames(hvf.info), n = 2000)
    }
    sel.features <- selectIntegrationFeature(count.list,var.features)
    return (list(norm.list,sel.features))}

#' scaling function 
scale_data <- function(count.list,norm.list,hvg.features){
    scale.list <- lapply(norm.list,function(mat){
        Seurat:::ScaleData.default(object = mat, features = hvg.features)})
    scale.list <- lapply(1:length(count.list),function(i){
        return (scale.list[[i]][na.omit(match(rownames(count.list[[i]]),rownames(scale.list[[i]]))),])})
    return (scale.list)}

#' select HVG genes
selectIntegrationFeature <- function(count.list,var.features,nfeatures = 2000){
    var.features1 <- unname(unlist(var.features))
    var.features2 <- sort(table(var.features1), decreasing = TRUE)
    for (i in 1:length(count.list)) {
        var.features3 <- var.features2[names(var.features2) %in% rownames(count.list[[i]])]}    
    tie.val <- var.features3[min(nfeatures, length(var.features3))]
    features <- names(var.features3[which(var.features3 > tie.val)])
    if (length(features) > 0) {
        feature.ranks <- sapply(features, function(x) {
            ranks <- sapply(var.features, function(y) {
                if (x %in% y) {
                    return(which(x == y))
                }
                return(NULL)
            })
            median(unlist(ranks))
        })
        features <- names(sort(feature.ranks))
    }
    features.tie <- var.features3[which(var.features3 == tie.val)]
    tie.ranks <- sapply(names(features.tie), function(x) {
        ranks <- sapply(var.features, function(y) {
            if (x %in% y) {return(which(x == y))}
            return(NULL)
        })
        median(unlist(ranks))
    })
    features <- c(features, names(head(sort(tie.ranks), nfeatures - length(features))))
    return(features)
}

#' select variable genes
select_feature <- function(data,label,nf=2000){
    M <- nrow(data); new.label <- label[,1]
    pv1 <- sapply(1:M, function(i){
        mydataframe <- data.frame(y=as.numeric(data[i,]), ig=new.label)
        fit <- aov(y ~ ig, data=mydataframe)
        summary(fit)[[1]][["Pr(>F)"]][1]})
    names(pv1) <- rownames(data)
    pv1.sig <- names(pv1)[order(pv1)[1:nf]]
    egen <- unique(pv1.sig)
    return (egen)
}

#' This function takes pseudo-spatail and real-spatial data to identify variable genes
data_process <- function(st_count,st_label,anova){
    if (anova){
        sel.features <- select_feature(st_count[[1]],st_label[[1]])
        st_count_new <- list(st_count[[1]][sel.features,],st_count[[2]][sel.features,]) } else {
        st_count_new <- st_count }
    res1 <- normalize_data(st_count)
    st_norm <- res1[[1]]; variable_gene <- res1[[2]]; 
    st_scale <- scale_data(st_count_new,st_norm,variable_gene)
    return (list(st_count_new,st_norm,st_scale,variable_gene))
}

#' @param count.list list of pseudo-spatial data and real-spatial data, of which rows are genes and columns are cells
#' @param label.list list of pseudo-spatail label and real-spatial label (if any)
#' @return This function returns files saved in folders "Datadir" & "Infor_Data"
Convert_Data <- function(count.list,label.list,anova=TRUE){
    step1 <- data_process(st_count=count.list,st_label=label.list,anova)
    st.count <- step1[[1]];
    st.norm <- step1[[2]];
    st.scale <- step1[[3]];
    variable.genes <- step1[[4]]

    #' create data folders
    dir.create('Datadir'); dir.create('Output'); dir.create('DSTG_Result')
    inforDir <- 'Infor_Data'; dir.create(inforDir)
    
    #' save counts data to certain path: 'Datadir'
    write.csv(t(st.count[[1]]),file='Datadir/Pseudo_ST1.csv',quote=F,row.names=T)
    write.csv(t(st.count[[2]]),file='Datadir/Real_ST2.csv',quote=F,row.names=T)
    
    #' save scaled data to certain path: 'Infor_Data'
    write.csv(variable.genes,file=paste0(inforDir,'/Variable_features.csv'),quote=F,row.names=F)

    if (!dir.exists(paste0(inforDir,'/ST_count'))){dir.create(paste0(inforDir,'/ST_count'))}
    if (!dir.exists(paste0(inforDir,'/ST_label'))){dir.create(paste0(inforDir,'/ST_label'))}
    if (!dir.exists(paste0(inforDir,'/ST_norm'))){dir.create(paste0(inforDir,'/ST_norm'))}
    if (!dir.exists(paste0(inforDir,'/ST_scale'))){dir.create(paste0(inforDir,'/ST_scale'))}
    
    for (i in 1:2){
        write.csv(st.count[[i]],file=paste0(inforDir,'/ST_count/ST_count_',i,'.csv'),quote=F)
        write.csv(label.list[[i]],file=paste0(inforDir,'/ST_label/ST_label_',i,'.csv'),quote=F)
        write.csv(st.norm[[i]],file=paste0(inforDir,'/ST_norm/ST_norm_',i,'.csv'),quote=F)
        write.csv(st.scale[[i]],file=paste0(inforDir,'/ST_scale/ST_scale_',i,'.csv'),quote=F)
    }
}





