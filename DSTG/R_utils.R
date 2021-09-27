#' evalaute performance
JSD_performance <- function(spots_true_composition, spots_predicted_composition) 
{
    suppressMessages(require(philentropy))
    jsd_matrix <- matrix(nrow = nrow(spots_true_composition), ncol = 1)
                            
    for (i in seq_len(nrow(spots_true_composition))) {
        x <- rbind(spots_true_composition[i, ], spots_predicted_composition[i, ])
        if (sum(spots_predicted_composition[i, ]) > 0) {
            jsd_matrix[i, 1] <- suppressMessages(JSD(x = x, unit = "log2", est.prob = "empirical"))
        }
        else { jsd_matrix[i, 1] <- 1 } }
    
    quants_jsd <- round(quantile(matrixStats::rowMins(
                                                  jsd_matrix, 
                                                  na.rm = TRUE), c(0.25, 0.5, 0.75)), 5)
    cat(sprintf("The following JSD quantiles are obtained:
               %s [%s - %s]", 
              quants_jsd[[2]], quants_jsd[[1]], quants_jsd[[3]]), sep = "\n")
    return(list(JSD = jsd_matrix))
}

#' normalize function
normalize_data <- function(count.list){
    norm.list <- vector('list')
    var.features <- vector('list')
    for ( i in 1:length(count.list)){
        norm.list[[i]] <- as.matrix(Seurat:::NormalizeData.default(count.list[[i]],verbose=FALSE))
        hvf.info <- Seurat:::FindVariableFeatures.default(count.list[[i]],selection.method='vst',verbose=FALSE)
        hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 0), ]
        hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
        var.features[[i]] <- head(rownames(hvf.info), n = 2000)
    }
    sel.features <- selectIntegrationFeature(count.list,var.features)
    return (list(norm.list,sel.features))}

#' scaling function 
scale_data <- function(count.list,norm.list,hvg.features){
    scale.list <- lapply(norm.list,function(mat){
        Seurat:::ScaleData.default(object = mat, features = hvg.features,verbose=FALSE)})
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
        st_count_new <- list(st_count[[1]][sel.features,],st_count[[2]][sel.features,])

        colnames(st_label[[1]]) <- 'subclass'
        tem.t1 <- Seurat::CreateSeuratObject(counts = st_count_new[[1]],meta.data=st_label[[1]]);
        Seurat::Idents(object = tem.t1) <- tem.t1@meta.data$subclass

        #' convert scRNA-seq data to pseudo-spatial data                                                                                                                          
        test.spot.ls1<-SPOTlight::test_spot_fun(se_obj=tem.t1,clust_vr='subclass',n=1000);
        test.spot.counts1 <- as.matrix(test.spot.ls1[[1]])
        colnames(test.spot.counts1)<-paste("mixt",1:ncol(test.spot.counts1),sep="_");
        metadata1 <- test.spot.ls1[[2]]
        test.spot.metadata1 <- do.call(rbind,lapply(1:nrow(metadata1),function(i){metadata1[i,]/sum(metadata1[i,])}))
        st_counts <- list(test.spot.counts1,st_count_new[[2]])

        st_label[[1]] <- test.spot.metadata1
        N1 <- ncol(st_counts[[1]]); N2 <- ncol(st_counts[[2]])
        label.list2 <- do.call("rbind", rep(list(st_label[[1]]), round(N2/N1)+1))[1:N2,]
        st_labels <- list(st_label[[1]],label.list2)
    } else {
        st_counts <- st_count; st_labels=st_label }
    res1 <- normalize_data(st_counts)
    st_norm <- res1[[1]]; variable_gene <- res1[[2]];
    st_scale <- scale_data(st_counts,st_norm,variable_gene)
    return (list(st_counts,st_labels,st_norm,st_scale,variable_gene))
}

#' @param count.list list of pseudo-spatial data and real-spatial data, of which rows are genes and columns are cells
#' @param label.list list of pseudo-spatail label and real-spatial label (if any)
#' @return This function returns files saved in folders "Datadir" & "Infor_Data"
Convert_Data <- function(count.list,label.list,anova=TRUE){
    step1 <- data_process(st_count=count.list,st_label=label.list,anova)
    st.count <- step1[[1]];
    st.label <- step1[[2]];
    st.norm <- step1[[3]];
    st.scale <- step1[[4]];
    variable.genes <- step1[[5]]

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
        write.csv(st.label[[i]],file=paste0(inforDir,'/ST_label/ST_label_',i,'.csv'),quote=F)
        write.csv(st.norm[[i]],file=paste0(inforDir,'/ST_norm/ST_norm_',i,'.csv'),quote=F)
        write.csv(st.scale[[i]],file=paste0(inforDir,'/ST_scale/ST_scale_',i,'.csv'),quote=F)
    }
}





