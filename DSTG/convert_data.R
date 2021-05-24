#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' This functions takes sythetic data to test DSTG's performance

#' @return This function returns files saved in folders "Datadir" & "Infor_Data"
#' @export: all files are saved in current path
#' @examples: load data from folder "syntheic_data"

source('R_utils.R')
#' if you have the scRNA-seq data and want to decompose ST data, plese run below:

if (length(args)==0) {
    message('run synthetic data...')
    synthetic.count <- readRDS('./synthetic_data/example_data.RDS')
    synthetic.label <- readRDS('./synthetic_data/example_label.RDS')
    Convert_Data(synthetic.count,synthetic.label,anova=FALSE)
} else if (length(args)==3) {  
    message('run real data...')
    sc.count <- readRDS(args[1])
    st.count <- readRDS(args[2])
    intersect.genes <- intersect(rownames(sc.count),rownames(st.count))
    sc.count <- sc.count[intersect.genes,]
    st.count <- st.count[intersect.genes,]
    count.list <- list(sc.count,st.count)
    label.list <- list(data.frame(readRDS(args[3]),stringsAsFactors=F))
    Convert_Data(count.list,label.list,anova=TRUE)
}

