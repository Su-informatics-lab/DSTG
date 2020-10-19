#' This functions takes sythetic data to test DSTG's performance

#' @return This function returns files saved in folders "Datadir" & "Infor_Data"
#' @export: all files are saved in current path
#' @examples: load data from folder "syntheic_data"

source('R_utils.R')
#' run sythetic data 
synthetic.count <- readRDS('./synthetic_data/example_data.RDS')
synthetic.label <- readRDS('./synthetic_data/example_label.RDS')

Convert_Data(synthetic.count,synthetic.label,anova=FALSE)

#' if you have the scRNA-seq data and want to decompose ST data, plese run below:

