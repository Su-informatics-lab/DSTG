
#' check JSD score
true <- read.csv('./DSTG_Result/true_output.csv',header=F)
predict <- read.csv('./DSTG_Result/predict_output.csv',header=F)

source('R_utils.R')
jsd.score <- synthetic_performance(
    test_spots_metadata_mtrx = as.matrix(true),
    spot_composition_mtrx = as.matrix(predict))



