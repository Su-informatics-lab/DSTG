
#' check JSD score
true <- read.csv('./DSTG_Result/true_output.csv',header=F)
predict <- read.csv('./DSTG_Result/predict_output.csv',header=F)

source('R_utils.R')

jsd.score <- JSD_performance(
    spots_true_composition = as.matrix(true),
    spots_predicted_composition = as.matrix(predict))



