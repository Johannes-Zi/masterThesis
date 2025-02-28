library(glmnet)
library(ggplot2)


# Path to the segmentation output file
segmentation_path_1 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000125629_10/Segmentation_ENSG00000125629_10_Pearson.txt"
# Path to the elastic net model file
elnet_path_1 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000125629_10/Elasticnet_Regression_Model_Segmentation_ENSG00000125629_10_Pearson.RData"

segmentation_path_2 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000156232_10/Segmentation_ENSG00000156232_10_Pearson.txt"
elnet_path_2 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000156232_10/Elasticnet_Regression_Model_Segmentation_ENSG00000156232_10_Pearson.RData"

segmentation_path_3 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000249092_10/Segmentation_ENSG00000249092_10_Pearson.txt"
elnet_path_3 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000249092_10/Elasticnet_Regression_Model_Segmentation_ENSG00000249092_10_Pearson.RData"

segmentation_path_4 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000259056_10/Segmentation_ENSG00000259056_10_Pearson.txt"
elnet_path_4 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000259056_10/Elasticnet_Regression_Model_Segmentation_ENSG00000259056_10_Pearson.RData"

segmentation_path_5 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000211788_10/Segmentation_ENSG00000211788_10_Pearson.txt"
elnet_path_5 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000211788_10/Elasticnet_Regression_Model_Segmentation_ENSG00000211788_10_Pearson.RData"


load(elnet_path_4)
segmentation_data<-read.table(segmentation_path_4,header=TRUE,sep="",row.names=1)
print(segmentation_data)

segmentation_data <- unique(segmentation_data)
write.table(segmentation_data, file = "segmentation_normalized_ENSG00000259056.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

segmentation_data_df <- data.frame(segmentation_data)

# Log2 transformation to normalize the data (+1 to avoid log of 0)
segmentation_data_df <- log2(segmentation_data_df+1)

segmentation_data_df <- data.frame(scale(segmentation_data_df,center=TRUE, scale=TRUE))
print(segmentation_data_df)

write.table(segmentation_data_df, file = "segmentation_scaled_centered_ENSG00000259056.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Removes RNAseq column from the dataframe
segmentation_data_df_normalized <- segmentation_data_df[, -ncol(segmentation_data_df)]

lambda_min <- elasticnet_model$model$lambda.min
#coefficients <- coef(elasticnet_model$model, s = elasticnet_model$model$lambda.min)

predictions <- predict(elasticnet_model$model, 
                       newx = as.matrix(segmentation_data_df_normalized),
                       s = lambda_min)
print(predictions)

rna_seq_data = segmentation_data_df[, ncol(segmentation_data_df)]

df_combined = data.frame(rna_seq_data, predictions)
print(df_combined)


ggplot(df_combined, aes(x = rna_seq_data, y = predictions)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    xlab("normalized and scaled original expression") +
    ylab("normalized and scaled predicted expression") +
    ggtitle("Predicted vs Original Expression Values ENSG00000211788") +
    theme(plot.title = element_text(hjust = 0.5))

ggsave("PredictedVsActualExpression_ENSG00000211788.png", plot = last_plot(),
       device = "png", path = getwd(), width = 5.5, height = 6.5)



extract_elasticnet_model_coefficients <- function(elasticnet_file_path){
  load(elasticnet_file_path)
  print(elasticnet_model$model$lambda.min)
  coefficients <- coef(elasticnet_model$model, s = elasticnet_model$model$lambda.min)
  print(coefficients)
  return(coefficients)
}

elasticnet_path <- "C:/Users/johan/Desktop/local_master_thesis_data/regression/LOneOCV_regression/regression_output/Elasticnet_Regression_Model_Segmentation_ENSG00000002587_10_Pearson.RData"
extract_elasticnet_model_coefficients(elasticnet_path)
