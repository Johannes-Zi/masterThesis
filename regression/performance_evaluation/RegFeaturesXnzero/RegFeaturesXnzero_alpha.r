
input_file_path <- paste(
  "C:/Users/johan/Desktop/standard_regression/regression_output/regression_output/",
  "Elasticnet_Regression_Model_Segmentation_ENSG00000001167_10_Spearman.RData",
  sep = ""
)

cat(sprintf("File path for import:\n%s\n", input_file_path))

#load(file.choose())
load(input_file_path)

lambda_min_value <- elasticnet_model$model$lambda.min
cat("lambda min of the model:", lambda_min_value)

lambda_min_index <- which(elasticnet_model$model$lambda == lambda_min_value)
cat("lambda min index:", lambda_min_index)

# corresponding nzero value
nzero_value <- elasticnet_model$model$nzero[lambda_min_index]
number_of_features <- elasticnet_model$model$glmnet.fit$dim[1]

print("overview for lambda.min")
cat("Number of features used in the model:", number_of_features, "\n")
cat("Number of features with coefficients unequal to 0:", nzero_value, "\n")

