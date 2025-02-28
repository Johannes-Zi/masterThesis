# Define parameters
input_data_dir="./combined_segmentation_output/"	# Combined output of segmentation step
output_data_dir="./regression_output/"		# ...
response_param="Expression"		# Expression column in segmentation file (Stichit sprecific)
core_param=10			# Performance
alpha_param=0.05		# Elastic net regression - grid search value for regularization
textsize_param=0.2		# percentage for testdataset of outer cross validation
regularization_param="E"		# Elastic net or eg. L for lasso
innerCV_param=6			# Number of folds for inner cross validation
outerCV_param=6			# ...
performance_param="TRUE"		# Performance
leaveOneOutCV_param="FALSE"		# Specific type of cross validation
asRData_param="FALSE"		# irrelevant old param
randomise_param="FALSE"		# not relevant
logResponse_param="TRUE"		# log transformation
ftest_param="FALSE"		# ...
coefP_param=1			# filter for pvalue of coefficients

# Regression Script call
Rscript 'Two_level_learning_elasticNet.R' --dataDir=$input_data_dir --outDir=$output_data_dir --response=$response_param --cores=$core_param --alpha=$alpha_param --testsize=$textsize_param --regularisation=$regularization_param --innerCV=$innerCV_param --outerCV=$outerCV_param --performance=$performance_param --leaveOneOutCV=$leaveOneOutCV_param --asRData=$asRData_param  --randomise=$randomise_param  --logResponse=$logResponse_param --ftest=$ftest_param  --coefP=$coefP_param
