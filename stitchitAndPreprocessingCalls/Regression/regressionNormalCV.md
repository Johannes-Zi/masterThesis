# Regression normal cross validation


# Description

Performs regression based on output of segmentation.

encountered an error - used updated Two_level_learning_elasticNet.R version

# File overview

- input
    - Output of segmentation in form of gene specific Pearson and Spearman correlations (11701 genes represented)
        - *_Pearson.txt (11291 files)
        - *_Spearman.txt (11367 files)
        
        ```r
        /projects/pulmonary_hypertension/work/analysis/clinical_data_stichit_run/segmentation/combined_segmentation_output
        ```
        
- output
    - path:
        
        ```bash
        /projects/pulmonary_hypertension/work/analysis/clinical_data_stichit_run/regression/regression_output/
        ```
        
    - perform_regression.sh.log (by me - saves complete stdout)
    - regression_output/
        - *.RData (**Elasticnet Models** for each gene)
            
            detailed file strucutre overview:
            
            [Files Structure](https://www.notion.so/Files-Structure-9ff7d501490948d69422cb345a35be19?pvs=21)
            
            - *_Pearson.RData  (**9159** files)
            - *_Spearman.RData  (**9326** files)
        - *.bed (genomic coordinates/ **OLS Models**)
            - Script to check  if all coeff OLS are zero
                
                ```bash
                # Commandline input
                ols_coeff_file_path="$1"
                
                # Check if all coeffincients are zero
                all_coeff_are_zero=true
                
                # Iterate over the lines
                while IFS=$'\t' read -r col1 col2 col3 col4 colX
                do
                        # Check if coordinate line
                        if [[ "$col1" == "chr"* ]]; then
                                # Check if coefficient is 0
                                if [[ "$col4" != 0 ]]; then
                                        all_coeff_are_zero=false
                                        break
                                fi
                        fi
                done < $ols_coeff_file_path
                
                # Retrun result
                echo $all_coeff_are_zero
                ```
                
            - *_Pearson.bed (**9152** files)
                - all OLS coefficients are zero - **2767**
                - at least one OLS coeff ≠ zero - **6385**
            - *_Spearman.bed (**9320** files)
                - all OLS coefficients are zero - **2826**
                - at least one OLS coeff ≠ zero - **6494**
        - failed_genes.log (contains information about regression runs of genes that failed)
            - **failed executions** (**2958 total files**, **2395 uniqe genes** represented)
                - failed *Segmentation_<gene_id>_??_Pearson.txt (**1521 files**)
                - failed *Segmentation_<gene_id>_??_Spearman.txt (**1437 files**)
            - **error types**
                - **2134** ‘x should be a matrix with 2 or more columns’
                - **729** ‘x has missing values; consider using makeX() to impute them’
                - **8** ‘elastic net model contains only zero coefficients’’
                - **14** ‘Error: from glmnet C++ code (error code 7777); All used predictors have zero variance’
                - **23** ‘Error in elnet(xd, is.sparse, y, weights, offset, type.gaussian, alpha, : y is constant; gaussian glmnet fails at standardization step’
                - **57** ‘y is constant; gaussian glmnet fails at standardization step’
        - Performance_Overview.txt
            - overview correlation measures
                
                
                | Pearson r
                (follows distribution) | PearsonVar () | Spearman (non-parametric - follows no distribution) | SpearmanVar R^2 ??? | MSE (Mean Squared Error) | MSEVar | pVal | qVal |  |
                | --- | --- | --- | --- | --- | --- | --- | --- | --- |
                |  measures the strength and direction of the linear relationship between two variables.
                
                1 indicates a perfect positive linear relationship,
                    -1 indicates a perfect negative linear relationship, and
                    0 indicates no linear relationship. |  | The correlation coefficient, denoted as "ρ" (rho), measures the strength and direction of the monotonic relationship between two variables. It ranges from -1 to 1:
                
                    1 indicates a perfect positive monotonic relationship,
                    -1 indicates a perfect negative monotonic relationship, and
                    0 indicates no monotonic relationship. |  | measure used to evaluate the performance of a predictive model, particularly in regression analysis. It quantifies the average squared difference between the actual values (observed outcomes) and the predicted values produced by the model.
                
                A lower MSE indicates better agreement between the predicted and actual values, implying better model performance. |  | The p-value indicates the probability of observing the correlation coefficient (or one more extreme) under the null hypothesis that there is no correlation between the variables. It helps assess the statistical significance of the observed correlation coefficient. | False Discovery Rate (FDR): The FDR is the expected proportion of false discoveries among the rejected hypotheses. In other words, it's the rate at which you're willing to tolerate false positives.
                
                q-value: The q-value is a measure of the FDR. It represents the minimum FDR at which a hypothesis may be considered significant. In practical terms, a q-value of 0.05 indicates that you expect no more than 5% of the rejected hypotheses to be false positives. |  |
                |      |  |  |  |  |  |  |  |  |
            
    - Correlation coefficients
    
- script
    - Download Two_level_learning_elasticNet.R
        
        ```bash
        # Downloaded on 14.03.2024; commit dde55a8
        wget https://raw.githubusercontent.com/SchulzLab/SingleCellStitchit/main/bin/Two_level_learning_elasticnet.R
        ```
        
    
    ```bash
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
    
    ```
    

# Detailed workflow

1. Adaption of minimum number of input files for regression
    
    [Regression - minimum file number adaption](https://www.notion.so/Regression-minimum-file-number-adaption-cc4ef762df844673aa031d4cad460f98?pvs=21)
    
2. Installed missing R packages
    1. glmnet 
    2. doMC 1.3.8.
3. Error and updated Two_level_learning_elasticNet.R (error_1/)
    
    I encountered an error when running the Two_level_learning_elasticNet.R  script
    
    - error:
        
        ```bash
        Call:  cv.glmnet(x = x_train, y = y_train, nfolds = as.numeric(argsL$innerCV),      keep = TRUE, alpha = x) 
        
        Measure: Mean-Squared Error 
        
            Lambda Index Measure     SE Nonzero
        min  0.583    12   1.151 0.4000       3
        1se  1.622     1   1.225 0.3339       0
        [1] 26
        [1] 4
        [1] "Learning sample  872"
        [1] "Processing sample matrix. can take a few minutes. Please wait."
        [1] "/projects/pulmonary_hypertension/work/analysis/clinical_data_stichit_run/segmentation/combined_segmentation_output/Segmentation_ENSG00000086289_10_Pearson.txt"
        [1] 49
        [1] "Outer cross validation fold:  1"
        [1] "constant response y in innerCV -> skip sample"
        [1] "Learning OLS model on reduced feature space for sample 872"
        NULL
        integer(0)
        Error in best_models[[min.cv.err.ind]] : 
          attempt to select less than one element in get1index
        In addition: There were 50 or more warnings (use warnings() to see the first 50)
        Execution halted
        ```
        
    
    Regarding to Laura, the error occurred because the program hasn’t captured a skipped regression for a single gene. The regression was skipped, because the inner cross validation generated a constant regression vector for one of the folds (happends when the expression in a fold is constant - e.g. 0 - is a result of the sparsed dataset). The folds are subsets of the samples. Laura suggested to skip a problematic fold exactly once, if it appears for a gene, to reduce the number of genes that get lost due to this constant expression in the fold. But the selection of the folds is done with the Monte Carlo strategy(random folds), thus the occurance of the error is also random. Also the errors occur not only at the constant fold expression, but also there is sometimes no variance in the predictors or missing values in the table.
    
    Laura: “The scaling procedure was also corrected - the scaling is now performed only on the training data for each outer CV fold then the scaling parameters from the training data to scale the test data instead of scaling the whole feature matrix before starting the CV-procedure”
    
4. Additional error fix (error_2/)
    
    ```r
    Error in model.frame.default(formula = Expression ~ ., data = ols_Data,  : 
      'data' must be a data.frame, environment, or list
    Calls: lm -> eval -> eval -> <Anonymous> -> model.frame.default
    In addition: There were 50 or more warnings (use warnings() to see the first 50)
    Execution halted
    
    ```
    
    - Details:
        
        @laurarumpf2 I was able to reproduce the same error by performing 1000 iterations with ENSG00000133597 and ENSG00000133612.
        
        The input, output and the scripts are located at:
        
        /projects/pulmonary_hypertension/work/analysis/clinical_data_stichit_run/regression/error_2/error_reproduction/
        
        The stdout log of the affected iteration is:
        
        execution_iteration_873_stdout.log
        
        The whole stdout log is:
        
        loop_regression_[command.sh](http://command.sh/).log
        
        @johanneszieres1 thanks for the log file. I was able to identify the error: we only have 5 regions that were selected as features with the segmentation for the gene “ENSG00000133597“. After generating the elastic net model we have no non-zero coefficients, i.e. we could not learn anything from the input data with only 25 samples. I caught this error and also write a corresponding entry to the error log file.
        
5. Run regression
    
    ```jsx
    nohup bash perform_regression.sh > perform_regression.sh.log 2>&1 &
    ```