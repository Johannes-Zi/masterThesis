# Regression leave one out cross validation 


# Description

Performs regression based on output of segmentation.

Version with cross validation in form of leave one out outer cross validation instead of using a random subset for the cross validation.

This was done due to the small size of the avaliable clinical dataset - excluding multiple samples for the cross validation probably schrinks down the training dataset to much

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
        /projects/pulmonary_hypertension/work/analysis/clinical_data_stichit_run/regression_leave_one_out_CV/regression_output/
        ```
        
    - perform_regression.sh.log (by me - saves complete stdout)
    - regression_output/
        - *.RData (**Elasticnet Models** for each gene)
            - *_Pearson.RData  (**9059** files)
            - *_Spearman.RData  (**9277** files)
        - *.bed (genomic coordinates/ **OLS Models**)
            - Script to check  if all OLS coeff are zero
                
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
                
            - *_Pearson.bed (**9057** files)
            - *_Spearman.bed (**9277** files)
        - failed_genes.log (contains information about regression runs of genes that failed)
            - **failed executions** (**2958 total files**, **2395 uniqe genes** represented)
                - failed *Segmentation_<gene_id>_??_Pearson.txt (**1521 files**)
                - failed *Segmentation_<gene_id>_??_Spearman.txt (**1437 files**)
            - **error types**
                - **2110** ‘x should be a matrix with 2 or more columns’
                - **664** ‘x has missing values; consider using makeX() to impute them’
                - 0 ‘elastic net model contains only zero coefficients’’
                - **25** ‘Error: from glmnet C++ code (error code 7777); All used predictors have zero variance’
                - **43** ‘Error in elnet(xd, is.sparse, y, weights, offset, type.gaussian, alpha, : y is constant; gaussian glmnet fails at standardization step’
                - **43** ‘y is constant; gaussian glmnet fails at standardization step’
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
    
- script
    - Download Two_level_learning_elasticNet.R
        
        ```bash
        # Downloaded on 14.03.2024; commit dde55a8
        wget https://raw.githubusercontent.com/SchulzLab/SingleCellStitchit/main/bin/Two_level_learning_elasticnet.R
        ```
        
    
    ```bash
    # Define parameters
    input_data_dir="./combined_segmentation_output/"        # Combined output of segmentation step
    output_data_dir="./regression_output/"          # ...
    response_param="Expression"             # Expression column in segmentation file (Stichit sprecific)
    core_param=32                   # Performance
    alpha_param=0.05                # Elastic net regression - grid search value for regularization
    textsize_param=0.2              # percentage for testdataset of outer cross validation
    regularization_param="E"                # Elastic net or eg. L for lasso
    innerCV_param=6                 # Number of folds for inner cross validation
    outerCV_param=6                 # ...
    performance_param="TRUE"                # Performance
    leaveOneOutCV_param="TRUE"              # Specific type of cross validation
    asRData_param="FALSE"           # irrelevant old param
    randomise_param="FALSE"         # not relevant
    logResponse_param="TRUE"                # log transformation
    ftest_param="FALSE"             # ...
    coefP_param=1                   # filter for pvalue of coefficients
    
    # Regression Script call
    Rscript 'Two_level_learning_elasticnet.R' --dataDir=$input_data_dir --outDir=$output_data_dir --response=$response_param --cores=$core_param --alpha=$alpha_param --testsize=$textsize_param --regularisation=$regularization_param --innerCV=$innerCV_param --outerCV=$outerCV_param --performance=$performance_param --leaveOneOutCV=$leaveOneOutCV_param --asRData=$asRData_param  --randomise=$randomise_param  --logResponse=$logResponse_param --ftest=$ftest_param  --coefP=$coefP_param
    ```
    

# Detailed workflow

For bugfixes and script adaption see documentation of regression version with normal cross validation

Command to run regression

```jsx
nohup bash perform_regression.sh > perform_regression.sh.log 2>&1 &
```