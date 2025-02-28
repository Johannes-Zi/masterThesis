# BigWig - reduced dataset

# Description

excluded bigWig files from patient id, which were not represented in the RNA Seq expression tables

- location:
    
    ```bash
    /projects/pulmonary_hypertension/work/clinical_data_stichit_run/segmentation_bigWig_input/bigWig_files_shortened_idents_reduced_dataset/
    ```
    

# File overview

- input
    - shortened idents bigWig files
        
        ```r
        ../bigWig_files_shortened_idents/*.bw 
        ```
        
    - output
        - reduced bigWig dataset
            
            ```bash
            bigWig_files_shortened_idents_reduced_dataset/
            ```
            
- script - manual removed:
```
healthy-7129.bw, pah-6785.bw, ph-heart.bw.bw, pah-7167.bw, pah-7193.bw, pah-7203.bw, pah-7209.bw, pah-7213.bw, pah-7236.bw, ph-heart-6750.bw, ph-heart-7246.bw, ph-heart-7247.bw, ph-lung-7228.bw, cteph-6649.bw, cteph-6923.bw, cteph-7195.bw, cteph-7220.bw, cteph-7231.bw
```