# BigWig - rename patient identifiers

# Description

Renamed patient identifiers to match them to discretized RNASeq data identifiers - removed suffix of each patient identifier.

- location:
    
    ```bash
    /projects/pulmonary_hypertension/work/clinical_data_stichit_run/segmentation_bigWig_input/bigWig_files_shortened_idents/
    ```
    

# File overview

- input_1
    - path:
        
        ```r
        ../bigWig_files/
        ```
        
- output
    - BigWig files
        
        ```bash
        bigWig_files_shortened_idents/*.bw 
        ```
        
    
- script - rename_bigWig_filenames_3.sh
```
# substring that has to be removed
remove_substring="_sorted"

# Loop over files
for file in *bw; do
	# Create new filename
	new_filename="${file/$remove_substring/}"
	# Rename file
	mv "$file" "$new_filename"
done
```