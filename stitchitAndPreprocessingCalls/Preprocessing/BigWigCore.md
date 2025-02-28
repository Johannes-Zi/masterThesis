# BigWig creation - core process

# Description

Sorted and index the .bam ATAC Seq files (initially provided index files were not matching to the -bam files). Transformed them to bigWig files for the segmenetation as input.

- location:
    
    ```bash
    /projects/pulmonary_hypertension/work/clinical_data_stichit_run/segmentation_bigWig_input/
    ```
    

# File overview

- input
    - ATAC Seq *.bam files
        
        ```bash
        /projects/pulmonary_hypertension/work/data/raw_data/ATAC/*.bam
        ```
        
    
- output
    - sorted_indexed_files/
        
        ```bash
        _sorted.bam
        _sorted.bam.bai
        ```
        
    - bigWig_files/
        
        ```bash
        _sorted.bw
        ```
        
- script
    - Sort and index files - sort_index_command.sh
        
        ```bash
        # input directory
        input_dir="/projects/pulmonary_hypertension/work/data/raw_data/ATAC/"
        output_dir="./sorted_indexed_files/"
        
        # Iterate over files in input directory
        for input_file in "$input_dir"*.bam
        do
        	# Check if file
        	if [ -f "$input_file" ]; then
        		# Extract filename
        		temp_filename=$(basename "$input_file" .bam)
        		# Create sorted file filepath
        		output_file="${output_dir}${temp_filename}_sorted.bam"
        		echo "input: ${input_file}"
        		echo "output sorted: ${output_file}"
        		# Create Sorted file
        		samtools sort $input_file -o $output_file
        		# Create index files for sorted files
        		samtools index $output_file
        	fi
        done
        ```
        
    - Create BigWig files - bigWig_creation_command.sh
    ```
    # input directory
    input_dir="./sorted_indexed_files/"
    output_dir="./bigWig_files/"

    # Iterate over files in input directory
    for input_file in "$input_dir"*.bam
    do
        # Check if file
        if [ -f "$input_file" ]; then
            temp_filename=$(basename "$input_file" .bam)
            output_file="${output_dir}${temp_filename}.bw"
            echo "input: ${input_file}"
            echo "output: ${output_file}"
            bamCoverage -b $input_file -o $output_file --normalizeUsing 'RPKM' --binSize 50 -p 5
        fi
    done
    ```