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