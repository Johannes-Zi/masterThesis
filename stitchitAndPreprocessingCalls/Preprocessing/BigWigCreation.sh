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