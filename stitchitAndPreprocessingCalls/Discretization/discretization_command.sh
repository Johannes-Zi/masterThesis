#Set parameters
input_expression_file_path="counts.matrix.norm.tsv"		# have to add tsv to it - or transform to csv
output_dir="./discretization_output/"
genome_annotation="gencode.v38.annotation.gtf"
output_dist_ge_matrix="${output_dir}sim_discretized_ge.tsv"	# Discretized gene expression marix
output_normalized_matrix="${output_dir}sim_tpm_normalized_ge_counts.tsv"	# Normalized(tpm, DESeq2) matrix (optional)
output_gene_lengths="${output_dir}sim_gene_lengths.txt"	# Gene lengths (optional)
# -n flag for normalized; -p for plot

# Discretization call (-b for binary option(0,1 instead of 0,1,2) is not used)
Rscript 'gene_expression_discretization.R' -f $input_expression_file_path -n -g $genome_annotation -p -o $output_dist_ge_matrix -t $output_normalized_matrix -l $output_gene_lengths