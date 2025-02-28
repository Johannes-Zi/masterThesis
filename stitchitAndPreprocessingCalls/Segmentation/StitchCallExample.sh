# Initialize parameters
gene_list_file="sub_gene_list_1.txt"

bigwig_dir="../bigWig_files_shortened_idents_reduced_dataset/"  # Path to dir with epigenetic signal
genome_annotation="../gencode.v38.annotation.gtf"       # Path to genome annotation
discretized_expression="../sim_discretized_ge_shortened_idents.tsv"     # Path to discretized expression data
continuous_expression="../sim_tpm_normalized_ge_counts_shortened_idents.tsv"    # Path to continuous expression data
chromosome_sizes="../GCF_000001405.39_GRCh38.p13_chr_sizes.txt" # Path to file that holds sizes of the chromosomes
extension_size=2000     # Extension size up and downstream around target gene
number_of_cores=20      # Number of cores used by run
p_value=0.05    # p-value threshold for regulatory element selection
merge_resolution=10     # Resolution to merge initial data to reduce runtime
output_path="output_dir/stichit_output_dir/"     # Path of output files
max_searchspace=1100000 # Maximum size of search space
segment_size=2000       # Segment size

single_stichit_call() {
	# Set local variable for first handed over parameter
	local target_gene_id="$1"
	# Perform stichit call
	/projects/pulmonary_hypertension/work/STITCHIT/build/core/STITCH -b $bigwig_dir -a $genome_annotation -d $discretized_expression -o $continuous_expression -s $chromosome_sizes -w $extension_size -c $number_of_cores -p $p_value -g $target_gene_id -z $merge_resolution -f $output_path -r $max_searchspace -t $segment_size
}

summary_log_file="output_dir/summary_stichit_run_log_file.log"
loop_iterator=0		# Counts the number of loops that are already performed
# Iterate over gene ids represented in the RNA Seq expression table
for gene_id in $(head -n 6000 $gene_list_file); do 	# Processes up to tenthousand gene ids
	((loop_iterator++))	# Increase the loop iterator
	echo -e "\e[33mStichit run ${loop_iterator}\e[0m"
	printf "${gene_id}\n\n"
	# Capture output or crash and save it as log file
        log_file="output_dir/log_dir/log_${gene_id}.log"
        # Perform function call and capture output
        {
        single_stichit_call "$gene_id" || printf "\nstichit call for ${gene_id} crashed!\n"
        } > $log_file 2>&1	# Redirect stderr and stdout to log file
	echo >> $log_file
	cat $log_file	# Print log file content
	cat $log_file >> $summary_log_file	# Append logs to summary log file
done