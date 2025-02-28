# substring that has to be removed
remove_substring="_sorted"

# Loop over files
for file in *bw; do
	# Create new filename
	new_filename="${file/$remove_substring/}"
	# Rename file
	mv "$file" "$new_filename"
done