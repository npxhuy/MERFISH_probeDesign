#!/bin/bash

temp_file="temp.txt"

# Loop over every directory in the current directory
ls -d */ | while read dir; do
    # Remove the trailing slash from the directory name
    dir=${dir%/}

    # Construct the file name
    file="${dir}/${dir}_gene_id.txt"

    # Check if the file exists
    if [[ -f "$file" ]]; then
        # Remove lines containing "No_Match" and lines where column 7 is 0.0
        awk -F'\t' '$7 == "0.0"' "$file" | grep -v "No_Match" > "$temp_file"

        # Count the number of lines in the file
        num_lines=$(wc -l < "$temp_file")

        # If the file has more than 51 lines
        if (( num_lines > 51 )); then
            # Remove the first line and randomly select 50 lines
            tail -n +2 "$temp_file" | shuf -n 50 > "$file"
        # If the file has between 21 and 51 lines (inclusive)
        elif (( num_lines >= 21 && num_lines <= 51 )); then
            # Remove the first line and take all lines
            tail -n +2 "$temp_file" > "$file"
        # If the file has less than 21 lines
        else
            # Discard the file
            rm "$file"
        fi
    fi
done

# Remove the temporary file
rm "$temp_file"