#!/bin/bash

# Loop through all *.png files
for png in *.png; do
  # Remove the file extension to isolate the gene_symbol
  base_name="${png%.png}"

  # Check if the corresponding .info file exists
  if [ ! -f "${base_name}.info" ]; then
    # If the .info file does not exist, delete the .png file
    echo "Deleting $png as no corresponding .info file exists."
    rm "$png"
  fi
done
