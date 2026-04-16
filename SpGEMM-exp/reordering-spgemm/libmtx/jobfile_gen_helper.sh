#!/bin/bash

input_file="$1"
output_file="$2"

if [[ -z "$input_file" || -z "$output_file" ]]; then
  echo "Usage: $0 <input_file> <output_file>"
  exit 1
fi

> "$output_file"  # Clear or create the output file

while IFS= read -r line; do
  if [[ "$line" =~ ^/global/homes/r/raqib/mtspgemmlib/bin/Shuffle_hw\ text\ /pscratch/sd/r/raqib/dataset-spgemm/([^/]+)/([^/]+)\.mtx$ ]]; then
    dir="${BASH_REMATCH[1]}"
    echo "echo \"~~~~~~~~~~~~~~~~~~~~~~~~~<$dir>~~~~~~~~~~~~~~~~~~~~~~~~~\"" >> "$output_file"
    echo "$line" >> "$output_file"
    echo "" >> "$output_file"
  else
    echo "$line" >> "$output_file"
  fi
done < "$input_file"
