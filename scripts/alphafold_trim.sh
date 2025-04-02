#!/bin/bash

in_dir=../test_data/alphafold_raw
out_dir=../test_data/alphafold_trimmed

mkdir $out_dir
for file_path in  ${in_dir}/*.nef
do
   filename=$(basename -- "$file_path" | sed -e 's/_raw//')      #  test/wibble.flap
   extension="${filename##*.}"                                   #  flap or .flap
   filename_root="${filename%.*}"                                #  wibble
   directory=$(dirname   "$file_path")                           #  test

   new_filename="${out_dir}/${filename_root}_trimmed.${extension}"

   echo "Processing $file_path -> $new_filename"
     nefl stream $file_path           \
   | nefl chains align --verbose      \
   | nefl loops trim                  \
   > "${new_filename}"

done
