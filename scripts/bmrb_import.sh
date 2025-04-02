#!/bin/bash

IFS=$'\n'
lines=($( cat ../test_data/shiftx2_normalised_filenames.txt))

echo number lines: ${#lines[@]}

line_number=1

IFS=" \n\t"
for line in "${lines[@]}"
do
   read -r pdb chain bmrb res <<< "$line"

   if test "$pdb" = "#" || test "$pdb" = "---"
   then
     echo skipped $line
     line_number=$(expr $line_number + 1)
     continue
   fi



   file_name="../test_data/raw/bmrb${bmrb}_${pdb}_${chain}_raw.nef"

   if { [ ! -f "${file_name}" ] || [ ! -s "${file_name}" ]; }; then
     echo [${line_number}/${#lines[@]}] "PDB: $pdb, RES: $res, BMRB: $bmrb, CHAIN: $chain"
     echo File not found: "${file_name}", downloading
       nefl nmrstar import project $bmrb                                   \
     | nefl frames rename $bmrb default                                    \
     | nefl shiftx2 import shifts $pdb --source-chain $chain --chain A     \
     > "${file_name}"
   fi

  line_number=$(expr $line_number + 1)
done


