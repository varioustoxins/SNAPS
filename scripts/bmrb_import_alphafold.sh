#!/bin/bash

IFS=$'\n'
lines=($( cat ../test_data/shiftx2_normalised_filenames.txt ))

echo number lines: ${#lines[@]}

line_number=1

IFS=" \n\t"
for line in "${lines[@]}"
do
   read -r pdb chain bmrb res alphafold <<< "$line"

   if test "$pdb" = "#" || test "$pdb" = "---"
   then
     line_number=$(expr $line_number + 1)
     continue
   fi

   echo [${line_number}/${#lines[@]}] "PDB: $pdb, RES: $res, BMRB: $bmrb, CHAIN: $chain"

   file_name="../test_data/alphafold_raw/bmrb${bmrb}_${pdb}_${chain}_alphafold_raw.nef"
   echo "$file_name"
   if [ ! -f "${file_name}" ] && [ "$alphafold" = "yes" ];  then
     echo File not found: "${file_name}" downloading
       nefl nmrstar import project $bmrb                                                       \
     | nefl frames rename $bmrb default                                                        \
     | nefl shiftx2 import shifts $pdb --source-chain $chain --chain A  --alphafold  --verbose \
     > "${file_name}"
   fi

  line_number=$(expr $line_number + 1)
done
