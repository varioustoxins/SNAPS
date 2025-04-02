#!/bin/zsh


SNAPS_DIR="${HOME}/Dropbox/git/SNAPS_garyt"
SNAPS=$SNAPS_DIR/python/SNAPS.py
DATA=$SNAPS_DIR/test_data

## statistical analysis
header="BMRB_code PDB chain bad total %bad"
echo $header
for file in bmrb*_snaps_out.nef;
do

  val=$("${SNAPS_DIR}"/scripts/snaps_nef_statistics.sh "$file")

  bmrb_code_PDB_chain=$(echo $file | tr '_' ' '| awk '{print $1, $2, $3}')
  length=$(echo $val | awk '{print $2}')
  bad=$(echo $val | awk '{print $5}')
  percentage=$(echo "scale=2; ($bad / $length) * 100" | bc)

  echo  $bmrb_code_PDB_chain $bad $length $percentage

done