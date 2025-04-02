cd ../python

filename=$(basename -- "$1" | sed -e 's/_raw//')              #  test/wibble.flap
extension="${filename##*.}"                                   #  flap or .flap
filename_root="${filename%.*}"                                #  wibble
directory=$(dirname   "$1")

mkdir "${directory}_out" &> /dev/null
log_file="${directory}_out/${filename_root}_log.txt"
out_file="${directory}_out/${filename_root}_snaps_out.nef"
echo log in $log_file
echo output in $out_file
./SNAPS.py                                                          \
--shift_type nef                                                    \
"$1"                                                                \
shiftx2                                                             \
"${out_file}"                                                       \
-c ../config/config_yaml.txt                                        \
-l ${log_file}                                                      \
> "$out_file"

