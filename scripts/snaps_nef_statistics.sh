#!/bin/bash


test=$(awk 'NF > 2 && $2 != "." {test = ($3-$5) ; test = sqrt(test^2)> 0 ; print test}' $1)


count=$( echo $test | count )
length=$( echo $test | add )
length=${length%.*}

echo length $count of bad $length



