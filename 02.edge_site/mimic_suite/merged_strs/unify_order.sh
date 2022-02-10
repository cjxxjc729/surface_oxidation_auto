#!/bin/bash
ls -l | grep cif  | awk -F ' ' '{print $9}' > file.list
home_dir=$(pwd)
nfile=$(cat file.list | wc -l)


for i in `seq 1 $nfile`
do
file=$(cat file.list | head -$i | tail -1)
echo "-------------------$file begins ---------------------"

prefix=$(echo $file | awk -F '.' '{print $1}' )
nline=$(cat $file | wc -l)
nline_1=$(echo "$nline-1" | bc -l)
nline_2=$(echo "$nline-2" | bc -l)
sed -i "${nline}s/O/S/" $file
sed -i "${nline_1}s/C/Si/" $file
sed -i "${nline_2}s/O/S/" $file


echo "-------------------$file ends ---------------------"
done


