#!/bin/bash


for file in ./data2/*.mol2
do
    for ((i = 1; i < 11; i++))
    do
    echo "WORKING"

	for ((s = 2; s < 21; s++))
		do
			name=${file##*/}
			base=${name%.mol2}
			for((g = 4;g < 7;g++))
			do 
				python3 main.py True "$file"  "${base}#s${s}#it${i}#g${g}insp"  ${i} ${s}  False  True  True  ${g}  2  &
				echo "Launched ${base}#s${s}#it${i}#g${g}insp"
			done
			wait
	done
    done
done
echo "SCRIPT IS DONE w simp"





