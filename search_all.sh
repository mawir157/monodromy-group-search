#!/bin/bash
array=( 3 4 5 6 7 8 9 10 12 14 15 18 20 21 24 27 28 30 32 36 40 42 48 54 60 90 )
for i in "${array[@]}"
do
	for j in "${array[@]}"
	do
		if [ "$i" -lt "$j" ]
		then
			./mono -l $i $j
		fi
	done
done
