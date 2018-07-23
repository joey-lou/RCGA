#!/bin/bash
g++ RCGA.c -o RCGA

for i in 0 1 2
do
	for j in 0 1
	do
		for k in 0 1
		do
			echo "$i,$j,$k"
			./RCGA $i $j $k 0
			sleep 1
			./RCGA $i $j $k 0
			sleep 1
			./RCGA $i $j $k 0
		done
	done
done
