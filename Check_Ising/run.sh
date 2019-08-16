#!/bin/bash

# Function which runs the Ising model for a number of results
module load cports openmpi


rm -f temp_res.txt

# Do it for a number of temperature
for i in $(seq 1 0.1 4)
do 
	echo "Calculating for temp = $i"
	mpirun -n 4 ./ising -t $i -n 1000000 -m 32 -p >> temp_res.txt
	echo "finished"

done
