#!/bin/bash
#SBATCH -n 1            # number of cores (8 cores = 1 physical node)
#SBATCH -p compute      # compute queue
#SBATCH -t 1-00:10:00   # time (ddd-hh:mm:ss)
#SBATCH -J Basic_Ising
#SBATCH --reservation=application # optionally use the reservation (only 7 nodes)
rm -f temp_res128-Basic_ising.txt

#Do it for a number of temperatures
for i in $(seq 1.5 0.05 3 ) 
do
	echo "Calculating for temp  = $i"
	./ISING -t $i -n 100000 -m 64 -p >> temp_res128-Basic_ising.txt
done

