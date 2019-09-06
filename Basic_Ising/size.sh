#!/bin/bash
#SBATCH -n 1            # number of cores (8 cores = 1 physical node)
#SBATCH -p compute      # compute queue
#SBATCH -t 08:00:00   # time (ddd-hh:mm:ss)
#SBATCH -J Basic_Ising
#SBATCH --reservation=application # optionally use the reservation (only 7 nodes)
rm -f size_res-basic.txt

#Do it for a number of temperatures
for i in 8 16 32 64 128 
do
	echo "Calculating for temp  = $i"
	./ISING -t 2 -n 100000 -m $i -p >> size_res-basic.txt
done

