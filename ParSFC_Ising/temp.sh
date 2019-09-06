#!/bin/bash
#SBATCH -n 8             # number of cores (8 cores = 1 physical node)
#SBATCH -p compute      # compute queue
#SBATCH -t 10:00:00   # time (ddd-hh:mm:ss)
#SBATCH -J ParSFC_Ising
rm -f temp_res128-ParSFC_ising.txt

module load cports openmpi

#Do it for a number of temperatures
for i in $(seq 1.5 0.05 3)
do
	echo "Number of processors  = $i"
	mpirun -n 8 ./ising -t $i -n 100000 -m 128 -p -c >> temp_res128-ParSFC_ising.txt

done

