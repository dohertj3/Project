#!/bin/bash
#SBATCH -n 23            # number of cores (8 cores = 1 physical node)
#SBATCH -p compute      # compute queue
#SBATCH -t 02:00:00   # time (ddd-hh:mm:ss)
#SBATCH -J check_Ising
#SBATCH --reservation=application # optionally use the reservation (only 7 nodes)
rm -f proc_res16-check_ising.txt

module load cports openmpi

#Do it for a number of temperatures
for i in 2 4 8 10 15 17 20 23
do
	echo "Number of processors  = $i"
	mpirun -n $i ./ising -t 2 -n 100000 -m 16 -p -c >> proc_res16-check_ising.txt
done

