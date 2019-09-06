#!/bin/bash
#SBATCH -n 128            # number of cores (8 cores = 1 physical node)
#SBATCH -p compute      # compute queue
#SBATCH -t 1-00:10:00   # time (ddd-hh:mm:ss)
#SBATCH -J Check_Ising
rm -f check_ising-proc_res.txt

module load cports openmpi

#Do it for a number of temperatures
for i in 2 4 8 10 16 20 32 64 128 
do
	echo "Number of processors  = $i"
	mpirun -n $i ./ising -t 2 -n 100000 -m 64 -p >> check_ising-proc_res.txt
done

