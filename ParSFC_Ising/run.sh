
#!/bin/bash
#SBATCH -n 8            # number of cores (8 cores = 1 physical node)
#SBATCH -p compute      # compute queue
#SBATCH -t 1-00:10:00   # time (ddd-hh:mm:ss)
#SBATCH -J Par_SFC
#SBATCH --reservation=application # optionally use the reservation (only 7 nodes)
rm -f Par_SFC-temp_res.txt

module load cports openmpi

#Do it for a number of temperatures
for i in $(seq 1.5 0.1 4 ) 
do
	echo "Calculating for temp  = $i"
	mpi_run -n 4 .ising -t $i -n 100000 -m 64 -p -c >> Par_SFC-temp_res.txt
done

