import matplotlib.pyplot as plt
import numpy as np
import sys

def main():
	# Load in Serial results
	Basic_res = np.loadtxt("size_res-basic.txt");
	SFC_res = np.loadtxt("size_res-SFC.txt");
	
	# Load in ParSFC results
	SFC_8 = np.loadtxt("proc_res8-ParSFC_ising.txt") 
	SFC_16 = np.loadtxt("proc_res16-ParSFC_ising.txt") 
	SFC_32 = np.loadtxt("proc_res32-ParSFC_ising.txt") 
	SFC_64 = np.loadtxt("proc_res64-ParSFC_ising.txt") 

	# Load in check results
	check_8 = np.loadtxt("proc_res8-check_ising.txt") 
	check_16 = np.loadtxt("proc_res16-check_ising.txt") 
	check_32 = np.loadtxt("proc_res32-check_ising.txt") 
	check_64 = np.loadtxt("proc_res64-check_ising.txt") 

	# Load in basic time values
	basic_time = Basic_res[:, 3]
	SFC_time = SFC_res[:, 3]

	# Load in Parallel time values
	SFC8_times = np.transpose(SFC_8[:, 5])
	SFC16_times = np.transpose(SFC_16[:, 5])
	SFC32_times = np.transpose(SFC_32[:, 5])
	SFC64_times = np.transpose(SFC_64[:, 5])

	check8_times = np.transpose(check_8[:, 5])
	check16_times = np.transpose(check_16[:, 5])
	check32_times = np.transpose(check_32[:, 5])
	check64_times = np.transpose(check_64[:, 5])

	nprocs8 = SFC_8[:, 3]
	nprocs16 = SFC_16[:, 3]
	nprocs32 = SFC_32[:, 3]
	nprocs64 = SFC_64[:, 3]


	# Create graphs and values
	plt.title("Speedup vs number of processors(Check)")
	plt.xlabel("Number of processors")
	plt.ylabel("Speedup")

	plt.scatter(nprocs64, SFC_time[3]/SFC64_times, color = 'red', label = "64")
	plt.plot(nprocs64, SFC_time[3]/SFC64_times, color = 'red', label = '_nolegend_')

	plt.scatter(nprocs32, SFC_time[2]/SFC32_times, color = 'green', label = "32")
	plt.plot(nprocs32, SFC_time[2]/SFC32_times, color = 'green', label = '_nolegend_')

	plt.scatter(nprocs16, SFC_time[1]/SFC16_times, color = 'blue', label = "16")
	plt.plot(nprocs16, SFC_time[1]/SFC16_times, color = 'blue', label = '_nolegend_')

	plt.scatter(nprocs8, SFC_time[0]/SFC8_times, color = 'pink', label = "8")
	plt.plot(nprocs8, SFC_time[0]/SFC8_times, color = 'pink', label = '_nolegend_')

	plt.legend(title = "size")
	
	plt.plot((0, 50), (1, 1), color = 'orange')

	plt.show()
	

if __name__ == "__main__":
	main()
