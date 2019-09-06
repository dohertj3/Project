import matplotlib.pyplot as plt
import numpy as np
import sys

def jackknife(matrix):
	x, y = matrix.shape

	averages = np.zeros((x, y))

	tot_sum = np.zeros(x)
	
	for i in range(x):
		tot_sum[i] = np.sum(matrix[i, :])

	for i in range(x):
		for j in range(y):
			averages[i, j] = (tot_sum[i] - matrix[i, j]) / (y - 1)


#	for i in range(x):
#		for j in range(y):
#			 averages[i, j] = np.mean(matrix[i, j:]) 

	variances = np.var(averages, axis = 1)

	return variances

def correlation(array, tau):

	array_av = np.mean(array)

	sum_tau = 0	

	sum_av = 0

	for i in range(array.size - tau - 1):
		sum_tau += (array[i] - array_av) * (array[i+tau] - array_av)

	for i in range(array.size):
		sum_av += (array[i] - array_av)**2

	return sum_tau/sum_av

def main():
	# Load Basic results
	Basic128 = np.loadtxt("temp_res128-Basic_ising.txt");
	SFC128 = np.loadtxt("temp_res128-SFC_ising.txt");
	check128 = np.loadtxt("temp_res128-check_ising.txt");
	par128 = np.loadtxt("temp_res128-ParSFC_ising.txt");

	#Take results
	temps = Basic128[:, 0]
	nsteps = Basic128[0, 1]

	# Get All magnetization results
	tot_basic = Basic128[:, 4:]
	tot_SFC = SFC128[:, 4:]
	tot_check = check128[:, 6:]
	tot_par = par128[:, 6:]
	
	# Get Magnetization results
	magnet_basic = Basic128[:, 1000:]
	magnet_SFC = SFC128[:, 1000:]
	magnet_check = check128[:, 1000:]
	magnet_par = par128[:, 1000:]


	# Calculate some correlations
	tot_tau = 30

	cor_basic = np.zeros(tot_tau)
	cor_SFC = np.zeros(tot_tau)
	cor_check = np.zeros(tot_tau)
	cor_Par = np.zeros(tot_tau)


	for i in range(1, tot_tau + 1):
		cor_basic[i-1] = correlation(magnet_basic[16, :], i)
		cor_SFC[i-1] = correlation(magnet_SFC[16, :], i)
		cor_check[i-1] = correlation(magnet_check[16, :], i)
		cor_Par[i-1] = correlation(magnet_par[16, :], i)

	
	# Calculate error
	uncor_basic = magnet_basic[:, ::10]
	uncor_SFC = magnet_check[:, ::10]
	uncor_check = magnet_check[:, ::10]
	uncor_Par = magnet_par[:, ::10]


#	magnetization = np.mean(uncor_val, axis = 1)
	magz_basic = np.mean(uncor_basic, axis = 1)
	magz_SFC = np.mean(uncor_SFC, axis = 1)
	magz_check = np.mean(uncor_check, axis = 1)
	magz_Par = np.mean(uncor_Par, axis = 1) 

	error_basic = jackknife(uncor_basic)
	error_SFC = jackknife(uncor_SFC)
	error_check = jackknife(uncor_check)
	error_Par = jackknife(uncor_Par)


	# Plot results
	# Plot correlation 

	plt.title("The Autocorrelation for a number of Implementations t = 2.25kbT")
	plt.xlabel("Tau, the step length")
	plt.ylabel("Autocorrelation")

	plt.scatter(range(1, tot_tau + 1), cor_basic, color = 'red', label = "Basic")
	plt.scatter(range(1, tot_tau + 1), cor_SFC, color = 'pink', label = "SFC")
	plt.scatter(range(1, tot_tau + 1), cor_check, color = 'green', label = "Check")
	plt.scatter(range(1, tot_tau + 1), cor_Par, color = 'blue', label = "ParSFC")

	plt.legend()

	# Plot Magnetisation 
	plt.figure()	

	plt.title("Magnetization vs Temperature")
	plt.xlabel("Temperature(kBT)")
	plt.ylabel("Magnetization")

	plt.scatter(temps, abs(magz_basic), color = 'red', label = "Basic")
	plt.errorbar(temps, abs(magz_basic), yerr = error_basic, fmt = None, label = '_nolegend_')

	plt.scatter(temps, abs(magz_SFC), color = 'pink', label = "SFC")
	plt.errorbar(temps, abs(magz_SFC), yerr = error_SFC, fmt = None, label = '_nolegend_')

	plt.scatter(temps, abs(magz_check), color = 'green', label = "check")
	plt.errorbar(temps, abs(magz_check), yerr = error_check, fmt = None, label = '_nolegend_')

	plt.scatter(temps, abs(magz_Par), color = 'blue', label = "ParSFC")
	plt.errorbar(temps, abs(magz_Par), yerr = error_Par, fmt = None, label = '_nolegend_')

	plt.legend();

	#Plot of all magnetizations
	plt.figure()

	plt.title("Magnetisation at each sample")
	plt.xlabel("Iteration")
	plt.ylabel("Magneisation Samples")
	
	plt.scatter(range(0, tot_par.shape[1], 40), tot_par[10, ::40], color = 'green', label = '2')
	plt.scatter(range(0, tot_par.shape[1], 40), tot_par[16, ::40], color = 'blue', label = '2.25')
	plt.scatter(range(0, tot_par.shape[1], 40), tot_par[21, ::40], color = 'pink', label = '2.5')
	plt.scatter(range(0, tot_par.shape[1], 40), tot_par[25, ::40], color = 'magenta', label = '3')
	plt.scatter(range(0, tot_par.shape[1], 40), tot_par[13, ::40], color = 'orange', label = '2.1')

	plt.legend(title = 'Temp(kbT)')

	plt.show()

if __name__ == "__main__":
	main()
