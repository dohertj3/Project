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
	results_Basic = np.loadtxt("Basic_ising-temp_res.txt");
	results_SFC = np.loadtxt("SFC_ising-temp_res.txt");
	results_Check = np.loadtxt("check_ising-temp_res.txt");
	results_Par = np.loadtxt("ParSFC_ising-temp_res.txt");

	# Get Magnetization results
	sam_basic = results_Basic[:, 4:]
	sam_sfc = results_SFC[:, 4:]
	sam_check = results_Check[:, 5:]
	sam_par = results_Par[:, 5:]

	# Remove all points not at equilibrium
	sam_basic  = sam_basic[:, 1000:]
	sam_sfc  = sam_sfc[:, 1000:]
	sam_check  = sam_check[:, 1000:]
	sam_par  = sam_par[:, 1000:]

	magnet_basic = np.mean(sam_basic, axis = 1)
	magnet_sfc = np.mean(sam_sfc, axis = 1)
	magnet_check = np.mean(sam_check, axis = 1)
	magnet_par = np.mean(sam_par, axis = 1)

	# Calculate some correlations
	tot_tau = 30

	cor_basic = np.zeros(tot_tau)
	cor_sfc = np.zeros(tot_tau)
	cor_check = np.zeros(tot_tau)
	cor_par = np.zeros(tot_tau)

	temp = 16

	for i in range(1, tot_tau + 1):
		cor_basic[i-1] = correlation(sam_basic[temp, :], i)
		cor_sfc[i-1] = correlation(sam_sfc[temp, :], i)
		cor_check[i-1] = correlation(sam_check[temp, :], i)
		cor_par[i-1] = correlation(sam_par[temp, :], i)

	
	# Plot results

	plt.title("The Autocorrelation for different programs at t = 2")
	plt.xlabel("Tau, the step length")
	plt.ylabel("Autocorrelation")

	plt.scatter(range(1, tot_tau + 1), cor_basic, color = 'red', label = "Basic")
	plt.scatter(range(1, tot_tau + 1), cor_sfc, color = 'pink', label = "SFC")
	plt.scatter(range(1, tot_tau + 1), cor_check, color = 'green', label = "Check")
	plt.scatter(range(1, tot_tau + 1), cor_par, color = 'blue', label = "Par SFC")

	plt.legend()

	plt.show()

if __name__ == "__main__":
	main()
