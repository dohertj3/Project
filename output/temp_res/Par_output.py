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
	results = np.loadtxt("check_ising-temp_res.txt");

	#Take results
	
	temps = results[:, 0]
	nsteps = results[0, 1]
	size = results[0, 2]

	# Get Magnetization results
	magnet_sam = results[:, 5:]

	# Remove all points not at equilibrium
	magnet_sam  = magnet_sam[:, 100:]
	magnetization = np.mean(magnet_sam, axis = 1)

	# Calculate some correlations
	tot_tau = 30

	correlations15 = np.zeros(tot_tau)
	correlations22 = np.zeros(tot_tau)
	correlations23 = np.zeros(tot_tau)
	correlations24 = np.zeros(tot_tau)
	correlations30 = np.zeros(tot_tau)


	for i in range(1, tot_tau + 1):
		correlations15[i-1] = correlation(magnet_sam[0, :], i)
		correlations22[i-1] = correlation(magnet_sam[16, :], i)
		correlations23[i-1] = correlation(magnet_sam[17, :], i)
		correlations24[i-1] = correlation(magnet_sam[18, :], i)
		correlations30[i-1] = correlation(magnet_sam[31, :], i)

	
	# Calculate error
	uncor_val = magnet_sam[:, ::10]


	magnetization = np.mean(uncor_val, axis = 1)

	error = jackknife(uncor_val)

	# Plot results
	plt.scatter(range(1, tot_tau + 1), correlations15, color = 'red', label = "SFC t = 1.5")
	plt.scatter(range(1, tot_tau + 1), correlations22, color = 'pink', label = "SFC t = 2.2")
	plt.scatter(range(1, tot_tau + 1), correlations23, color = 'green', label = "SFC t = 2.3")
	plt.scatter(range(1, tot_tau + 1), correlations24, color = 'blue', label = "SFC t = 2.4")
	plt.scatter(range(1, tot_tau + 1), correlations30, color = 'magenta', label = "SFC t = 3.0")

	plt.legend()

	plt.figure()	
	plt.scatter(temps, abs(magnetization))
	plt.errorbar(temps, abs(magnetization), yerr = error, fmt = None)
	plt.show()

if __name__ == "__main__":
	main()
