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
	results = np.loadtxt(open(sys.argv[1]));

	#Take results
	temps = results[:, 0]
	nsteps = results[0, 1]
	size = results[0, 2]

	# Get Magnetization results
	magnet_sam = results[:, 4:]

	# Remove all points not at equilibrium
	magnet_sam  = magnet_sam[:, 1000:]
	magnetization = np.mean(magnet_sam, axis = 1)

	# Calculate some correlations
	tot_tau = 30

	correlations15 = np.zeros(tot_tau)
	correlations225 = np.zeros(tot_tau)
	correlations20 = np.zeros(tot_tau)
	correlations25 = np.zeros(tot_tau)
	correlations30 = np.zeros(tot_tau)


	for i in range(1, tot_tau + 1):
		correlations15[i-1] = correlation(magnet_sam[0, :], i)
		correlations225[i-1] = correlation(magnet_sam[16, :], i)
		correlations20[i-1] = correlation(magnet_sam[11, :], i)
		correlations25[i-1] = correlation(magnet_sam[21, :], i)
		correlations30[i-1] = correlation(magnet_sam[31, :], i)

	
	# Calculate error
	uncor_val = magnet_sam[:, ::10]


	magnetization = np.mean(uncor_val, axis = 1)

	error = jackknife(uncor_val)


	# Plot results

	plt.title("The Autocorrelation for a number of steplengths at varying temperatures")
	plt.xlabel("Tau, the step length")
	plt.ylabel("Autocorrelation")

	plt.scatter(range(1, tot_tau + 1), correlations15, color = 'red', label = "SFC t = 1.5")
	plt.scatter(range(1, tot_tau + 1), correlations225, color = 'pink', label = "SFC t = 2.25")
	plt.scatter(range(1, tot_tau + 1), correlations20, color = 'green', label = "SFC t = 2.0")
	plt.scatter(range(1, tot_tau + 1), correlations25, color = 'blue', label = "SFC t = 2.5")
	plt.scatter(range(1, tot_tau + 1), correlations30, color = 'magenta', label = "SFC t = 3.0")

	plt.legend()


	plt.figure()	

	plt.title("Magnetization vs Temperature")
	plt.xlabel("Temperature(kBT)")
	plt.ylabel("Magnetization")

	plt.scatter(temps, abs(magnetization))
	plt.errorbar(temps, abs(magnetization), yerr = error, fmt = None)
	plt.show()

if __name__ == "__main__":
	main()
