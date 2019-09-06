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

def plot_stuff(files):

	# Load Basic results
	results = np.loadtxt(files);

	#Take results
	
	temps = results[:, 0]
	nsteps = results[0, 1]
	size = results[0, 2]

	# Get Magnetization results
	magnet_sam = results[:, 5:]

	# Remove all points not at equilibrium
	magnet_sam  = magnet_sam[:, 1000:]
	magnetization = np.mean(magnet_sam, axis = 1)


	
	# Calculate error
	uncor_val = magnet_sam[:, ::10]


	magnetization = np.mean(uncor_val, axis = 1)

	error = jackknife(uncor_val)

	return magnetization, error, temps

def main():
	# Load Basic results
	res_64, error_64, temps = plot_stuff("ParSFC_ising-temp_res.txt")
	res_32, error_32, temps = plot_stuff("ParSFC32_ising-temp_res.txt")
	res_16, error_16, temps = plot_stuff("ParSFC16_ising-temp_res.txt")
	res_8, error_8, temps = plot_stuff("ParSFC8_ising-temp_res.txt")


	# Plot results
	plt.figure()	
	plt.title("Magnetization vs Temperature for different matrix sizes")
	plt.xlabel("Temperature(kBT)")
	plt.ylabel("Magnetization")

	plt.scatter(temps, abs(res_64), color = 'red', label = "64")
	plt.errorbar(temps, abs(res_64), yerr = error_64, fmt = None, label = '_nolegend_')

	plt.scatter(temps, abs(res_32), color = 'green', label = "32")
	plt.errorbar(temps, abs(res_32), yerr = error_32, fmt = None, label = '_nolegend_')

	plt.scatter(temps, abs(res_16), color = 'blue', label = "16")
	plt.errorbar(temps, abs(res_16), yerr = error_16, fmt = None, label = '_nolegend_')

	plt.scatter(temps, abs(res_8), color = 'purple', label = "8")
	plt.errorbar(temps, abs(res_8), yerr = error_8, fmt = None, label = '_nolegend_')

	plt.legend(title = "Side Length")

	plt.show()

if __name__ == "__main__":
	main()
