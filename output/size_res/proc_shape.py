import matplotlib.pyplot as plt
import numpy as np
import sys


def decomp_2D(nprocs):
	
	x = int(np.sqrt(nprocs))
	

	while( (nprocs % x) != 0):
		x += 1

	y = nprocs / x

	return x, y



def main():
	# Load results For Hilbert curve
	results = np.loadtxt(open(sys.argv[1]))

	nprocs = int(results[0, results.shape[0]-1]) + 1

	fig = plt.figure()

	plt.subplot(1, 2, 1)

	plt.xlim(0, 64)
	plt.ylim(0, 64)
	plt.title("Hilbert Decomp")

	plt.pcolormesh(results)

	plt.subplot(1, 2, 2)


	# Create standard division
	size = 64

	matrix = np.zeros((size, size))

	x, y = decomp_2D(nprocs)

	x_div = size/x
	x_rem = size % x

	y_div = size/y
	y_rem = size % y

	y_proc = -1 * x_div
	x_proc = -1

	y_size = 0
	x_size = 0

	for i in range(size):
		if(y_size == i):
			y_proc += x_div
			y_size += y_div
			if (y_proc < y_rem*x_div):
				y_size += 1

		for j in range(size):
			if(x_size == j):
				x_proc += 1
				x_size += x_div
				if (x_proc < x_rem):
					x_size += 1

			matrix[i][j] = x_proc + y_proc

		x_size = 0
		x_proc = -1

	plt.xlim(0, 64)
	plt.ylim(0, 64)

	plt.title("Basic Decomp")

	plt.pcolormesh(matrix)

	fig.suptitle("Decomposition of 64 * 64 grid into {} processors".format(nprocs))

	plt.show()





if __name__ == "__main__":
	main()
