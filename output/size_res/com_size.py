import matplotlib.pyplot as plt
import numpy as np
import sys


def decomp_2D(nprocs):
	
	x = int(np.sqrt(nprocs))
	

	while( (nprocs % x) != 0):
		x += 1

	y = nprocs / x

	return x, y


def per_area(nprocs, size):

	x, y  = decomp_2D(nprocs)

	xsize = np.zeros(x)
	ysize = np.zeros(y)

	area = np.zeros(nprocs)
	perimeter = np.zeros(nprocs)

	x_size = size / x
	x_rem = size % x
	
	y_size = size/y
	y_rem = size % y

	for i in range(x):
		xsize[i] = x_size
		if i < x_rem:
			xsize[i] += 1

	for i in range(y):
		ysize[i] = y_size
		if i < y_rem:
			ysize[i] += 1

	for i in range(y):
		for j in range(x):
			area[j + i*x] = xsize[j] * ysize[i]

			if(ysize.size != 1):
				perimeter[j + i*x] += 2*xsize[j]

			if(xsize.size != 1):
				perimeter[j + i*x] += 2*ysize[i]

	return area, perimeter

def main():

	nprocs = np.array([2, 4, 8, 10, 15, 17, 20, 23, 25, 31, 35, 43, 50])



	# Load in results for file
	f = open("ParSFC_ising-neighbour.txt", 'r')

	serial_res = np.zeros(nprocs.size)

	par_res = np.zeros(nprocs.size)


	for i in range(nprocs.size):

		hil_size = np.zeros(nprocs[i])

		#calculate the size of each area in Hilbert curve
		div = 64 / nprocs[i]
		rem = 64 % nprocs[i]

		for j in range(nprocs[i]):
			hil_size[j] = div
			if j < rem:
				hil_size[j] += 1
	

		# Calculate perimeter for regular decomp method	
		size , perimeter = per_area(nprocs[i], 64)

		serial_res[i] = np.amax(size) - np.amin(size)

		current_line = f.readline()

		array = np.fromstring(current_line[:-1], dtype=np.int, sep=' ')
	
		par_res[i] = np.amax(hil_size) - np.amin(hil_size)

	plt.title("Average perimeter vs number of processors")
	plt.xlabel("Number of processors")
	plt.ylabel("Average perimenter of each processor")
	plt.scatter(nprocs, serial_res, label = "Basic")
	plt.scatter(nprocs, par_res, color = 'red', label = "SFC")

	plt.legend()

	plt.show()

	
	
	



if __name__ == "__main__":
	main()
