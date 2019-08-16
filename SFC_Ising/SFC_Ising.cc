#include "sfc.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Function to print a matrix
void print_mat(double** mat, int n, int m);

// Function which runs a metropolis algorithm to calculate the Ising model
// n = number of iterations
// t = temp
// res = average magnetization at the end
void SFC_Ising(SFC* a, int n, double t, double &res);

int main(int argc, char* argv[]){
	// Size of Matrix
	int n = 4;
	int rank;

	// Temperature of matrix
	double temp = 0.2;

	int nsteps = 1000;

	int p = 0;	

	char opt;

	// Use getopt to take command line arguments
	while((opt = getopt(argc, argv, "nmtp"))!= -1){
		switch(opt){
			case 'n':
				nsteps = atoi(argv[optind]);
				break;

			case'm':
				n = atoi(argv[optind]);
				break;

			case't':
				temp = atof(argv[optind]);
				break;

			case'p':
				p = 1;
				break;

		}
	}

	// Create Matrix
	double** mat = new double*[n];	
	double* array = new double[n*n];

	for(int i=0; i<n; i++){
		mat[i] = &array[i*n];
	}

	srand48((int)time(NULL));

	// Randomly set values
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(drand48() < .5){
				mat[i][j] = 1;
			} else {
				mat[i][j] = -1;
			}	
		}
	}


	rank = (int) (log(n)/log(2));

	hilbert test(mat, rank , n, n);

	// Now start Metropolis algorithm
	double res;
	SFC_Ising(&test, nsteps, temp, res);

	if(p == 0){
		printf(" Global Average = %f, Temp = %f, nsteps = %d, size = %d\n", res, temp, nsteps, n);
	} else {
		printf(" %f %f %d %d\n", res, temp, nsteps, n);
	}

	delete [] mat;
	delete [] array;

}

/* Metropolis function
 * This function goes updates the grid in the order of a SFC
 * n is the number of iterations
 *
 */
void SFC_Ising(SFC* a, int n, double t, double &res){

	// Initalise values
	double energy, flip_energy, d_energy;
	int x, y;
	int up, down, left, right;

	res = 0;

	// Now iterate through the grid

	for(int i = 0; i<n; i++){
		// Now iterate along SFC
		for(int j=0; j< a->mat_size(); j++){
			a->index_2D(j, x, y);

			a->index_1D(up, x, y+1);
			a->index_1D(down, x, y - 1);
			a->index_1D(left, x + 1, y);
			a->index_1D(right, x - 1, y);

			energy = -1 * (*a)[j] * ((*a)[up] + (*a)[down] + (*a)[left] + (*a)[right]);

			flip_energy = (*a)[j] * ((*a)[up] + (*a)[down] + (*a)[left] + (*a)[right]);
			
			d_energy = flip_energy - energy;	

			if(d_energy <= 0){
				a->set(-1 * (*a)[j], j);
			} else if(drand48() <  exp(-1 * d_energy / t)) {
				a->set(-1 * (*a)[j], j);
			}
		}


	}

	// Calculate total Magnetism
	for(int i = 0; i< a->mat_size(); i++){
		res += (*a)[i];
	}
	
	res /= a->mat_size();
}


// Function to print a matrix
void print_mat(double** mat, int n, int m){
	int i, j;

	printf("\n------------Printing matrix---------------\n");
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			printf(" %f", mat[i][j]);

		}
		printf("\n");
	}
	printf("\n");
}
