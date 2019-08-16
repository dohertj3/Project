#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

// Function to Run the Ising Model
void Ising(double** mat, int nsteps, int m, double t, double &res);

// Boundary conditions
int bound(int n, int m);

// Function to print Matrix
void print_mat(double** mat, int n, int m);

int main(int argc, char* argv[]){
	// Size of Matrix
	int n = 10;
	
	// Temperature of Matrix
	double temp = 2.3;

	// Number of steps
	int nsteps = 1000;

	// print in a format readable by my plotting script
	int p=0;

	// Use getopt to parse command line arguments
	char opt;

	while((opt = getopt(argc, argv, "nmtp"))!= -1){
		switch(opt){
			case 'n':
				nsteps = atoi(argv[optind]);
				break;
			case 'm':
				n = atoi(argv[optind]);
				break;
			case 't':
				temp = atof(argv[optind]);
				break;
			case 'p':
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

	double res;

	// Now start Ising model
	Ising(mat, nsteps, n, temp, res);

}



void Ising(double** mat, int nsteps, int m, double t, double &res){

	
	int energy, flip_energy, d_energy;

	res = 0;

	// Now iterate through grid nsteps times;
	for(int iter=0; iter<nsteps; iter++){
		for(int i=0; i<m; i++){
			for(int j=0; j<m; j++){
				energy = -1* mat[i][j] * ( mat[bound(i+1, m)][bound(j, m)] + mat[bound(i-1, m)][bound(j, m)] + mat[bound(i, m)][bound(j+1, m)] + mat[bound(i, m)][bound(j-1, m)] );
				flip_energy = mat[i][j] * ( mat[bound(i+1, m)][bound(j, m)] + mat[bound(i-1, m)][bound(j, m)] + mat[bound(i, m)][bound(j+1, m)] + mat[bound(i, m)][bound(j-1, m)] );
				
				d_energy = flip_energy - energy;

				if(d_energy <=0){
					mat[i][j] *= -1;	
				} else if(drand48() < exp(-1 * d_energy / t)){
					mat[i][j] *= -1;
				}

			}
		}
	}

	// Calculate the global average
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			res += mat[i][j];
		}
	}

	res /= (m*m);

}

int bound(int n, int m){
	n = n % m + m % m;
	return n;

}
