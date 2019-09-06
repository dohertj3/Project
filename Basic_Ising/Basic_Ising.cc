#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <chrono>

using namespace std::chrono;

/* Function which calculates the Ising model different temperatures and size
 */


// Function to Run the Ising Model
void Ising(double** mat, int nsteps, int m, double t, double* res);

// Boundary conditions
int bound(int n, int m);

// Function to print Matrix
void print_mat(double** mat, int n, int m);

// Function send sides
void send_side(double** mat, int n);

int main(int argc, char* argv[]){
	// Size of Matrix
	int n = 10;
	
	// Temperature of Matrix
	double temp = 2.3;

	// Number of steps
	int nsteps = 10;

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

	// Input handling
	if(n < 4){
		printf("This function calculates the Ising model.\n Flags change the follow:\n");
		printf("t - the temperature of the lattice\n");
		printf("n - the number of iterations\n");
		printf("m - the size of the lattice\n");
		printf("p - prints in a way readable by certain functions\n");
		printf("\nm value must be greater that 4\n");
		return 0;
	
	}
	if(temp < 0){
		printf("This function calculates the Ising model.\n Flags change the follow:\n");
		printf("t - the temperature of the lattice\n");
		printf("n - the number of iterations\n");
		printf("m - the size of the lattice\n");
		printf("p - prints in a way readable by certain functions\n");
		printf("\nt value must be greater that 0\n");
		return 0;
	
	}
	if(nsteps < 1){
		printf("This function calculates the Ising model.\n Flags change the follow:\n");
		printf("t - the temperature of the lattice\n");
		printf("n - the number of iterations\n");
		printf("m - the size of the lattice\n");
		printf("p - prints in a way readable by certain functions\n");
		printf("\nn value must be greater that 1\n");
		return 0;
	
	}



	// Array to store results
	double* res = new double[nsteps]();

	// Create Matrix
	double** mat = new double*[n+2];
	double* array = new double[n*(n+2)];

	// Move index along
	mat++;

	// Initialise array
	for(int i=0; i<n; i++){
		mat[i] = &array[2 + i*(n+2)];
	}

	// Create boundary conditions row wise
	mat[-1] = mat[n-1];
	mat[n] = mat[0];

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

	// test
	send_side(mat, n);

	// Now start Ising model
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	Ising(mat, nsteps, n, temp, res);

	high_resolution_clock::time_point t2 = high_resolution_clock::now();


	duration<double> func_time = duration_cast<duration<double>>(t2-t1);

	// P flag is used to output results in a format so to be analysed plotting functions
	if(p == 1){
		printf(" %f %d %d ", temp, nsteps, n);
		printf(" %f ", func_time.count());

		for(int i=0; i<nsteps; i++){
			printf(" %f", res[i]);
		}
		printf("\n");
	} else {
		printf("Temp = %f, nsteps =  %d, size = %d, time = %f\n", temp, nsteps, n, func_time.count());
		printf("Magnetization of last Configuration = %f\n", res[nsteps-1]);
		

	}

	// Clean up
	mat--;
	delete[] res;
	delete[] mat;
	delete[] array;

}


/* Function which used to calculate the Ising model
 * mat, matrix which holds the lattice of spins
 * nsteps, number of iterations to go through
 * m, size of the lattice
 * t, temperature of the matrix
 * res, array which stores the results
 *
 */
void Ising(double** mat, int nsteps, int m, double t, double* res){
	
	double energy, flip_energy, d_energy;

	// Now iterate through grid nsteps times;
	for(int iter=0; iter<nsteps; iter++){
		for(int i=0; i<m; i++){
			for(int j=0; j<m; j++){

				// Calculate the energy
				energy = -1* mat[i][j] * ( mat[i+1][j] + mat[i-1][j] + mat[i][j+1] + mat[i][j-1] );

				// Propose change
				flip_energy = mat[i][j] * ( mat[i-1][j] + mat[i+1][j] + mat[i][j-1] + mat[i][j+1] );
				

				// Calculate change in energy
				d_energy = flip_energy - energy;

				// Decide to accept change
				if(d_energy <=0){
					mat[i][j] *= -1;	
				} else if(drand48() < exp(-1 * d_energy / t)){
					mat[i][j] *= -1;
				}

			}
		}

		// Update edges with new values
		send_side(mat, m);


	
		// Calculate the global average
		for(int i=0; i<m; i++){
			for(int j=0; j<m; j++){
				res[iter] += mat[i][j];
			}
		}
		res[iter] /= (m*m);

	}
}

// Function used for testing purposes
void print_mat(double** mat, int n, int m){
	for(int i=-1; i<n+1; i++){
		for(int j=-1; j<m+1; j++){
			printf(" %f", mat[i][j]);
		}
		printf("\n");
	}

}

// Updates edges with new values
void send_side(double** mat, int n){
	for(int i =0; i<n; i++){
		mat[i][-1] = mat[i][n-1];
		mat[i][n] = mat[i][0];
	}
}
