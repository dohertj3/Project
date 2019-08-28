#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <chrono>

#include "mpi.h"

using namespace std::chrono;

/*
 * This script uses the check model to parallelise the Ising mode
 *
 */

void ising(double** mat, MPI_Comm comm, double temp, int step, int nsteps, int n, int m, double* measurements);

void print_rank(double** mat, int n, int m, int rank, int p_rank, MPI_Comm comm);

int main(int argc, char* argv[]){
	// Iterables
	int i, j;
	int nsteps = 100;
	
	// Value which tell a processor whether to start on evens or odds
	int step;
	int x_step = 0;
	int y_step = 0;
	double temp = 500;
	int length = 10;

	int p=0;

	int c = 0;


	// Use get opt to parse agtuments
	char opt;

	while((opt = getopt(argc, argv, "tnmpc"))!= -1){
		switch(opt){
			case 't':
				temp = atof(argv[optind]);
				break;

			case 'n':
				nsteps = atoi(argv[optind]);
				break;

			case 'm':
				length = atoi(argv[optind]);
				break;

			case 'p':
				p=1;
				break;
			case 'c':
				c=1;
				break;

		}
	}

	int rank, nprocs;

	int nsamples = nsteps;


	// begin MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


	MPI_Barrier(MPI_COMM_WORLD);
	high_resolution_clock::time_point setup_start = high_resolution_clock::now();

	// Find dimensions of cartesian communicator
	int* psize = new int[2]();

	MPI_Dims_create(nprocs, 2, psize);

	// Create boundary conditions
	int bound[3];
	bound[0] = 1;
	bound[1] = 1;
	bound[2] = 1;

	// Create cartehamstring strainsian co-ordinates
	MPI_Comm cart_comm;
	MPI_Cart_create(MPI_COMM_WORLD, 2, psize, bound, 0, &cart_comm);

	// Calculate the number of points on each processor
	int y_length;
	int x_length;

	int y_remainder = length % psize[0];
	int y_div = length / psize[0];

	int x_remainder = length % psize[1];
	int x_div = length / psize[1];

	// Get cartesion co-ordinate
	int coords[2];

	MPI_Cart_coords(cart_comm, rank, 2, coords);

	// Loop through each processor, this needs to be done in a loop so each processor 
	for(i=0; i<psize[1]; i++){

		x_length = x_div;

		// If processor rank is less than remainder, add an extra point
		if(i < x_remainder){
			x_length++;
		}
		
		// Finish loop if done your processor
		if(coords[1] == i){
			break;
		}

		// If the length is odd, add one to step counter so it knows when to start and stop
		if(x_length % 2 == 1){
			x_step++;
		}
	}
	
	// Loop through each processor, this needs to be done in a loop so each processor knows the step if ither processors
	for(i=0; i<psize[0]; i++){

		y_length = y_div;

		// If processor rank is less than remainder, add an extra point
		if(i < y_remainder){
			y_length++;
		}

		// Finish loop if done your processor
		if(coords[0] == i){
			break;
		}

		// If the length is odd, add one to step counter so it knows when to start and stop
		if(y_length % 2 == 1){
			y_step++;
		}
		
	}
	step = x_step + y_step;
	step = step % 2;
	
	// Create matrix on each processor
	// Need extra rows and columns for halo points
	double** mat = new double*[y_length + 2]();
	double* array = new double[(y_length + 2)*(x_length + 2)]();


	// Now set up matrix
	// Move mat so that 
	for(i=0; i<(y_length+2); i++){
		mat[i] = &(array[i*(x_length+2)]) + 1;
	}

	srand48(time(NULL)+ rank);

	// Now move mat point so first row isn't in the halo;
	mat++;


	for(i=0; i<y_length; i++){
		for(j=0; j<x_length; j++){
			if( drand48() < 0.5){
				mat[i][j] = 1;
			} else {
				mat[i][j] = -1;
			}
		}
	}

	MPI_Barrier(cart_comm);
	high_resolution_clock::time_point setup_end = high_resolution_clock::now();

	duration<double> setup_time = duration_cast<duration<double>>(setup_end - setup_start);

	// Print x_length and y_length

	// Now begin the Ising model for each rank

	// Create array that will store the results
	double* measurements = new double[nsamples]();
	
	MPI_Barrier(cart_comm);	

	high_resolution_clock::time_point func_start = high_resolution_clock::now();

	ising(mat, cart_comm, temp, step, nsteps, y_length, x_length, measurements);

	MPI_Barrier(cart_comm);

	high_resolution_clock::time_point func_end = high_resolution_clock::now();

	duration<double> func_time = duration_cast<duration<double>>(func_end - func_start);

	// Now Gather local Average to create 
	double* global_av = new double[nsamples];

	MPI_Allreduce(measurements, global_av, nsamples, MPI_DOUBLE, MPI_SUM, cart_comm);

	// Now get the average
	for( i=0; i<nsamples; i++){
		global_av[i] /= length*length;	
	}

	if(rank == 0){
		if(p == 0){
			printf("Global average = %f, Temp = %f, nsteps = %d, length = %d, nprocs = %d\n", global_av[nsamples-1], temp, nsteps, length, nprocs);
		} else if(p == 1){
			printf("%f %d %d %d ", temp, nsteps, length, nprocs);

			printf("times= %f %f ", setup_time.count(), func_time.count());

			printf("results=");
			for(i=0; i<nsamples; i++){
				printf(" %f", global_av[i]);
			}
			printf("\n");

		}
	}

	// Move mat back so it can be properly deleted
	mat--;

	delete[] psize;

	delete[] mat;
	delete[] array;

	delete[] measurements;
	delete[] global_av;


	
	MPI_Finalize();
	return 0;
}


// Function which Calculuates the Ising model on each processor
void ising(double** mat, MPI_Comm comm, double temp, int step, int nsteps, int n, int m, double* measurements){
	int iter, i, j;

	double energy, flip_energy, d_energy;

	// First find which processors are top bottom left and right
	int tp, bp, lp, rp;
	int odd_even;

	int rank;
	MPI_Comm_rank(comm, &rank);

	MPI_Cart_shift( comm, 0, 1, &tp, &bp);
	MPI_Cart_shift( comm, 1, 1, &lp, &rp);

	// Create column
	MPI_Datatype column;
	MPI_Type_vector(n, 1, m+2, MPI_DOUBLE, &column);
	MPI_Type_commit(&column);
	MPI_Request* req;
	req = new MPI_Request[4];


	// Begin iterations;
	for( iter=0; iter<nsteps; iter++){

		odd_even = step;		

		// Send halo points
		// Send Top recv bottom
		MPI_Isend(&(mat[0][0]), m, MPI_DOUBLE, tp, 0, comm, &req[0]);
		MPI_Irecv(&(mat[n][0]), m, MPI_DOUBLE, bp, 0, comm, &req[0]);

		// Send bottom recv top
		MPI_Isend(&(mat[n-1][0]), m, MPI_DOUBLE, bp, 1, comm, &req[1]);
		MPI_Irecv(&(mat[-1][0]), m, MPI_DOUBLE, tp, 1, comm, &req[1]);

		// Send right recv left
		MPI_Isend(&(mat[0][m-1]), 1, column, rp, 2, comm, &req[2]);
		MPI_Irecv(&(mat[0][-1]), 1, column, lp, 2, comm, &req[2]);	

		// Send left recv right
		MPI_Isend(&(mat[0][0]), 1, column, lp, 3, comm, &req[3]);
		MPI_Irecv(&(mat[0][m]), 1, column, rp, 3, comm, &req[3]);

		MPI_Waitall(4, req, MPI_STATUS_IGNORE);


		// Now begin ising model
		for(i = 0; i<n; i++){
			for(j = odd_even; j<m; j+=2){

				energy = -1 * mat[i][j] * (mat[i+1][j] + mat[i-1][j] + mat[i][j-1] + mat[i][j+1]);
				d_energy = -2*energy;
				if(d_energy <= 0){
					mat[i][j] *= -1;
				} else if ( drand48() < exp(-1 *  d_energy / temp)){
					mat[i][j] *= -1;
				}
			}

			odd_even = !odd_even;

		}

		// Now Send halo points again
		// Send Top recv bottom
		MPI_Isend(&(mat[0][0]), m, MPI_DOUBLE, tp, 0, comm, &req[0]);
		MPI_Irecv(&(mat[n][0]), m, MPI_DOUBLE, bp, 0, comm, &req[0]);

		// Send bottom recv top
		MPI_Isend(&(mat[n-1][0]), m, MPI_DOUBLE, bp, 1, comm, &req[1]);
		MPI_Irecv(&(mat[-1][0]), m, MPI_DOUBLE, tp, 1, comm, &req[1]);

		// Send right recv left
		MPI_Isend(&(mat[0][m-1]), 1, column, rp, 2, comm, &req[2]);
		MPI_Irecv(&(mat[0][-1]), 1, column, lp, 2, comm, &req[2]);

		// Send left recv right
		MPI_Isend(&(mat[0][0]), 1, column, lp, 3, comm, &req[3]);
		MPI_Irecv(&(mat[0][m]), 1, column, rp, 3, comm, &req[3]);

		MPI_Waitall(4, req, MPI_STATUS_IGNORE);

		odd_even = !step;

		// Now begin ising model
		for(i = 0; i<n; i++){
			for(j = odd_even; j<m; j+=2){
				energy = -1 * mat[i][j] * (mat[i+1][j] + mat[i-1][j] + mat[i][j-1] + mat[i][j+1]);
				flip_energy = mat[i][j] * (mat[i+1][j] + mat[i-1][j] + mat[i][j-1] + mat[i][j+1]);

				d_energy = flip_energy - energy;

				if(d_energy <= 0){
					mat[i][j] *= -1;
				} else if ( drand48() < exp(-1 * d_energy / temp)){
					mat[i][j] *= -1;
				}
			}

			odd_even = !odd_even;

		}

		// Once decorrelation has occured, measure magnetization
		for(i=0; i< n; i++){
			for(j=0; j<m; j++){
				measurements[iter] += mat[i][j];
			}
		}
	}

	delete[] req;

}

// Function which print the matrix on a certain rank
void print_rank(double** mat, int n, int m, int rank, int p_rank, MPI_Comm comm){
	int i, j;
	MPI_Barrier(comm);


	if(rank == p_rank){
		printf("\n-----------Printing matrix from rank %d------------\n", rank);
		for(i=-1; i<n+1; i++){
			for(j=-1; j<m+1; j++){
				printf(" %f", mat[i][j]);
			}
			printf("\n");
		}
		printf("---------------------------------------------------\n");
	}
	MPI_Barrier(comm);
}
