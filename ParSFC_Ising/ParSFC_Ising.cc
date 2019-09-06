#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <chrono>

#include "ParSFC.h"

#include "mpi.h"

using namespace std::chrono;

/*
 * This script used an SFC to calculate the Ising model
 */

// Test function used to print a list
void print_list(std::list<int> a);

// Function which sends the Halo points between processors
void send_edges(SFC* curve, int nprocs, MPI_Datatype* data_send, MPI_Datatype* data_recv, MPI_Comm comm, std::vector<int> s_proc) ;

// Function which follows the ising model
void ParSFC_Ising(SFC* curve, double temp, int nsteps, int rank, int a, int b, int nprocs, MPI_Datatype* data_send, MPI_Datatype* data_recv, MPI_Comm comm, double* &res, std::vector<int> s_proc);

int main(int argc, char* argv[]){

	// Set up initial values
	double temp = 0.05;
	int nsteps = 100;
	int p = 0;
	int c = 0;

	// Size of matrix
	int n = 8;

	// MPI values
	int rank, nprocs;

	// Begin Getopt
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
				n = atoi(argv[optind]);
				break;

			case 'p':
				p = 1;
				break;
			case 'c':
				c = 1;
				break;
		}
	}


	// Test if n is a power of 2
	int test = n;
	while(test % 2 == 0){
		test /= 2;
		if(test == 1){
			break;
		}
	}

	double* res = new double[nsteps]();

	// Begin MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Exit if incorrect values are found
	if(test != 1){
		if(rank == 0){
			printf("This function calculates the Ising model.\n Flags change the follow:\n");
			printf("t - the temperature of the lattice\n");
			printf("n - the number of iterations\n");
			printf("m - the size of the lattice\n");
			printf("p - prints in a way readable by certain functions\n");
			printf("c - returns timings of the setup and function\n");
			printf("\nm value must be a power of  2\n");
		}
		MPI_Finalize();
		return 0;
		
	}
	if(temp < 0){
		if(rank == 0){
		printf("This function calculates the Ising model.\n Flags change the follow:\n");
		printf("t - the temperature of the lattice\n");
		printf("n - the number of iterations\n");
		printf("m - the size of the lattice\n");
		printf("p - prints in a way readable by certain functions\n");
		printf("c - returns timings of the setup and function\n");
		printf("\nt value must be greater that 0\n");
		}
		MPI_Finalize();
		return 0;
	
	}
	if(nsteps < 1){
		if(rank == 0){
			printf("This function calculates the Ising model.\n Flags change the follow:\n");
			printf("t - the temperature of the lattice\n");
			printf("n - the number of iterations\n");
			printf("m - the size of the lattice\n");
			printf("p - prints in a way readable by certain functions\n");
			printf("c - returns timings of the setup and function\n");
			printf("\nn value must be greater that 1\n");
		}
		MPI_Finalize();
		return 0;
	}


	high_resolution_clock::time_point setup_start = high_resolution_clock::now();

	// Create matrix
	double** mat = new double*[n];
	double* array = new double[n*n];

	for(int i=0; i<n; i++){
		mat[i] = &array[i*n];
	}

	// Now use SFC on matrix
	int order = (int) (log(n)/log(2));

	MPI_Barrier(MPI_COMM_WORLD);

	hilbert Hil_SFC(mat, order, n, n);

	// Divide up SFC
	int a, b, rem, div, size;
	div = Hil_SFC.mat_size() / nprocs;
	rem = Hil_SFC.mat_size() % nprocs;
	a = 0;

	// Get relevant part of SFC
	for(int i=0; i< nprocs; i++){
		size = div;
		if(i < rem){
			size++;
		}	
		b = a + size;

		if(rank == i){
			break;
		}
		a = b;
	}

	// Vector which willl hold which points that need to be received
	std::vector< std::list<int> > recv_proc;
	std::vector< std::list<int> > send_proc;

	// Function which finds points that need to be send and what ones need to be received
	SFC_Neighbours(&Hil_SFC, rank, nprocs, a, b, recv_proc, send_proc);	

	// Get total number of neighbours
	int tot_neighbours = 0;
	for(int i=0; i<nprocs; i++){
		tot_neighbours += recv_proc[i].size();
	} 

	// Get total boundary
	int tot_boundary = 0;
	for(int i=0; i<nprocs; i++){
		tot_boundary += send_proc[i].size();
	}

	// Create MPI_Datatype
	// First change vector of lists into a matrix
	// Need to create a contiguous matrix so it can be sent
	int** send_mat = new int*[send_proc.size()];
	int* send_array = new int[tot_boundary]();

	int** recv_mat = new int*[recv_proc.size()];
	int* recv_array = new int[tot_neighbours]();

	int** block_send = new int*[send_proc.size()];
	int* b_send_array = new int[tot_boundary]();

	int** block_recv = new int*[recv_proc.size()];
	int* b_recv_array = new int[tot_neighbours]();

	// Initialise sending arrays	
	int j = 0;
	for(int i=0; i<send_proc.size(); i++){
		send_mat[i] = &send_array[j];
		block_send[i] = &b_send_array[j];
		j +=send_proc[i].size();
	}

	// Initailise receiving arrays
	j = 0;
	for(int i=0; i<recv_proc.size(); i++){
		recv_mat[i] = &recv_array[j];
		block_recv[i] = &b_recv_array[j];
		j +=recv_proc[i].size();
	}


	// Set up values in matrix using iterator through each list
	std::list<int>::iterator it; 

	// Iterate through for receiving arrays
	for(int i=0; i<recv_proc.size(); i++){

		j = 0;
		for(it=recv_proc[i].begin(); it!=recv_proc[i].end(); it++){
			recv_mat[i][j] = *it;
			block_recv[i][j] = 1;
			j++;
		}
	}

	// Iterate through for sending arrays
	for(int i=0; i<send_proc.size(); i++){
		j = 0;
		for(it=send_proc[i].begin(); it!=send_proc[i].end(); it++){
			send_mat[i][j] = *it;
			block_send[i][j] = 1;
			j++;
		}
	}

	// Now Create the Datatype sending and recieving datatypes
	MPI_Datatype* send_type = new MPI_Datatype[nprocs];
	MPI_Datatype* recv_type = new MPI_Datatype[nprocs];

	// Test Recv_mat and send_mat

	for(int i=0; i<nprocs; i++){
		MPI_Type_indexed(recv_proc[i].size(), block_recv[i], recv_mat[i], MPI_DOUBLE, &recv_type[i]);
		MPI_Type_commit(&recv_type[i]);
		MPI_Type_indexed(send_proc[i].size(), block_send[i], send_mat[i], MPI_DOUBLE, &send_type[i]);
		MPI_Type_commit(&send_type[i]);
	}

	// Set up values
	srand48((int)time(NULL) + rank);
	for(int i=a; i<b; i++){
		if(drand48() < .5){
			Hil_SFC.set(1, i);
		} else {
			Hil_SFC.set(-1, i);
		}
	}

	std::vector<int> s_proc;

	for(int i=0; i<send_proc.size(); i++){
		if(send_proc[i].size() != 0){
			s_proc.push_back(i);
		}
	}


	MPI_Barrier(MPI_COMM_WORLD);
	high_resolution_clock::time_point setup_end = high_resolution_clock::now();

	duration<double> setup_time = duration_cast<duration<double>>(setup_end - setup_start);
	

	// Begin Function
	MPI_Barrier(MPI_COMM_WORLD);
	high_resolution_clock::time_point func_start = high_resolution_clock::now();

	ParSFC_Ising(&Hil_SFC, temp, nsteps, rank,  a, b, nprocs, send_type, recv_type, MPI_COMM_WORLD, res, s_proc);

	MPI_Barrier(MPI_COMM_WORLD);
	high_resolution_clock::time_point func_end = high_resolution_clock::now();

	duration<double> func_time = duration_cast<duration<double>>(func_end - func_start);

	// Now gather all the results
	double* global_res = new double[nsteps];

	MPI_Allreduce(res, global_res, nsteps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// Now get the average
	for(int i=0; i<nsteps; i++){
		global_res[i] /= n*n;	
	}

	// Print results
	if(rank == 0){
		if(p == 1){
			printf(" %f %d %d %d ", temp, nsteps, n, nprocs);

			if(c == 1){
				printf(" %f %f ", setup_time.count(), func_time.count());
			}
			
			
			for(int i=0; i<nsteps; i++){
				printf(" %f", global_res[i]);
			}
			printf("\n");
		} else {
			printf("Temperature =  %f, nsteps = %d, size = %d, nprocs = %d\n", temp, nsteps, n, nprocs);
			printf("Final Magnetization  = %f\n", global_res[nsteps - 1]);

			if(c == 1){
				printf("Setup Time = %f, Function Time = %f\n", setup_time.count(), func_time.count());
			}
		}

	}

	// Clean up
	delete[] send_mat;
	delete[] send_array;

	delete[] recv_mat;
	delete[] recv_array;
	
	delete[] block_send;
	delete[] b_send_array;

	delete[] block_recv;
	delete[] b_recv_array;

	delete[] send_type;
	delete[] recv_type;

	delete[] array;
	delete[] mat;

	delete[] res;
	delete[] global_res;

	// Finish MPI
	MPI_Finalize();

	return 0;
}

// Testing function used to print a list
void print_list( std::list<int> a){
	std :: list<int> ::iterator it;
	for(it = a.begin(); it != a.end(); it++){
		std::cout <<' ' <<*it;
	}

	std::cout << '\n';

}

// Send Halo points between each SFC
void send_edges(SFC* curve, int no_proc, MPI_Datatype* data_send, MPI_Datatype* data_recv, MPI_Comm comm, std::vector<int> s_proc){

	MPI_Request* send_req = new MPI_Request[s_proc.size()];
	MPI_Request* recv_req = new MPI_Request[s_proc.size()];

	for(int i=0; i<s_proc.size(); i++){
		MPI_Isend(curve->start(), 1, data_send[s_proc[i]], s_proc[i], 0, comm, &send_req[i]);
	}

	// Now receive
	for(int i=0; i<s_proc.size(); i++){
		MPI_Irecv(curve->start(), 1, data_recv[s_proc[i]], s_proc[i], 0, comm, &recv_req[i]);
	}

	MPI_Waitall(s_proc.size(), send_req, MPI_STATUSES_IGNORE);
	MPI_Waitall(s_proc.size(), recv_req, MPI_STATUSES_IGNORE);

	delete[] send_req;
	delete[] recv_req;

}

/* This function calculates the metropolis-hastings model for the Ising model 
 * Curve, an SFC class which stores the values of a matrix in a one dimensional array
 * Temp, temperature of the latticec
 * nsteps, number of iterations of the model
 * rank, rank of processor which runs the function
 * a, start of section in SFC
 * b, end of section in SFC
 * nprocs, number of processors
 * Send type, array of MPI_datatypes that must be send to adjacent processors 
 * Recv type, array of MPI_datatype that must be received from adjacent processors
 * comm, communicator used to send data between processors
 * res, array which stores the sum of the spins on each point
 *
 * This function iterates though through an SFC and updates values based on their adjacent neighbours
 * After each iteration the sum of all spins are saved into the res array
 */
void ParSFC_Ising(SFC* curve, double temp, int nsteps, int rank, int a, int b, int nprocs, MPI_Datatype* data_send, MPI_Datatype* data_recv, MPI_Comm comm, double* &res, std::vector<int> s_proc){
	
	// Initalise values
	double energy;
	int up, down, left, right, it;
	it = 0;

	int x, y;
	

	int even_start, odd_start;
	even_start = a;

	// Calculate colour of starting site
	if ((a % 2) == 1){
		even_start++;
		odd_start = a;	
	} else {
		odd_start = a+1;
	}
	

	// Set up values
	send_edges(curve, nprocs, data_send, data_recv, comm, s_proc);

	// Now iterate along SFC
	for(int i=0; i< nsteps; i++){
		// Iterate along even points in the SFC
		for(int j=even_start; j<b; j+=2){
			curve->index_2D(j, x, y);
			
			curve->index_1D(up, x, y + 1);
			curve->index_1D(down, x, y - 1);
			curve->index_1D(left, x + 1, y);
			curve->index_1D(right, x - 1, y);

			energy = -1 * (*curve)[j] *( (*curve)[up] + (*curve)[down] + (*curve)[left] + (*curve)[right]);

			if( -2 * energy <= 0){
				curve->set(-1 * (*curve)[j], j);
			} else if(drand48() < exp((2*energy)/temp)){
				curve->set(-1 * (*curve)[j], j);
			}

		}

		// Send edges again
		send_edges(curve, nprocs, data_send, data_recv, comm, s_proc);

		// Iterate along the odd points of the SFC
		for(int j=odd_start; j<b; j+=2){
			curve->index_2D(j, x, y);
			
			curve->index_1D(up, x, y + 1);
			curve->index_1D(down, x, y - 1);
			curve->index_1D(left, x + 1, y);
			curve->index_1D(right, x - 1, y);

			energy = -1 * (*curve)[j] *( (*curve)[up] + (*curve)[down] + (*curve)[left] + (*curve)[right]);
	
			if( -2 * energy <= 0){
				curve->set(-1 * (*curve)[j], j);
			} else if(drand48() < exp((2*energy)/temp)){
				curve->set(-1 * (*curve)[j], j);
			}

		}
		send_edges(curve, nprocs, data_send, data_recv, comm, s_proc);

		// Calculate the magnetization
		for(int j=a; j<b; j++){
			res[i] += (*curve)[j];			
		}


	}

}









