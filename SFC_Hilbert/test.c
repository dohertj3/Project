#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hilbert.h"

// Script to test hilbert functions
int main(){
	int i, j, k;	

	// Number of Rows and Columns
	int n, m;
	n = 8;
	m = 8;

	// Create basic matrix
	double** mat = calloc(n, sizeof(double*));
	double* arr = calloc(n*m, sizeof(double));
	
	// Initialise values
	for(i=0; i<n*m; i++){
		arr[i] = i;
	}

	for(i=0; i<n; i++){
		mat[i] = &arr[i*m];
	}

	// Create hilbert struct
	sfc_hilbert test;
	test.rank = 3;
	test.array = calloc(n*m, sizeof(double));

	// Test Hilbert Initialisation
	hilbert_init(mat, &test);	


	// Print out matrix
	printf("\n\n------------------Matrix----------------------\n\n");

	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			printf(" %f", mat[i][j]);
		}
		printf("\n");
	}

	//Print out curve
	printf("\n\n-------------Hilbert Curve------------------\n\n");

	for(i=0; i<n*m; i++){
		printf(" %f", test.array[i]);
	}
	printf("\n");

	// Free Data
	free(arr);
	free(mat);
	free(test.array);

	return 0;
}
