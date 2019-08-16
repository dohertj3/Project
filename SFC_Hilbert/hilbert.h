#include "stdlib.h"
#include "stdio.h"
#include "math.h"

/*
 * Struct that holds the information of the hilbert curve
 */
typedef struct{
	int rank;
	double* array;
	int m, n;
} sfc_hilbert;

// Function to  Create Hilbert curve
void hilbert_init(double** mat, sfc_hilbert* sfc);

// Function that references a hilbert curve based on the coordinates of the matrix
double hilbert_ref(sfc_hilbert* sfc, int n, int m);
