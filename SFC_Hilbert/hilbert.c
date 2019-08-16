#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "hilbert.h"

/*
 * This File Holds all the functions used to create and reference the values in the Hilbert curve 
 *
 */

/*
 * This Function takes a matrix of doubles and updates an array which will store the hilbert ordering of the curve
 *
 * 
 */ 
void next_rank(double** mat, int rank, double* array, int direction, int* i, int* j, int* h);
enum {
	UP,
	DOWN,
	LEFT,
	RIGHT,
};


void hilbert_init(double** mat, sfc_hilbert* sfc){
	// Create iterables which will go through the matrix and the curve;
	int i, j, h;
	
	i=0;
	j=0;
	h=0;

	// Set up first value
	sfc->array[0] = mat[0][0];

	// Now call recursive function
	next_rank(mat, sfc->rank, sfc->array, DOWN, &i, &j, &h);
}

/* 
 * This is a helper function for the hilbert init function;
 * It is a recursive function that calls itself until its at the 1st order
 */
void next_rank(double** mat, int rank, double* array, int direction, int* i, int* j, int* h){
	// If level is 1;
	if(rank == 1){	
		switch (direction) {
			case LEFT:
				// Move right
				(*j)++;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				// Move up
				(*i)++;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				// Move left
				(*j)--;
				(*h)++;
				array[*h] = mat[*i][*j];
				break;
				
	
			case RIGHT:			
				// Move left
				(*j)--;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				// Move down
				(*i)--;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				// Move right
				(*j)++;
				(*h)++;
				array[*h] = mat[*i][*j];
				break;			
			
			case DOWN:
				// Move up
				(*i)++;
				(*h)++;
				array[*h] = mat[*i][*j];

				// Move right
				(*j)++;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				// Move down
				(*i)--;
				(*h)++;
				array[*h] = mat[*i][*j];
				break;
		
			case UP:
				// Move down
				(*i)--;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				// Move left
				(*j)--;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				// Move up
				(*i)++;
				(*h)++;
				array[*h] = mat[*i][*j];
				break;
		}

	} else {
		switch (direction) {
			case LEFT:
				next_rank(mat, rank-1, array, DOWN, i, j, h);

				// Move right
				(*j)++;
				(*h)++;
				array[*h] = mat[*i][*j];	

				next_rank(mat, rank-1, array, LEFT, i, j, h);

				// Move up
				(*i)++;
				(*h)++;
				array[*h] = mat[*i][*j];

				next_rank(mat, rank-1, array, LEFT, i, j, h);
	
				// Move left
				(*j)--;
				(*h)++;
				array[*h] = mat[*i][*j];

				next_rank(mat, rank-1, array, UP, i, j, h);

				break;

			case RIGHT:			
				next_rank(mat, rank-1, array, UP, i, j, h);

				// Move left
				(*j)--;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				next_rank(mat, rank-1, array, RIGHT, i, j, h);

				// Move down
				(*i)++;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				next_rank(mat, rank-1, array, RIGHT, i, j, h);

				// Move right
				(*j)++;
				(*h)++;
				array[*h] = mat[*i][*j];

				next_rank(mat, rank-1, array, DOWN, i, j, h);

				break;			
	
			case DOWN:
				next_rank(mat, rank-1, array, LEFT, i, j, h);

				// Move up
				(*i)++;
				(*h)++;
				array[*h] = mat[*i][*j];

				next_rank(mat, rank-1, array, DOWN, i, j, h);

				// Move right
				(*j)++;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				next_rank(mat, rank-1, array, DOWN, i, j, h);

				// Move down
				(*i)--;
				(*h)++;
				array[*h] = mat[*i][*j];

				next_rank(mat, rank-1, array, RIGHT, i, j, h);

				break;

			case UP:
				next_rank(mat, rank-1, array, RIGHT, i, j, h);

				// Move down
				(*i)--;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				next_rank(mat, rank-1, array, UP, i, j, h);

				// Move left
				(*j)++;
				(*h)++;
				array[*h] = mat[*i][*j];
	
				next_rank(mat, rank-1, array, UP, i, j, h);

				// Move up
				(*i)--;
				(*h)++;
				array[*h] = mat[*i][*j];

				next_rank(mat, rank-1, array, LEFT, i, j, h);

				break;
		}
	}
}
	


/* 
 * This function takes in the x and y values of a point in the matrix and finds there coresponding point in the hilbert curve
 */
/*
 * Create a struct used to find out which shape the section of the hilbert curve we are on is
 */


typedef struct hilbert_shape{
	hilbert_shape** next_shape;
	int* order;
} hilbert_shape;

double hilbert_ref(sfc_hilbert* sfc, int n, int m){
	// Set up hilbert_shapes
	hilbert_shape left, right, up, down;

	// Set up arrays
	left.order = calloc(4, sizeof(int));


}
