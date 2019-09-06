#include "sfc.h"


// ------------------------- Base SFC Class -------------------
// Base Class Constructor
SFC::SFC(int n, int m, int t_rank):
	rank(t_rank),
	size(n*m),
	side(n),
	array(new double[n*m]){
}

SFC::~SFC(){
	delete[] array;
}

double SFC::operator[](int i){
	return array[i];
}


// Rank
int SFC::order(){
	return rank;
}

// Side length
int SFC::length(){
	return side;
}

// Size
int SFC::mat_size(){
	return size;
}

// Set value
void SFC::set(double val, int index){
	array[index] = val;
}

// Print 1D SFC
void SFC::print_1D(){
	printf("\n-----------Printing 1D curve-----------\n");
	for(int i=0; i<size; i++){
		printf(" %f", array[i]);
	}
	printf("\n");
}

// ------------------------ Hilbert Derived Class ----------------

void next_rank(double** mat, int t_rank, double* array, int direction, int &i, int &j, int &h);

// Enum used for different directions
enum {
	UP,
	DOWN,
	LEFT,
	RIGHT,
};

// Constructor
hilbert::hilbert(double** mat, int t_rank, int n, int m):
	SFC(n, m, t_rank)
	{
	
	// Create some iterables which will wund through matrix and curve
	int i, j, h;

	i=0; j=0; h=0;

	// Set First value
	array[0] = mat[0][0];

	// Now Call recursive function
	next_rank(mat, rank, array, DOWN, i, j, h);
} 

void next_rank(double** mat, int t_rank, double* array, int direction, int &i, int &j, int &h){
	// If Level is 1
	if(t_rank == 1){
		switch (direction) {
			case LEFT:
				// Move right
				(j)++;
				(h)++;
				array[h] = mat[i][j];
	
				// Move up
				(i)++;
				(h)++;
				array[h] = mat[i][j];
	
				// Move left
				(j)--;
				(h)++;
				array[h] = mat[i][j];
				break;
				
	
			case RIGHT:			
				// Move left
				(j)--;
				(h)++;
				array[h] = mat[i][j];
	
				// Move down
				(i)--;
				(h)++;
				array[h] = mat[i][j];
	
				// Move right
				(j)++;
				(h)++;
				array[h] = mat[i][j];
				break;			
			
			case DOWN:
				// Move up
				(i)++;
				(h)++;
				array[h] = mat[i][j];

				// Move right
				(j)++;
				(h)++;
				array[h] = mat[i][j];
	
				// Move down
				(i)--;
				(h)++;
				array[h] = mat[i][j];
				break;
		
			case UP:
				// Move down
				(i)--;
				(h)++;
				array[h] = mat[i][j];
	
				// Move left
				(j)--;
				(h)++;
				array[h] = mat[i][j];
	
				// Move up
				(i)++;
				(h)++;
				array[h] = mat[i][j];
				break;
		}

	} else {
		switch (direction) {
			case LEFT:

				next_rank(mat, t_rank-1, array, DOWN, i, j, h);

				// Move right
				(j)++;
				(h)++;
				array[h] = mat[i][j];	

				next_rank(mat, t_rank-1, array, LEFT, i, j, h);

				// Move up
				(i)++;
				(h)++;
				array[h] = mat[i][j];

				next_rank(mat, t_rank-1, array, LEFT, i, j, h);
	
				// Move left
				(j)--;
				(h)++;
				array[h] = mat[i][j];

				next_rank(mat, t_rank-1, array, UP, i, j, h);

				break;

			case RIGHT:			
				next_rank(mat, t_rank-1, array, UP, i, j, h);

				// Move left
				(j)--;
				(h)++;
				array[h] = mat[i][j];
	
				next_rank(mat, t_rank-1, array, RIGHT, i, j, h);

				// Move down
				(i)--;
				(h)++;
				array[h] = mat[i][j];
	
	
				next_rank(mat, t_rank-1, array, RIGHT, i, j, h);

				// Move right
				(j)++;
				(h)++;
				array[h] = mat[i][j];


				next_rank(mat, t_rank-1, array, DOWN, i, j, h);

				break;			
	
			case DOWN:
				next_rank(mat, t_rank-1, array, LEFT, i, j, h);

				// Move up
				(i)++;
				(h)++;
				array[h] = mat[i][j];

				next_rank(mat, t_rank-1, array, DOWN, i, j, h);

				// Move right
				(j)++;
				(h)++;
				array[h] = mat[i][j];
	
				next_rank(mat, t_rank-1, array, DOWN, i, j, h);

				// Move down
				(i)--;
				(h)++;
				array[h] = mat[i][j];

				next_rank(mat, t_rank-1, array, RIGHT, i, j, h);

				break;

			case UP:
				next_rank(mat, t_rank-1, array, RIGHT, i, j, h);

				// Move down
				(i)--;
				(h)++;
				array[h] = mat[i][j];
	
				next_rank(mat, t_rank-1, array, UP, i, j, h);

				// Move left
				(j)--;
				(h)++;
				array[h] = mat[i][j];
	
				next_rank(mat, t_rank-1, array, UP, i, j, h);

				// Move up
				(i)++;
				(h)++;
				array[h] = mat[i][j];

				next_rank(mat, t_rank-1, array, LEFT, i, j, h);

				break;
		}
	}
}

// Function which returns the 2D index of the 1D index
// Helper function for index_2D
void rotate( int n, int &x, int &y, int rx, int ry);


void hilbert::index_2D(int index, int &x, int &y){
	int rx, ry, s;
	int t = index;
	x = 0;
	y = 0;

	for( s=1; s<side; s =s*2 ){
		rx = (1 & (t/2) );
		ry = (1 & (t ^ rx) );
		rotate(s, x, y, rx, ry);
		x = x + s * rx;
		y = y + s * ry;
		t = t/4;
	}
}

void rotate( int n, int &x, int &y, int rx, int ry){
	int t;
	if(ry == 0){
		if(rx == 1){
			x = n - 1 - x;
			y = n - 1 - y;
		}
		t = x;
		x = y;
		y = t;	

	}
}

// Function which returns the index in the hilbert curve of the 2D co-ordinate
void hilbert::index_1D(int &d, int x, int y){
	d = 0;
	int rx, ry, s;
	
	for(s = side/2; s>0; s = s/2){
		rx = (x & s) > 0;
		ry = (y & s) > 0;

		d = d + s * s * ((3 * rx) ^ ry);
		rotate(s, x, y, rx, ry);

	}
	
}

void hilbert::print_2D(){
	int d;

	printf("\n\n-----------Printing 2D matrix------------\n");
	for(int i=0; i<side; i++){
		for(int j=0; j<side; j++){
			this->index_1D(d, j, i);

			printf(" %f", (*this)[d]);
		}
		printf("\n");
	}
	printf("---------------Finished----------------------\n\n");

}

// Destructor
hilbert::~hilbert(){
}
