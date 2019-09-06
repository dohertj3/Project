#include "ParSFC.h"


// ------------------------- Base SFC Class -------------------
// Base Class Constructor
SFC::SFC(int n, int m, int t_rank):
	rank(t_rank),
	size(n*m),
	side(n),
	array(new double[n*m]){
}

// Deep Copy Constructor
SFC::SFC(const SFC &SFC_in) :
	rank(SFC_in.rank),
	size(SFC_in.size),
	side(SFC_in.size),
	array(new double[size]() ){
	for(int i=0; i<size; i++){
		array[i] = SFC_in.array[i];	

	}
}

// Destructor
SFC::~SFC(){
	delete[] array;
}

// Reference operator
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

// Start
double* SFC::start(){
	return array;
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

// Helper function
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

// Deep Copy Constructor


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

// Function which prints the Hilbert curve as a matrix
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



// Function which checks which processes need to be talked to and points need to be sent
void SFC_Neighbours(SFC* test, int rank, int nproc, int a, int b, std::vector< std::list<int> > &com_r, std::vector< std::list<int> > &com_s){
	

	// Resize vectors
	com_r.resize(nproc);
	com_s.resize(nproc);

	// Set up variables to check for neighbours
	int n = test->length();

	int div = test->mat_size() / nproc;
	int rem = test->mat_size() % nproc;
	
	// Have each processor run though
	int temp_size = 0;

	double proc = -1;

	// Assign values to the matrix based on what processors they are on
	for(int i=0; i<test->mat_size(); i++){


		if(temp_size == i){
			proc++;
			temp_size += div;
			if(proc < rem ){
				temp_size++;
			}
		} 
		test->set(proc, i);
	}
	
	// Now begin checking neighbours
	int up, down, left, right;
	int x, y;

	// Loop through local section to find what processors are close to each other
	for(int i=a; i < b; i++){
		test->index_2D(i, x, y);

		test->index_1D(up, x, y+1);
		test->index_1D(down, x, y-1);
		test->index_1D(right, x+1, y);
		test->index_1D(left, x-1, y);

		// Used to get Neighbours, or what values will be received
		com_r[(int)(*test)[up]].push_back(up);
		com_r[(int)(*test)[down]].push_back(down);
		com_r[(int)(*test)[right]].push_back(right);
		com_r[(int)(*test)[left]].push_back(left);

		// Used to get own borders, or what values will be sent
		com_s[(int)(*test)[up]].push_back(i);
		com_s[(int)(*test)[down]].push_back(i);
		com_s[(int)(*test)[right]].push_back(i);
		com_s[(int)(*test)[left]].push_back(i);

	}

	// Get rid of points on processor
	com_r[rank].erase(com_r[rank].begin(), com_r[rank].end());
	com_s[rank].erase(com_s[rank].begin(), com_s[rank].end());
	
	// Now group in size order and remove duplicates
	for(int i=0; i<com_r.size(); i++){
		com_r[i].sort();
		com_r[i].unique();

		com_s[i].sort();
		com_s[i].unique();
		
	}


}
