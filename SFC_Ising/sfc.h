#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * This Defines an SFC base class
 *
 */ 

class SFC{
	protected:
	double* array;
	int rank;
	int size;
	int side;

	public:

	// Constructor
	SFC(int n, int m, int t_rank);

	// Destructor
	virtual ~SFC();

	// Rank
	int order();

	// side length
	int length();

	// Size
	int mat_size();

	// Set
	void set(double val, int index);

	// Print 1D
	void print_1D();

	// Returns value at point in curve
	double operator[](int i);

	// Get 1D index
	virtual void index_1D(int &d, int i, int j)=0;

	// Get 2D index
	virtual void index_2D(int index, int &res_x, int &res_y)=0;

	// Print 2D
	virtual void print_2D()=0;

};

/*
 * This is the Hilbert derived class
 * It has its own constructor and assignment operator
 */

class hilbert : public SFC{
	public:
	
	// Hilbert Constructor
	hilbert(double** mat, int t_rank, int n, int m);

	// Hilbert Destructor
	~hilbert();

	// Get 1D index
	void index_1D(int &d, int i, int j);

	// Get 2D index
	void index_2D(int index, int &res_x, int &res_y);

	// print 2D
	void print_2D();

};
