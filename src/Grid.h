#ifndef GRID_H_
#define GRID_H_

#include<iostream>

using namespace std;
using std::size_t;

class Grid
{
private:
	// parameters
	int nx_, ny_, nghost_;

	// data
	double* data_;
	double* ghost_w_;
	double* ghost_e_;
	double* ghost_s_;
	double* ghost_n_;

public:
	Grid();
	Grid(size_t nx, size_t ny, size_t nghost);
	~Grid();

	void Init(size_t nx, size_t ny, size_t nghost);

	double& operator()(int i, int j);
	void operator=(const Grid& grid);
	Grid operator+(const Grid& grid);
	Grid operator-(const Grid& grid);
	Grid operator*(const Grid& grid);
	Grid operator/(const Grid& grid);
};


#endif