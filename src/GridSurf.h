#ifndef GRIDSURF_H_
#define GRIDSURF_H_

#include"Lib.h"

class GridSurf
{
private:
	// parameters
	int nx_, ny_;

	// data
	double* data_;

public:
	GridSurf();
	GridSurf(size_t nx, size_t ny);
	~GridSurf();

	void Init(size_t nx, size_t ny);

	double& w(int i, int j);
	double& e(int i, int j);
	double& s(int i, int j);
	double& n(int i, int j);
};

#endif