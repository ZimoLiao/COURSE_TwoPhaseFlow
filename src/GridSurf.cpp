#include "GridSurf.h"

GridSurf::GridSurf()
{
	nx_ = 1;
	ny_ = 1;

	data_ = new double[4];

	// TEST
	data_[0] = 0.0;
	data_[1] = 1.0;
	data_[2] = 2.0;
	data_[3] = 3.0;
}

GridSurf::GridSurf(size_t nx, size_t ny)
{
	nx_ = nx;
	ny_ = ny;

	data_ = new double[nx * ny * 2 + nx + ny];
}

GridSurf::~GridSurf()
{
	delete[] data_;
}

void GridSurf::Init(size_t nx, size_t ny)
{
	nx_ = nx;
	ny_ = ny;

	data_ = new double[nx * ny * 2 + nx + ny];

	// velocity initialization
	for (int i = 0; i != nx_; i++) {
		for (int j = 0; j != ny_; j++) {
			

			this->w(i,j)=u()
		}
	}
}

double& GridSurf::w(int i, int j)
{
	return data_[(size_t(nx_) + ny_ + 1) * size_t(i) + j];
}

double& GridSurf::e(int i, int j)
{
	return data_[(size_t(nx_) + ny_ + 1) * (size_t(i) + 1) + j];
}

double& GridSurf::s(int i, int j)
{
	return data_[(size_t(nx_) + ny_ + 1) * size_t(i) + j + ny_];
}

double& GridSurf::n(int i, int j)
{
	return data_[(size_t(nx_) + ny_ + 1) * size_t(i) + j + ny_ + 1];
}
