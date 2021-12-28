#include "Grid.h"

Grid::Grid()
{
	nx_ = 1;
	ny_ = 1;
	nghost_ = 0;

	data_ = new double[nx_*ny_];
}

Grid::Grid(size_t nx, size_t ny, size_t nghost)
{
	nx_ = nx;
	ny_ = ny;
	nghost_ = nghost;

	data_ = new double[nx*ny];
	ghost_w_ = new double[nghost*ny];
	ghost_e_ = new double[nghost*ny];
	ghost_s_ = new double[nx*nghost];
	ghost_n_ = new double[nx*nghost];
}

Grid::~Grid()
{
	delete[] data_;
	delete[] ghost_w_;
	delete[] ghost_e_;
	delete[] ghost_s_;
	delete[] ghost_n_;
}

// TODO: 越界检查
// 怎么改成debug和release模式不一样？宏？
double & Grid::operator()(int i, int j)
{
	bool ic = (i >= 0 && i < nx_);
	bool jc = (j >= 0 && j < ny_);

	if (ic && jc) {
		return data_[size_t(nx_)*i + j];
	}
	else if (ic)
	{
		bool js = (fabs(-j) <= nghost_);
		if (js) {
			return ghost_s_[size_t(nx_)*i + j + nghost_];
		}
		bool jn = (fabs(j - ny_ + 1) <= nghost_);
		if (jn) {
			return ghost_n_[size_t(nx_)*i + j - ny_];
		}

		cout << "wrong | index out of range\n";
	}
	else if (jc)
	{
		bool iw = (fabs(-i) <= nghost_);
		if (iw) {
			return ghost_w_[size_t(nx_)*(i + nghost_) + j];
		}
		bool ie = (fabs(i - nx_ + 1) <= nghost_);
		if (ie) {
			return ghost_e_[size_t(nx_)*(i - nx_) + j];
		}

		cout << "wrong | index out of range\n";
	}
	else
	{
		cout << "wrong | index out of range\n";
	}
}

void Grid::operator=(const Grid & grid)
{
	if (nx_ == grid.nx_ && ny_ == grid.ny_ && nghost_ == grid.nghost_) {
		for (int i = 0; i != nx_ * ny_; i++) {
			data_[i] = grid.data_[i];
		}
		if (nghost_) {
			for (int i = 0; i != nghost_ * ny_; i++) {
				ghost_w_[i] = grid.ghost_w_[i];
				ghost_e_[i] = grid.ghost_e_[i];
			}
			for (int i = 0; i != nx_ * nghost_; i++) {
				ghost_s_[i] = grid.ghost_s_[i];
				ghost_n_[i] = grid.ghost_n_[i];
			}
		}
	}
	else {
		cout << "wrong | size don't match"
	}
}

Grid Grid::operator+(const Grid & grid)
{
	return Grid();
}

Grid Grid::operator-(const Grid & grid)
{
	return Grid();
}

Grid Grid::operator*(const Grid & grid)
{
	return Grid();
}

Grid Grid::operator/(const Grid & grid)
{
	return Grid();
}
