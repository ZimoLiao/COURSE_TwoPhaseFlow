#include "Grid.h"

Grid::Grid()
{
	nx_ = 1;
	ny_ = 1;
	nghost_ = 0;

	data_ = new double[size_t(nx_) * ny_];

	ghost_w_ = new double[1];
	ghost_e_ = new double[1];
	ghost_s_ = new double[1];
	ghost_n_ = new double[1];
}

Grid::Grid(size_t nx, size_t ny, size_t nghost)
{
	nx_ = nx;
	ny_ = ny;
	nghost_ = nghost;

	data_ = new double[nx * ny];
	ghost_w_ = new double[nghost * ny];
	ghost_e_ = new double[nghost * ny];
	ghost_s_ = new double[nx * nghost];
	ghost_n_ = new double[nx * nghost];
}

Grid::~Grid()
{
	delete[] data_;
	delete[] ghost_w_;
	delete[] ghost_e_;
	delete[] ghost_s_;
	delete[] ghost_n_;
}

void Grid::Init(size_t nx, size_t ny, size_t nghost)
{
	nx_ = nx;
	ny_ = ny;
	nghost_ = nghost;

	data_ = new double[nx * ny];
	ghost_w_ = new double[nghost * ny];
	ghost_e_ = new double[nghost * ny];
	ghost_s_ = new double[nx * nghost];
	ghost_n_ = new double[nx * nghost];
}

// TODO: Խ����
// ��ô�ĳ�debug��releaseģʽ��һ�����ꣿ
double& Grid::operator()(int i, int j)
{
	bool ic = (i >= 0 && i < nx_);
	bool jc = (j >= 0 && j < ny_);

	if (ic && jc) {
		return data_[size_t(ny_) * i + j];
	}
	else if (ic)
	{
		bool js = (j < 0 && fabs(-j) <= nghost_);
		if (js) {
			return ghost_s_[size_t(nghost_) * i + j + nghost_];
		}
		bool jn = (j >= ny_ && fabs(j - ny_ + 1) <= nghost_);
		if (jn) {
			return ghost_n_[size_t(nghost_) * i + j - ny_];
		}

		cout << "wrong | index out of range\n";
	}
	else if (jc)
	{
		bool iw = (i < 0 && fabs(-i) <= nghost_);
		if (iw) {
			return ghost_w_[size_t(ny_) * (size_t(i) + nghost_) + j];
		}
		bool ie = (i >= nx_ && fabs(i - nx_ + 1) <= nghost_);
		if (ie) {
			return ghost_e_[size_t(ny_) * (size_t(i) - nx_) + j];
		}

		cout << "wrong | index out of range\n";
	}
	else
	{
		cout << "wrong | index out of range\n";
	}
}

void Grid::operator=(const Grid& grid)
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
		cout << "wrong | size don't match\n";
	}
}

Grid Grid::operator+(const Grid& grid)
{
	return Grid();
}

Grid Grid::operator-(const Grid& grid)
{
	return Grid();
}

Grid Grid::operator*(const Grid& grid)
{
	return Grid();
}

Grid Grid::operator/(const Grid& grid)
{
	return Grid();
}
