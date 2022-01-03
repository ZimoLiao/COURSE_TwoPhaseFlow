#include "ArrayStaggered.h"

ArrayStaggered::ArrayStaggered()
{
	ni_ = 1;
	nj_ = 1;
	sizei_ = 2;
	sizej_ = 2;

	datai_.Init(2, 1, 0);
	datai_.Init(1, 2, 0);
}

ArrayStaggered::ArrayStaggered(int ni, int nj)
{
	ni_ = ni;
	nj_ = nj;
	sizei_ = (ni_ + 1) * nj_;
	sizej_ = ni_ * (nj_ + 1);

	datai_.Init(ni_ + 1, nj, 0);
	dataj_.Init(ni_, nj + 1, 0);
}

ArrayStaggered::ArrayStaggered(int ni, int nj, double d)
{
	ni_ = ni;
	nj_ = nj;
	sizei_ = (ni_ + 1) * nj_;
	sizej_ = ni_ * (nj_ + 1);

	datai_.Init(ni_ + 1, nj, 0, d);
	dataj_.Init(ni_, nj + 1, 0, d);
}

ArrayStaggered::ArrayStaggered(const ArrayStaggered& as)
{
	ni_ = as.ni_;
	nj_ = as.nj_;
	sizei_ = as.sizei_;
	sizej_ = as.sizej_;

	datai_.Init(as.datai_);
	dataj_.Init(as.dataj_);
}

ArrayStaggered::~ArrayStaggered()
{
}

void ArrayStaggered::Init(int ni, int nj)
{
	ni_ = ni;
	nj_ = nj;
	sizei_ = (ni_ + 1) * nj_;
	sizej_ = ni_ * (nj_ + 1);

	datai_.Init(ni_ + 1, nj, 0);
	dataj_.Init(ni_, nj + 1, 0);
}

void ArrayStaggered::Init(int ni, int nj, double d)
{
	ni_ = ni;
	nj_ = nj;
	sizei_ = (ni_ + 1) * nj_;
	sizej_ = ni_ * (nj_ + 1);

	datai_.Init(ni_ + 1, nj, 0, d);
	dataj_.Init(ni_, nj + 1, 0, d);
}

void ArrayStaggered::Init(const ArrayStaggered& as)
{
	ni_ = as.ni_;
	nj_ = as.nj_;
	sizei_ = as.sizei_;
	sizej_ = as.sizej_;

	datai_.Init(as.datai_);
	dataj_.Init(as.dataj_);
}

double ArrayStaggered::x(int i, int j)
{
	return 0.5 * (datai_(i, j) + datai_(i + 1, j));
}

double ArrayStaggered::y(int i, int j)
{
	return 0.5 * (dataj_(i, j) + dataj_(i, j + 1));
}

double& ArrayStaggered::w(int i, int j)
{
	return datai_(i, j);
}

double& ArrayStaggered::e(int i, int j)
{
	return datai_(i + 1, j);
}

double& ArrayStaggered::s(int i, int j)
{
	return dataj_(i, j);
}

double& ArrayStaggered::n(int i, int j)
{
	return dataj_(i, j + 1);
}
