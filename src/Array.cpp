#include "Array.h"

Array::Array()
{
	ni_ = 1;
	nj_ = 1;
	nghost_ = 0;
	sizei_ = 1;
	sizej_ = 1;
	size_ = 1;

	data_ = new double[1];
	data_[0] = 0.0;
}

Array::Array(int ni, int nj, int nghost)
{
	ni_ = ni;
	nj_ = nj;
	nghost_ = nghost;

	sizei_ = ni_ + 2 * nghost_;
	sizej_ = nj_ + 2 * nghost_;
	size_ = sizei_ * sizej_;
	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = 0.0; }
}

Array::Array(int ni, int nj, int nghost, double d)
{
	ni_ = ni;
	nj_ = nj;
	nghost_ = nghost;

	sizei_ = ni_ + 2 * nghost_;
	sizej_ = nj_ + 2 * nghost_;
	size_ = sizei_ * sizej_;
	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = d; }
}

Array::Array(const Array& a)
{
	ni_ = a.ni_;
	nj_ = a.nj_;
	nghost_ = a.nghost_;
	sizei_ = a.sizei_;
	sizej_ = a.sizej_;
	size_ = a.size_;

	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = a.data_[i]; }
}

Array::~Array()
{
	delete[] data_;
}

void Array::Init(int ni, int nj, int nghost)
{
	ni_ = ni;
	nj_ = nj;
	nghost_ = nghost;

	sizei_ = ni_ + 2 * nghost_;
	sizej_ = nj_ + 2 * nghost_;
	size_ = sizei_ * sizej_;
	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = 0.0; }
}

void Array::Init(int ni, int nj, int nghost, double d)
{
	ni_ = ni;
	nj_ = nj;
	nghost_ = nghost;

	sizei_ = ni_ + 2 * nghost_;
	sizej_ = nj_ + 2 * nghost_;
	size_ = sizei_ * sizej_;
	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = d; }
}

void Array::Init(const Array& a)
{
	ni_ = a.ni_;
	nj_ = a.nj_;
	nghost_ = a.nghost_;
	sizei_ = a.sizei_;
	sizej_ = a.sizej_;
	size_ = a.size_;

	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = a.data_[i]; }
}

void Array::ExtrapGhost()
{
	int ic, jc;
	int iw, ie, js, jn;

	// linear extrapolation
	// corners
	for (int ig = 1; ig <= nghost_; ig++) {
		ic = nghost_ - ig;
		jc = nghost_ - ig;
		data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic + 1) + (jc + 1)] \
			- data_[sizej_ * (ic + 2) + (jc + 2)];

		ic = nghost_ - ig;
		jc = nj_ + nghost_ - 1 + ig;
		data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic + 1) + (jc - 1)] \
			- data_[sizej_ * (ic + 2) + (jc - 2)];

		ic = ni_ + nghost_ - 1 + ig;
		jc = nghost_ - ig;
		data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic - 1) + (jc + 1)] \
			- data_[sizej_ * (ic - 2) + (jc + 2)];

		ic = ni_ + nghost_ - 1 + ig;
		jc = nj_ + nghost_ - 1 + ig;
		data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic - 1) + (jc - 1)] \
			- data_[sizej_ * (ic - 2) + (jc - 2)];
	}

	// sides
	for (int ig = 1; ig <= nghost_; ig++) {
		for (int j = nghost_ - ig + 1; j != nj_ + nghost_ + ig - 1; j++) {
			iw = nghost_ - ig;
			ie = ni_ + nghost_ - 1 + ig;
			data_[sizej_ * iw + j] = 2.0 * data_[sizej_ * (iw + 1) + j] \
				- data_[sizej_ * (iw + 2) + j];
			data_[sizej_ * ie + j] = 2.0 * data_[sizej_ * (ie - 1) + j] \
				- data_[sizej_ * (ie - 2) + j];
		}
		for (int i = nghost_ - ig + 1; i != ni_ + nghost_ + ig - 1; i++) {
			js = nghost_ - ig;
			jn = nj_ + nghost_ - 1 + ig;
			data_[sizej_ * i + js] = 2.0 * data_[sizej_ * i + (js + 1)] \
				- data_[sizej_ * i + (js + 2)];
			data_[sizej_ * i + jn] = 2.0 * data_[sizej_ * i + (jn - 1)] \
				- data_[sizej_ * i + (jn - 2)];
		}
	}
}

void Array::ExtrapGhost(int order)
{
	int ic, jc;
	int iw, ie, js, jn;

	switch (order)
	{
	case 1: // linear extrapolation
		// corners
		for (int ig = 1; ig <= nghost_; ig++) {
			ic = nghost_ - ig;
			jc = nghost_ - ig;
			data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic + 1) + (jc + 1)] \
				- data_[sizej_ * (ic + 2) + (jc + 2)];

			ic = nghost_ - ig;
			jc = nj_ + nghost_ - 1 + ig;
			data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic + 1) + (jc - 1)] \
				- data_[sizej_ * (ic + 2) + (jc - 2)];

			ic = ni_ + nghost_ - 1 + ig;
			jc = nghost_ - ig;
			data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic - 1) + (jc + 1)] \
				- data_[sizej_ * (ic - 2) + (jc + 2)];

			ic = ni_ + nghost_ - 1 + ig;
			jc = nj_ + nghost_ - 1 + ig;
			data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic - 1) + (jc - 1)] \
				- data_[sizej_ * (ic - 2) + (jc - 2)];
		}

		// sides
		for (int ig = 1; ig <= nghost_; ig++) {
			for (int j = nghost_ - ig + 1; j != nj_ + nghost_ + ig - 1; j++) {
				iw = nghost_ - ig;
				ie = ni_ + nghost_ - 1 + ig;
				data_[sizej_ * iw + j] = 2.0 * data_[sizej_ * (iw + 1) + j] \
					- data_[sizej_ * (iw + 2) + j];
				data_[sizej_ * ie + j] = 2.0 * data_[sizej_ * (ie - 1) + j] \
					- data_[sizej_ * (ie - 2) + j];

			}
			for (int i = nghost_ - ig + 1; i != ni_ + nghost_ + ig - 1; i++) {
				js = nghost_ - ig;
				jn = nj_ + nghost_ - 1 + ig;
				data_[sizej_ * i + js] = 2.0 * data_[sizej_ * i + (js + 1)] \
					- data_[sizej_ * i + (js + 2)];
				data_[sizej_ * i + jn] = 2.0 * data_[sizej_ * i + (jn - 1)] \
					- data_[sizej_ * i + (jn - 2)];
			}
		}
		break;
	case 2:

		break;

	default: // linear extrapolation
		// corners
		for (int ig = 1; ig <= nghost_; ig++) {
			ic = nghost_ - ig;
			jc = nghost_ - ig;
			data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic + 1) + (jc + 1)] \
				- data_[sizej_ * (ic + 2) + (jc + 2)];

			ic = nghost_ - ig;
			jc = nj_ + nghost_ - 1 + ig;
			data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic + 1) + (jc - 1)] \
				- data_[sizej_ * (ic + 2) + (jc - 2)];

			ic = ni_ + nghost_ - 1 + ig;
			jc = nghost_ - ig;
			data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic - 1) + (jc + 1)] \
				- data_[sizej_ * (ic - 2) + (jc + 2)];

			ic = ni_ + nghost_ - 1 + ig;
			jc = nj_ + nghost_ - 1 + ig;
			data_[sizej_ * ic + jc] = 2.0 * data_[sizej_ * (ic - 1) + (jc - 1)] \
				- data_[sizej_ * (ic - 2) + (jc - 2)];
		}

		// sides
		for (int ig = 1; ig <= nghost_; ig++) {
			for (int j = nghost_ - ig + 1; j != nj_ + nghost_ + ig - 1; j++) {
				iw = nghost_ - ig;
				ie = ni_ + nghost_ - 1 + ig;
				data_[sizej_ * iw + j] = 2.0 * data_[sizej_ * (iw + 1) + j] \
					- data_[sizej_ * (iw + 2) + j];
				data_[sizej_ * ie + j] = 2.0 * data_[sizej_ * (ie - 1) + j] \
					- data_[sizej_ * (ie - 2) + j];

			}
			for (int i = nghost_ - ig + 1; i != ni_ + nghost_ + ig - 1; i++) {
				js = nghost_ - ig;
				jn = nj_ + nghost_ - 1 + ig;
				data_[sizej_ * i + js] = 2.0 * data_[sizej_ * i + (js + 1)] \
					- data_[sizej_ * i + (js + 2)];
				data_[sizej_ * i + jn] = 2.0 * data_[sizej_ * i + (jn - 1)] \
					- data_[sizej_ * i + (jn - 2)];
			}
		}
		break;
	}
}

double& Array::operator()(int i, int j)
{
	// TODO: 没有越界检查
	i += nghost_;
	j += nghost_;

	return data_[size_t(sizej_) * i + j];
}

void Array::operator=(const Array& a)
{
	for (int i = 0; i != size_; i++) { data_[i] = a.data_[i]; }
}

void Array::operator+=(const Array& a)
{
	for (int i = 0; i != size_; i++) { data_[i] += a.data_[i]; }
}

void Array::operator+=(const double& d)
{
	for (int i = 0; i != size_; i++) { data_[i] += d; }
}

void Array::operator-=(const Array& a)
{
	for (int i = 0; i != size_; i++) { data_[i] -= a.data_[i]; }
}

void Array::operator-=(const double& d)
{
	for (int i = 0; i != size_; i++) { data_[i] -= d; }
}

void Array::operator*=(const Array& a)
{
	for (int i = 0; i != size_; i++) { data_[i] *= a.data_[i]; }
}

void Array::operator*=(const double& d)
{
	for (int i = 0; i != size_; i++) { data_[i] *= d; }
}

void Array::operator/=(const Array& a)
{
	for (int i = 0; i != size_; i++) { data_[i] /= a.data_[i]; }
}

void Array::operator/=(const double& d)
{
	for (int i = 0; i != size_; i++) { data_[i] /= d; }
}

Array Array::operator+(const Array& a)
{
	Array an(ni_, nj_, nghost_);
	for (int i = 0; i != size_; i++) { an.data_[i] = this->data_[i] + a.data_[i]; }
	return an;
}

Array Array::operator+(const double& d)
{
	Array an(ni_, nj_, nghost_);
	for (int i = 0; i != size_; i++) { an.data_[i] = this->data_[i] + d; }
	return an;
}

Array Array::operator-(const Array& a)
{
	Array an(ni_, nj_, nghost_);
	for (int i = 0; i != size_; i++) { an.data_[i] = this->data_[i] - a.data_[i]; }
	return an;
}

Array Array::operator-(const double& d)
{
	Array an(ni_, nj_, nghost_);
	for (int i = 0; i != size_; i++) { an.data_[i] = this->data_[i] - d; }
	return an;
}

Array Array::operator*(const Array& a)
{
	Array an(ni_, nj_, nghost_);
	for (int i = 0; i != size_; i++) { an.data_[i] = this->data_[i] * a.data_[i]; }
	return an;
}

Array Array::operator*(const double& d)
{
	Array an(ni_, nj_, nghost_);
	for (int i = 0; i != size_; i++) { an.data_[i] = this->data_[i] * d; }
	return an;
}

Array Array::operator/(const Array& a)
{
	Array an(ni_, nj_, nghost_);
	for (int i = 0; i != size_; i++) { an.data_[i] = this->data_[i] / a.data_[i]; }
	return an;
}

Array Array::operator/(const double& d)
{
	Array an(ni_, nj_, nghost_);
	for (int i = 0; i != size_; i++) { an.data_[i] = this->data_[i] / d; }
	return an;
}
