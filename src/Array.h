#ifndef ARRAY_H_
#define ARRAY_H_

class Array
{
protected:
	/* parameters */
	int ni_, nj_, nghost_, sizei_, sizej_, size_;

	/* data */
	double* data_;

public:
	/* constructors & destructor */
	Array();
	Array(int ni, int nj, int nghost);
	Array(int ni, int nj, int nghost, double d);
	Array(const Array& a);
	~Array();

	/* initialization methods */
	void Init(int ni, int nj, int nghost);
	void Init(int ni, int nj, int nghost, double d);
	void Init(const Array& a);

	/* extrapolation for ghost layers */
	void ExtrapGhost();
	void ExtrapGhost(int order);

	/* operators */
	// indices
	//	-nghost:-1			0:n+nghost-1	n+nghost:n+2*nghost-1
	//	ghost layer left	data			ghost layer right
	double& operator()(int i, int j);

	void operator=(const Array& a);
	void operator=(const double& d);
	void operator+=(const Array& a);
	void operator+=(const double& d);
	void operator-=(const Array& a);
	void operator-=(const double& d);
	void operator*=(const Array& a);
	void operator*=(const double& d);
	void operator/=(const Array& a);
	void operator/=(const double& d);

	Array operator+(const Array& a);
	Array operator+(const double& d);
	Array operator-(const Array& a);
	Array operator-(const double& d);
	Array operator*(const Array& a);
	Array operator*(const double& d);
	Array operator/(const Array& a);
	Array operator/(const double& d);
};

#endif // !ARRAY_H_