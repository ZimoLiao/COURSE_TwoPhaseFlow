#ifndef LIB_H_
#define LIB_H_

constexpr auto pi = 3.14159265358979;

double u(double x, double y)
{
	double u = -2 * pi * cos(pi * (x - 0.5)) * sin(pi * (y - 0.5));
	if (abs(u) > 1e-10) {
		return u;
	}
	else {
		return 0.0;
	}
};

double v(double x, double y)
{
	double v = 2 * pi * sin(pi * (x - 0.5)) * cos(pi * (y - 0.5));
	if (abs(v) > 1e-10) {
		return v;
	}
	else {
		return 0.0;
	}
};


#endif