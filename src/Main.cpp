#include<iostream>
#include"Lib.h"
#include"Grid.h"
#include"GridSurf.h"

//using namespace std;


int main()
{
	/*
	const int nx = 10, ny = 10;
	double lx = 1.0, ly = 1.0;
	double h = lx / double(nx); // TODO: 改成两方向适用的

	double xc = 0.5, yc = 0.3, rc = 0.2;

	Grid phi(nx, ny, 2);

	// initialization
	double x, y;
	for (int i = -2; i != nx + 2; i++) {
		for (int j = -2; j != ny + 2; j++) {
			x = (i + 0.5) * h;
			y = (j + 0.5) * h;

			if ((i >= 0 && i < nx) || (j >= 0 && j < ny)) {
				phi(i, j) = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc)) - rc;
			}
		}
	}

	// advancing
	double dt = 0.01;

	double ue, uw, vs, vn;
	for (int n = 0; n != 10; n++) {

		// flow domain
		double dphi;
		for (int i = 0; i != nx; i++) {
			for (int j = 0; j != ny; j++) {


				dphi = 0.0;

				x = (i + 0.5) * h;
				y = (j + 0.5) * h;

				uw = u(x - 0.5 * h, y);
				ue = u(x + 0.5 * h, y);
				vs = u(x, y - 0.5 * h);
				vn = u(x, y + 0.5 * h);

				if (abs(uw) > 1e-10) {
					if (uw > 0.0) {
						dphi += (3.0 * phi(i, j) + 6.0 * phi(i - 1, j) - phi(i - 2, j)) * uw;
					}
					else {
						dphi += (3.0 * phi(i - 1, j) + 6.0 * phi(i, j) - phi(i + 1, j)) * uw;
					}
				}

				if (abs(ue) > 1e-10) {
					if (ue > 0.0) {
						dphi += (3.0 * phi(i + 1, j) + 6.0 * phi(i, j) - phi(i - 1, j)) * ue;
					}
					else {
						dphi += (3.0 * phi(i, j) + 6.0 * phi(i + 1, j) - phi(i + 2, j)) * ue;
					}
				}

				if (abs(vs) > 1e-10) {
					if (vs > 0.0) {
						dphi += (3.0 * phi(i, j) + 6.0 * phi(i, j - 1) - phi(i, j - 2)) * vs;
					}
					else {
						dphi += (3.0 * phi(i, j - 1) + 6.0 * phi(i, j) - phi(i, j + 1)) * vs;
					}
				}

				if (abs(vn) > 1e-10) {
					if (vn > 0.0) {
						dphi += (3.0 * phi(i, j + 1) + 6.0 * phi(i, j) - phi(i, j - 1)) * vn;
					}
					else {
						dphi += (3.0 * phi(i, j) + 6.0 * phi(i, j + 1) - phi(i, j + 2)) * vn;
					}
				}

				dphi *= -dt / 8.0;

				phi(i, j) += dphi;
			}
		}

		// ghost layers

	}
	*/

	GridSurf u;

	cout << u.w(0, 0) << endl;
	cout << u.s(0, 0) << endl;
	cout << u.n(0, 0) << endl;
	cout << u.e(0, 0) << endl;
	
}
