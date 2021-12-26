#include<iostream>
#include<fstream>
#include<iomanip>

#include"LevelSet.h"

using namespace std;


int main()
{
	const double pi = 3.141592653589793;

	/* pre-processing */
	// grid initialization
	const int nx = 100, ny = 100;
	double ux[nx][ny];
	double uy[nx][ny];
	double phi[nx][ny];

	// velocity field setting
	// "deformation of a circle in shear"
	double x, y;
	for (int i = 0; i != nx; i++) {
		for (int j = 0; j != ny; j++) {
			x = double(i + 0.5) / double(nx);
			y = double(j + 0.5) / double(ny);
			ux[i][j] = -2 * pi * cos((x - 0.5) * pi) * sin((y - 0.5) * pi);
			uy[i][j] = 2 * pi * sin((x - 0.5) * pi) * cos((y - 0.5) * pi);

			// TEST
			phi[i][j] = 0.0;
		}
	}

	// sign distance function initialization
	double xc = 0.5, yc = 0.3, rc = 0.2;
	for (int i = 0; i != nx; i++) {
		for (int j = 0; j != ny; j++) {
			x = double(i + 0.5) / double(nx);
			y = double(j + 0.5) / double(ny);

			phi[i][j] = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc)) - rc;
		}
	}



	// TEST
	// just advancing! without reinitialization
	double dt = 0.01;

	for (int t = 0; t != 100; t++) {


		double uw, ue, us, un, phin;
		for (int i = 0; i != nx; i++) {
			for (int j = 0; j != ny; j++) {
				// calculate exact velocity at boundaries
				x = double(i) / double(nx);
				y = double(j + 0.5) / double(ny);
				uw = -2 * pi * cos((x - 0.5) * pi) * sin((y - 0.5) * pi);

				x = double(i + 1.0) / double(nx);
				y = double(j + 0.5) / double(ny);
				ue = -2 * pi * cos((x - 0.5) * pi) * sin((y - 0.5) * pi);

				x = double(i + 0.5) / double(nx);
				y = double(j) / double(ny);
				us = 2 * pi * sin((x - 0.5) * pi) * cos((y - 0.5) * pi);

				x = double(i + 0.5) / double(nx);
				y = double(j + 1.0) / double(ny);
				un = 2 * pi * sin((x - 0.5) * pi) * cos((y - 0.5) * pi);


				// QUICK


			}
		}

	}



	/* post-processing */

	// TEST
	ofstream fout;
	fout.open("output.dat");

	// header
	fout << "TITLE     = \"Flow\"" << endl;
	fout << "FILETYPE  = FULL" << endl;
	fout << "VARIABLES = \"x\", \"y\", \"u\", \"v\", \"phi\"" << endl;
	fout << "ZONE    F = point" << endl;
	fout << "        I = " << nx << endl;
	fout << "        J = " << ny << endl;

	// write data
	for (int i = 0; i != nx; i++) {
		for (int j = 0; j != ny; j++) {

			fout << left << setw(8) << double(i + 0.5) / double(nx) << ' '\
				<< left << setw(8) << double(j + 0.5) / double(ny) << ' '\
				<< scientific << setprecision(12) << ux[i][j] << ' '\
				<< scientific << setprecision(12) << uy[i][j] << ' '\
				<< scientific << setprecision(12) << phi[i][j] \
				<< endl;

		}
	}

	fout.close();

	// deallocation

}