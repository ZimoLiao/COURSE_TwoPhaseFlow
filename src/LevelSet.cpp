#include"LevelSet.h"

#include<fstream>
#include<iomanip>

using namespace std;

void WriteData(int nx, int ny, double** ux, double** uy, double** phi)
{
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

			fout << left << setw(8) << i << ' '\
				<< left << setw(8) << j << ' '\
				<< scientific << setprecision(12) << ux[i][j] << ' '\
				<< scientific << setprecision(12) << uy[i][j] << ' '\
				<< scientific << setprecision(12) << phi[i][j] \
				<< endl;

		}
	}

	fout.close();
}
