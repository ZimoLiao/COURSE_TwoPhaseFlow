#include<iostream>
#include<string>

#include"LevelSetSolver.h"

using namespace std;

int main()
{
	int n = 256, s = 1, t = 3, r = 0;

	// 固定最大cfl数计算
	double umax = 2.0 * pi;
	double cfl = 0.8;
	double dt = cfl * 1.0 / n / umax, tmax = 2.0;
	int step_r = floor(0.1 / dt);

	LevelSetSolver ls(n, dt, tmax, s, t, r, step_r);
	ls.InitVelocity(3);
	ls.InitPhi(2);

	ls.Calculation();

	string fname = "n" + to_string(n) + "s" + to_string(s) + to_string(t) + to_string(r) \
		+ "t" + to_string(int(tmax * 10.0)) + ".dat";
	ls.WriteAll(fname);
}