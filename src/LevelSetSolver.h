#ifndef LEVELSETSOLVER_H_
#define LEVELSETSOLVER_H_

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<math.h>

#include"Array.h"
#include"ArrayStaggered.h"

using namespace std;

constexpr auto pi = 3.1415926535897932;

class LevelSetSolver
{
private:
	/* parameters */
	double dt_, tmax_;
	int n_;
	double l_ = 1.0, h_;
	int step_, stepmax_;

	/* numerical schemes */
	//	spatial discretization
	//	1	QUICK		(3rd-order)
	//	2	HJ-ENO
	int scheme_s_ = 1;
	//	time advancing
	//	1	Euler		(1st-order)
	//	2	TVD-RK2		(2nd-order)
	//	3	TVD-RK3		(3rd-order)
	int scheme_t_ = 3;
	//	reinitialization
	//	0	no reinitialization
	//	1	Simple method
	//	2	Subcell fix method
	int scheme_r_ = 2;
	int step_r_ = 1000;


	/* data */
	Array phin, phin1;
	Array f;
	Array Sigma, D; // for subcell fix reinit
	ArrayStaggered u;


	/* core solver */
	// spatial discretization
	void QUICK(Array& phi);
	void HJENO(Array& phi);

	void F(Array& phi);

	// time advancing
	void EULER();
	void TVDRK2();
	void TVDRK3();

	// reinitialization
	void SimpleReinit();
	void SubcellFix();

	void Reinit();

public:
	/* constructors & destructor */
	LevelSetSolver(int n, double dt, double tmax);
	LevelSetSolver(int n, double dt, double tmax, int scheme_s, int scheme_t, int scheme_r, int step_r);
	~LevelSetSolver();


	/* methods */
	void Init(int n, double dt, double tmax, int scheme_s, int scheme_t, int scheme_r, int step_r);
	void InitVelocity(int type_init);
	void InitPhi(int type_init);

	void Calculation();

	void WriteVelocity();
	void WriteVelocity(string fname);
	void WritePhi();
	void WritePhi(string fname);
	void WriteAll();
	void WriteAll(string fname);
};

#endif // !LEVELSETSOLVER_H_