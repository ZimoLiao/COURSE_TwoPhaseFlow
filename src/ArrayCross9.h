#ifndef ARRAYCROSS9_H_
#define ARRAYCROSS9_H_

#include "Array.h"

// 可能实际上不需要
// Array with "cross-9" computational molecule
class ArrayCross9 :
	public Array
{
protected:
	/* index for computational molecule */
	int indp = 0;
	int indw = 0;
	int inde = 0;
	int inds = 0;
	int indn = 0;
	int indww = 0;
	int indee = 0;
	int indss = 0;
	int indnn = 0;

public:
	using Array::Array;

	/* pointers */
	double& P(int i, int j);
	double& W(int i, int j);
	double& E(int i, int j);
	double& S(int i, int j);
	double& N(int i, int j);
	double& WW(int i, int j);
	double& EE(int i, int j);
	double& SS(int i, int j);
	double& NN(int i, int j);

	void Index(int i, int j);
	double& P();
	double& W();
	double& E();
	double& S();
	double& N();
	double& WW();
	double& EE();
	double& SS();
	double& NN();
};


#endif // !ARRAYCROSS9_H_