#ifndef ARRAYSTAGGERED_H_
#define ARRAYSTAGGERED_H_

#include "Array.h"

class ArrayStaggered
{
protected:
	/* parameters */
	int ni_, nj_, sizei_, sizej_;

	/* data */
	Array datai_;
	Array dataj_;

public:
	/* constructors & destructor */
	ArrayStaggered();
	ArrayStaggered(int ni, int nj);
	ArrayStaggered(int ni, int nj, double d);
	ArrayStaggered(const ArrayStaggered& as);
	~ArrayStaggered();

	/* initialization methods */
	void Init(int ni, int nj);
	void Init(int ni, int nj, double d);
	void Init(const ArrayStaggered& as);

	/* pointers */
	double x(int i, int j);
	double y(int i, int j);

	double& w(int i, int j);
	double& e(int i, int j);
	double& s(int i, int j);
	double& n(int i, int j);
};

#endif // !ARRAYSTAGGERED_H_