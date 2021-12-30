#ifndef LEVELSET_H_
#define LEVELSET_H_

#include"Grid.h"
#include"GridSurf.h"

class LevelSet
{
private:
	// parameters
	int nx_, ny_, nghost_;
	double h_;

	// data
	Grid phi_, phin_;
	GridSurf u_;

public:
	LevelSet(int nx, int ny, int nghost, double h);
	~LevelSet();

	void Advancing();
};

#endif