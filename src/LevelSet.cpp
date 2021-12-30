#include "LevelSet.h"

LevelSet::LevelSet(int nx, int ny, int nghost, double h)
{
	nx_ = nx;
	ny_ = ny;
	nghost_ = nghost;
	h_ = h;

	phi_.Init(nx, ny, nghost);
	phin_.Init(nx, ny, nghost);
}

LevelSet::~LevelSet()
{
}

void LevelSet::Advancing()
{

}
