#pragma once
#include "head.h"
#include "grid.h"
#include "solver.h"


class FEM: public GRID, public SOLVER {
public:
	void init(function1D &_u, function1D &_f, bool _isGridUniform, bool _isTimeUniform, int _condType, int _coefGrid, int _coefTime);
	void solve();
	


protected:
	double lambda0 = 1, lambda1 = 1;
	double sigma = 1;
	function1D f, u;

	void buildGlobalMatrixA();
	void buildGlobalVectorb();
	void printGlobalMatrixA();
	void printGlobalVectorb();

	void buildLocalMatrixG(int elemNumber);
	void buildLocalMatrixM(int elemNumber);
	void buildLocalmatrixA(int elemNumber);
	void buildLocalVectorb(int elemNumber);
	matrix2D GLocal, MLocal, ALocal;
	vector1D bLocal;

};
