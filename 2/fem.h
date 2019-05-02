#pragma once
#include "head.h"
#include "grid.h"
#include "solver.h"


class FEM: public GRID, public SOLVER {
public:
	void init(const function1D &_u, const function1D &_f, const function2D &_lambda, bool _isGridUniform, bool _isTimeUniform, int _condType, int _coefGrid, int _coefTime);
	void solve();
	


protected:
	double lambda0, lambda1;
	double sigma = 1;
	function1D f, u;
	function2D lambda;
	void buildGlobalMatrixA(double _dt);
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
