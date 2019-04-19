#pragma once
#include "head.h"
#include "grid.h"
#include "SNE.h"


class FEM: public GRID, public SNE {
public:
	void buildGlobalMatrixA();
	void buildGlobalVectorb();


private:
	double lambda0 = 1, lambda1 = 1;
	double sigma = 1;
	function1D f, u;

	void buildLocalMatrixG(int elemNumber);
	void buildLocalMatrixM(int elemNumber);
	void buildLocalmatrixA(int elemNumber);
	void buildLocalVectorb(int elemNumber);
	vector <vector <double>> GLocal, MLocal, ALocal;
	vector <double> bLocal;

};
