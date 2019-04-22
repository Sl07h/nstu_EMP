#pragma once
#include "head.h"

class SOLVER
{
public:
	void calcWithLUDecomposition();
	bool shouldCalc(int i);
	void LUdecomposition();
	void executeDirectTraversal();
	void executeReverseTraversal();
	void testSLAE();

protected:
	matrix2D A;
	vector1D di, al, au;
	vector1D b;
	vector1D q, qPrevTime;
	int maxiter;
	double E, delta;

	double calcNormE(const vector1D &x) { return sqrt(x*x); }
	vector1D multAonQ();
};

