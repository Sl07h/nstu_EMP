#pragma once
#include "head.h"
#include "grid.h"
#include "solver.h"


class FEM : public GRID, public SOLVER {
public:
	void init(const function2D &_u, const function2D &_f, const function2D &_lambda, double _sigma, bool _isGridUniform, bool _isTimeUniform, int _condType, int _coefGrid, int _coefTime);
	pair<int, double> solve();
	inline int getNodesCount() { return nodesCount; }


protected:
	double lambda0, lambda1;
	double sigma;
	double t;
	function2D f, u, lambda;
	double calcNormAtMainNodes(const vector1D &x) {
		double tmp = 0;
		for (size_t i = 0; i < x.size(); i++)
			tmp += pow((x[i] - u(nodes[i].x, t)), 2);
		return sqrt(tmp) / nodes.size();
	}

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
