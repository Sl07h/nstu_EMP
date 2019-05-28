#pragma once
#include "head.h"
#include "grid.h"
#include "solver.h"


class FEM : public GRID, public SOLVER {
public:

	void init(
		const function1D & _u_s,
		const function1D & _u_c,
		const function1D & _f_s,
		const function1D & _f_c,
		double _lambda,
		double _sigma,
		double _omega,
		double _hi,
		bool _isGridUniform,
		bool _isTimeUniform,
		int _condType,
		int _coefGrid,
		int _coefTime
	);
	void solve();
	inline int getNodesCount() { return nodesCount; }


protected:
	double p00, p01, c00, c01;
	function1D u_s, u_c, f_s, f_c;
	double lambda, sigma, omega, hi;
	
	double calcNormAtMainNodes(const vector1D &x) {
		double tmp = 0;
		for (size_t i = 0; i < x.size(); i++)
			if (i % 2 == 0)
				tmp += pow((x[i] - u_c(nodes[i].x)), 2);
			else
				tmp += pow((x[i] - u_s(nodes[i].x)), 2);
		return sqrt(tmp) / nodes.size();
	}

	void buildGlobalMatrixA();
	void buildGlobalVectorb();
	/*void printGlobalMatrixA();
	void printGlobalVectorb();*/

	void buildLocalmatrixA(int elemNumber);
	void buildLocalVectorb(int elemNumber);
	matrix2D ALocal;
	vector1D bLocal;
};
