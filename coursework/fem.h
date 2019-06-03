#pragma once
#include "head.h"
#include "grid.h"
#include "solver.h"


class FEM : public GRID, public SOLVER {
public:

	void init(
		function1D _u,
		function1D _f,
		double _lambda,
		double _gamma,
		double _sigma,
		double _chi,
		bool _isGridUniform,
		bool _isTimeUniform,
		int _coefGrid,
		int _coefTime
	);
	void solve();
	inline int getNodesCount() { return nodesCount; }
	void convAToDense();
	void outputA();
	void outputALocal();


protected:
	double	p00, p01, p10, p11,
			c00, c01, c10, c11;
	function1D u, f;
	double lambda, gamma, sigma, chi;

	/*double calcNormAtMainNodes(const vector1D &x) {
		double tmp = 0;
		for (size_t i = 0; i < x.size(); i++)
			if (i % 2 == 0)
				tmp += pow((x[i] - u_c(nodes[i].x)), 2);
			else
				tmp += pow((x[i] - u_s(nodes[i].x)), 2);
		return sqrt(tmp) / nodes.size();
	}*/

	void buildGlobalMatrixA();
	void buildGlobalVectorb();
	/*void printGlobalMatrixA();
	void printGlobalVectorb();*/

	void buildLocalmatrixA();
	void buildLocalVectorb(int elemNumber);
	void buildLocalmatrixG();
	void buildLocalmatrixM();

	double f1, f2, f3, f4;
	matrix2D ALocal, GLocal, MLocal, tmpLocal;
	vector1D bLocal;
};
