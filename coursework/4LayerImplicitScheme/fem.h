#pragma once
#include "head.h"
#include "grid.h"
#include "solver.h"


class FEM : public GRID, public SOLVER {
public:

	void init(
		function3D _u,
		function3D _f,
		double _lambda,
		double _gamma,
		double _sigma,
		double _chi,
		bool _isGridUniform,
		bool _isTimeUniform,
		int _coefGrid,
		int _coefTime
	);
	double solve();
	inline int getNodesCount() { return nodesCount; }
	void convAToDense();
	void outputA();
	void outputG();
	void outputM();


protected:
	double	p00, p01, p10, p11,
			c00, c01, c10, c11;
	function3D u, f;
	double lambda, gamma, sigma, chi;
	double t, dt;
	double t0, t1, t2, t3;
	double  t01, t02, t03,
			t23, t12, t13;
	double d1, d2, m1, m2;

	double calcNormAtMainNodes(const vector1D &x, int time) {
		double normE = 0;
		for (size_t i = 0; i < x.size(); i++)
			normE += pow((x[i] - u(nodes[i].x, nodes[i].y, time)), 2);
		return sqrt(normE) / nodes.size();
	}

	void implicitCheme4(int timeLayer);

	void buildGlobalMatrixG();
	void buildGlobalMatrixM();
	vector1D buildGlobalVectorb(int timeLayer);

	void buildLocalVectorb(int elemNumber);
	void buildLocalmatrixG(int elemNumber);
	void buildLocalmatrixM(int elemNumber);

	double f1, f2, f3, f4;
	matrix2D A, G, M;
	matrix2D GLocal, MLocal, tmpMatrix;
	vector1D bLocal, tmpVector;
};