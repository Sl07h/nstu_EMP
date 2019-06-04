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
	void solve();
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
			t20, t21, t23,
			t30, t31, t32,
			t10, t12, t13;
	double d1, d2, m1, m2;

	/*double calcNormAtMainNodes(const vector1D &x) {
		double tmp = 0;
		for (size_t i = 0; i < x.size(); i++)
			if (i % 2 == 0)
				tmp += pow((x[i] - u_c(nodes[i].x)), 2);
			else
				tmp += pow((x[i] - u_s(nodes[i].x)), 2);
		return sqrt(tmp) / nodes.size();
	}*/

	void CranckNicolson(int timeLayer);
	void implicitCheme(int timeLayer);

	void buildGlobalMatrixG();
	void buildGlobalMatrixM();
	vector1D buildGlobalVectorb(int timeLayer);
	/*void printGlobalMatrixA();
	void printGlobalVectorb();*/

	void buildLocalVectorb(int elemNumber);
	void buildLocalmatrixG(int elemNumber);
	void buildLocalmatrixM(int elemNumber);

	double f1, f2, f3, f4;
	matrix2D A, G, M;
	matrix2D GLocal, MLocal, tmpMatrix;
	vector1D bLocal, tmpVector;
};