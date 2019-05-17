#pragma once
#include "head.h"
#include "grid.h"
#include "solver.h"


class FEM: public GRID, public SOLVER {
public:
	void init(const function2D &_u, const function2D &_f, const function1D &_lambda, double _sigma, bool _isGridUniform, bool _isTimeUniform, int _condType, int _coefGrid, int _coefTime);
	void solve(std::ofstream& fout);
	


protected:
	double lambda0, lambda1;
	double sigma;
	double t;
	function2D f, u;
	function1D lambda;
	double calcNormAtMainNodes(const vector1D &x) {
		double tmp = 0;
		int count = 0;
		for (size_t i = 0; i < x.size(); i++)
			if (i % int(pow(2, coefGrid)) == 0) {
				tmp += x[i] * x[i];
				count++;
			}
		cout << count << endl;
		return sqrt(tmp);
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
