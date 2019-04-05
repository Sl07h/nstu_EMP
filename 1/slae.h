#pragma once
#include "head.h"

class SLAE {

public:
	void inputSLAEParameters();
	void generateInitualGuess() { x.clear(); x.resize(n, 0); }
	void writeXToFile(const string &filename);
	void writebToFile(const string &filename);

	void convMatrixToDense();
	void writeDenseMatrixToFile(const string &filename);
	void writeSecondDenseMatrixToFile(const string &filename);
	double multLine(vector<double>& line, int i, int mode);
	void mult();

	void Jacobi(double w);
	void GaussSeildel(double w);
	int calcIterative(int mode, double w);
	double findOptimalW(int mode);


protected:
	vector <vector <double>> A, A2;
	vector <double> di, au1, au2, al1, al2;
	vector <double> x, b;
	int n, m, maxiter;
	double E;


	double calcNormE(vector <double> &x);
	double calcRelativeDiscrepancy();
	double calcAii(int i);
};