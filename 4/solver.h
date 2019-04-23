#pragma once
#include "head.h"

class SOLVER
{
public:
	void initSLAE();
	void BiCG();

private:
	vector1D di, al, au, ia, ja;
	vector1D b;
	matrix2D A;

	vector1D x;		// решение на k итерации
	vector1D xPrev;	// решение на k-1 итерации
	vector1D r;		// вектор невязки
	vector1D z;		// вектор спуска
	vector1D p, s;	// вспомогательный вектор
	double alpha, beta;

	int n;			// размерность СЛАУ
	int maxiter;
	double E, delta;

	// Вспомогательные методы класса
	double calcNormE(const vector1D &x) { return sqrt(x*x); }
	void generateInitialGuess();
	vector1D multAOn(const vector1D &v);
	vector1D multAtOn(const vector1D &v);
	bool doStop(int i);
};

