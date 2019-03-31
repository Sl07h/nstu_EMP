#pragma once
#include "head.h"
#include "slae.h"
#include "grid.h"



// Класс МКР. L-образная область.
class FDM : public GRID, public SLAE
{
public:
	void init(function2D &_u, function2D &_f, int _testNumber, bool isGridUniform);
	void inputEquationParameters();
	void outputSLAE(const string &fileA, const string &fileB);
	void transformGridToSLAE();
	void calcAbsResidual(const string &filename);

private:
	double lambda = 1, gamma = 1;					// коэффициенты диффуров
	double C1 = 1, C2 = 1;
	int testNumber;

	function2D u, f;


	double calcFirstDerivativeX(double x, double y) {
		const double h = 0.0001;
		return (-f(x + 2 * h, y) + 8 * f(x + h, y) - 8 * f(x - h, y) + f(x - 2 * h, y)) / (12 * h);
	}


	double calcFirstDerivativeY(double x, double y) {
		const double h = 0.0001;
		return (-f(x, y + 2 * h) + 8 * f(x, y + h) - 8 * f(x, y - h) + f(x, y - 2 * h)) / (12 * h);
	}
};

