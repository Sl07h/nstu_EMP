#pragma once
#include "head.h"
#include "slae.h"
#include "grid.h"


// Класс МКР. L-образная область.
class FDM : public GRID, public SLAE
{
public:
	void init(function2D &_u, function2D &_f, bool isGridUniform, int _condType, int _coef);
	void inputEquationParameters();
	void outputSLAE(const string &fileA);
	void transformGridToSLAE();
	double calcAbsResidual(const string &filepath);

	void checkAnswer();



private:
	vector <double> xExp;
	double lambda = 1, gamma = 1;	// коэффициенты диффуров
	double C1 = 5, C2 = 5;
	
	function2D u, f;


	double calcFirstDerivativeX(double x, double y) {
		const double h = 1e-9;
		return (-u(x + 2 * h, y) + 8 * u(x + h, y) - 8 * u(x - h, y) + u(x - 2 * h, y)) / (12 * h);
	}


	double calcFirstDerivativeY(double x, double y) {
		const double h = 1e-9;
		return (-u(x, y + 2 * h) + 8 * u(x, y + h) - 8 * u(x, y - h) + u(x, y - 2 * h)) / (12 * h);
	}
};

