#include "fem.h"
#include <thread>

function1D calcFirstDerivative(const function1D& f) {
	return [f](double x) -> double {
		const double h = 0.00001;
		return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
	};
}

function1D calc_f(
	const function1D& u,
	double lambda,
	double gamma,
	double sigma,
	double chi
) {
	return [=](double x)->double {
		auto divGrad = calcFirstDerivative(calcFirstDerivative(u));
		return -lambda * divGrad(x) + gamma * u(x) + sigma * calcFirstDerivative(u)(x) + chi * calcFirstDerivative(calcFirstDerivative(u))(x);
	};
}



void main() {

	int coefGrid = 0;
	int coefTime = 0;
	bool isGridUniform = true;
	bool isTimeUniform = true;
	double lambda = 1;
	double gamma = 1;
	double	sigma = 1;
	double chi = 1e-11;

	cout << fixed << setprecision(2);

	vector <function1D> u(3), f;
	u[0] = { [](double x) -> double { return 3 * x; } };
	u[1] = { [](double x) -> double { return 2 * x * x; } };
	u[2] = { [](double x) -> double { return 2 * x * x * x; } };

	for (size_t i = 0; i < 3; i++)
		f[i] = calc_f(u[i], lambda, gamma, sigma, chi);

	int i = 0;
	FEM fem;
	fem.init(u[i], f[i], lambda, gamma, sigma, chi, isGridUniform, isTimeUniform, coefGrid, coefTime);
	fem.inputGrid();
	fem.buildGrid();
	//fem.inputTime();
	//fem.buildTimeGrid();
	fem.solve();
}