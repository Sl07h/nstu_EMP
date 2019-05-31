#include "fem.h"
#include <thread>

function1D calcFirstDerivative(const function1D& f) {
	return [f](double x) -> double {
		const double h = 0.00001;
		return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
	};
}

function1D calc_f_s(
	const function1D & u_s,
	const function1D & u_c,
	double lambda,
	double sigma,
	double omega,
	double hi
) {
	// f_s = - lambda * div(grad u_s) - omega * sigma * u_c - omega^2 * hi * u_s
	return [=](double x) -> double {
		using namespace placeholders;
		auto divGrad = calcFirstDerivative(calcFirstDerivative(u_s));
		return -lambda * divGrad(x) - omega * sigma * u_c(x) - omega * omega * hi * u_s(x);
	};
}


function1D calc_f_c(
	const function1D & u_s,
	const function1D & u_c,
	double lambda,
	double sigma,
	double omega,
	double hi
) {
	// f_c = - lambda * div(grad u_c) + omega * sigma * u_s - omega^2 * hi * u_c
	return [=](double x) -> double {
		using namespace placeholders;
		auto divGrad = calcFirstDerivative(calcFirstDerivative(u_c));
		return -lambda * divGrad(x) + omega * sigma * u_s(x) - omega * omega * hi * u_c(x);
	};
}


void main() {

	//vector <double>  lambda = { 0, 0.5, 5, 10, 100 };
	int coefGrid = 0;
	bool isGridUniform = true;
	bool isTimeUniform = true;
	double lambda = 1;
	double	sigma = 1;
	double omega = 1;
	double hi = 1;

	//cout << fixed << setprecision(2);

	vector <function1D> u_s(3), u_c(3), f_s(3), f_c(3);
	u_s[0] = { [](double x) -> double { return 3 * x; } };
	u_s[1] = { [](double x) -> double { return 2 * x * x; } };
	u_s[2] = { [](double x) -> double { return 2 * x * x * x; } };

	u_c[0] = { [](double x) -> double { return -10 * x; } };
	u_c[1] = { [](double x) -> double { return x * x; } };
	u_c[2] = { [](double x) -> double { return x * x* x; } };

	for (size_t i = 0; i < 3; i++)
	{
		f_s[i] = calc_f_s(u_s[i], u_c[i], lambda, sigma, omega, hi);
		f_c[i] = calc_f_c(u_s[i], u_c[i], lambda, sigma, omega, hi);
	}


	FEM fem;
	fem.init(u_s[0], u_c[0], f_s[0], f_c[0], lambda, sigma, omega, hi, isGridUniform, isTimeUniform, 1, coefGrid, 0);
	fem.inputGrid();
	fem.buildGrid();
	//fem.inputTime();
	//fem.buildTimeGrid();
	fem.solve();
}