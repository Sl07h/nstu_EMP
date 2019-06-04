#include "fem.h"
#include <thread>

function1D calcFirstDerivative(const function1D& f) {
	return [f](double x) -> double {
		const double h = 0.001;
		return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
	};
}


function1D calcSecondDerivative(const function1D& f) {
	return [f](double x) -> double {
		const double h = 0.001;
		return (-f(x + 2 * h) + 16 * f(x + h) - 30 * f(x) + 16 * f(x - h) - f(x - 2 * h)) / (12 * h*h);
	};
}


function3D calc_f(
	const function3D& u,
	double lambda,
	double gamma,
	double sigma,
	double chi
) {
	return [=](double x, double y, double t) -> double
	{
		using namespace placeholders;
		auto u_t = calcFirstDerivative(bind(u, x, y, _1));

		auto u_xx = calcSecondDerivative(bind(u, _1, y, t));
		auto u_yy = calcSecondDerivative(bind(u, x, _1, t));
		auto u_tt = calcSecondDerivative(bind(u, x, y, _1));

		return -lambda * (u_xx(x) + u_yy(y)) + gamma * u(x, y, t) + sigma * u_t(t) + chi * u_tt(t);
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
	double chi = 1;

	cout << fixed << scientific;
	//cout << fixed << setprecision(4);

	vector <function3D> u(3), f(3);
	u[0] = { [](double x, double y, double t) -> double { return pow(x,2) + pow(y,2) + pow(t,3); } };
	u[1] = { [](double x, double y, double t) -> double { return 2 * x * x; } };
	u[2] = { [](double x, double y, double t) -> double { return 2 * x * x * x; } };

	for (size_t i = 0; i < 3; i++)
		f[i] = calc_f(u[i], lambda, gamma, sigma, chi);

	string prefix = "";
	if (!isGridUniform)
		prefix = "Non";
	string gridFile = "grids/" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";
	string gridBorderFile = "grids/Border" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";


	int i = 0;
	FEM fem;
	fem.init(u[i], f[i], lambda, gamma, sigma, chi, isGridUniform, isTimeUniform, coefGrid, coefTime);
	fem.inputGrid();
	fem.buildGrid();
	fem.inputTime();
	fem.buildTimeGrid();
	//fem.solveParabolic();
	fem.solve();


	/*fem.saveGridAndBorder(gridFile, gridBorderFile);
	string runVisualisation = "python plot.py " + gridFile + " " + gridBorderFile;
	system(runVisualisation.c_str());*/
}