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


function3D sum_u(
	const function3D& u_x,
	const function3D& u_t,
	double lambda,
	double gamma,
	double sigma,
	double chi
) {
	return [=](double x, double y, double t) -> double
	{
		return (u_x(x, y, t) + u_t(x, y, t));
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

	cout << setprecision(3) << scientific;

	string prefix = "";
	if (!isGridUniform)
		prefix = "Non";
	string gridFile = "grids/" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";
	string gridBorderFile = "grids/Border" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";



	vector <function3D> u_x(8), u_t(8), u(64), f(64);
	u_x[0] = { [](double x, double y, double t) -> double { return 1; } };
	u_x[1] = { [](double x, double y, double t) -> double { return x + y; } };
	u_x[2] = { [](double x, double y, double t) -> double { return pow(x, 2) + pow(y, 2); } };
	u_x[3] = { [](double x, double y, double t) -> double { return pow(x, 3) + pow(y, 3); } };
	u_x[4] = { [](double x, double y, double t) -> double { return pow(x, 4) + pow(y, 4); } };
	u_x[5] = { [](double x, double y, double t) -> double { return pow(x, 5) + pow(y, 5); } };
	u_x[6] = { [](double x, double y, double t) -> double { return sin(x) + sin(y); } };
	u_x[7] = { [](double x, double y, double t) -> double { return exp(x) + exp(y); } };

	u_t[0] = { [](double x, double y, double t) -> double { return 1; } };
	u_t[1] = { [](double x, double y, double t) -> double { return t; } };
	u_t[2] = { [](double x, double y, double t) -> double { return pow(t, 2); } };
	u_t[3] = { [](double x, double y, double t) -> double { return pow(t, 3); } };
	u_t[4] = { [](double x, double y, double t) -> double { return pow(t, 4); } };
	u_t[5] = { [](double x, double y, double t) -> double { return pow(t, 5); } };
	u_t[6] = { [](double x, double y, double t) -> double { return sin(t); } };
	u_t[7] = { [](double x, double y, double t) -> double { return exp(t); } };

	vector <string> u_x_names = {
		"$1$",
		"$x+y$",
		"$x^2+y^2$",
		"$x^3+y^3$",
		"$x^4+y^4$",
		"$x^5+y^5$",
		"$sin(x)+sin(y)$",
		"$e^x+e^y$"
	};


	// Исследование сходимости метода на различных функциях по времени и пространству
	auto reseacrhConvergence = [&](bool isGridUniform, bool isTimeUniform) {
		ofstream table1("report/table_" + to_string(isGridUniform) + to_string(isTimeUniform) + ".txt");
		table1 << setprecision(2) << scientific;
		table1 << "a\t$1$\t$t$\t$t^2$\t$t^3$\t$t^4$\t$t^5$\t$sin(t)$\t$e^t$" << endl;
		for (size_t i = 0; i < 8; i++)
		{
			table1 << u_x_names[i] << "\t";
			for (size_t j = 0; j < 8; j++)
			{
				int k = 8 * i + j;
				u[k] = sum_u(u_x[i], u_t[j], lambda, gamma, sigma, chi);
				f[k] = calc_f(u[k], lambda, gamma, sigma, chi);

				FEM fem;
				fem.init(u[k], f[k], lambda, gamma, sigma, chi, isGridUniform, isTimeUniform, coefGrid, coefTime);
				fem.inputGrid();
				fem.buildGrid();
				fem.inputTime();
				fem.buildTimeGrid();
				if (j != 7)
					table1 << fem.solve() << "\t";
				else
					table1 << fem.solve() << endl;

				/*fem.saveGridAndBorder(gridFile, gridBorderFile);
				string runVisualisation = "python plot.py " + gridFile + " " + gridBorderFile;
				system(runVisualisation.c_str());*/
			}
		}
		table1.close();
	};
	reseacrhConvergence(true, true);
	reseacrhConvergence(true, false);
	reseacrhConvergence(false, true);
	reseacrhConvergence(false, false);

}