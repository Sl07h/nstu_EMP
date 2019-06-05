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

	pair <int, double> result;
	double lambda = 1;
	double	sigma = 1;
	double omega = 1;
	double hi = 1;

	cout << scientific << setprecision(3);

	vector <function1D> u_s(8), u_c(8), f_s(8), f_c(8);
	u_s[0] = { [](double x) -> double { return 1; } };
	u_s[1] = { [](double x) -> double { return x; } };
	u_s[2] = { [](double x) -> double { return pow(x, 2); } };
	u_s[3] = { [](double x) -> double { return pow(x, 3); } };
	u_s[4] = { [](double x) -> double { return pow(x, 4); } };
	u_s[5] = { [](double x) -> double { return pow(x, 5); } };
	u_s[6] = { [](double x) -> double { return sin(x); } };
	u_s[7] = { [](double x) -> double { return exp(x); } };

	u_c[0] = { [](double x) -> double { return -2; } };
	u_c[1] = { [](double x) -> double { return -2 * x; } };
	u_c[2] = { [](double x) -> double { return -2 * pow(x, 2); } };
	u_c[3] = { [](double x) -> double { return -2 * pow(x, 3); } };
	u_c[4] = { [](double x) -> double { return -2 * pow(x, 4); } };
	u_c[5] = { [](double x) -> double { return -2 * pow(x, 5); } };
	u_c[6] = { [](double x) -> double { return -2 * sin(x); } };
	u_c[7] = { [](double x) -> double { return -2 * exp(x); } };

	vector <string> u_s_names = {
		"$1$",
		"$x$",
		"$x^2$",
		"$x^3$",
		"$x^4$",
		"$x^5$",
		"$sin(x)$",
		"$e^x$",
	};

	ofstream foutTable("report/table.txt");
	foutTable << scientific << setprecision(2);
	foutTable << "a\t$1$\t$x$\t$x^2$\t$x^3$\t$x^4$\t$x^5$\t$sin(x)$\t$e^x$" << endl;
	cout << "Research 1: convergence with diferent u_s and u_c" << endl;
	for (size_t i = 0; i < u_s.size(); i++)
	{
		foutTable << u_s_names[i] << "\t";
		for (size_t j = 0; j < u_c.size(); j++)
		{
			std::cout << int(float(i*u_c.size() + j) * 100.0 / (u_s.size()*u_c.size())) << " %\r";
			f_s[i] = calc_f_s(u_s[i], u_c[i], lambda, sigma, omega, hi);
			f_c[j] = calc_f_c(u_s[i], u_c[i], lambda, sigma, omega, hi);

			FEM fem;
			fem.init(u_s[i], u_c[j], f_s[i], f_c[j], lambda, sigma, omega, hi, true, true, 1, 0, 0);
			fem.inputGrid();
			fem.buildGrid();
			result = fem.solve(1);

			if (j + 1 == u_c.size())
				foutTable << result.second << endl;
			else
				foutTable << result.second << "\t";
		}
	}
	cout << endl;



	// Исследование точности при дроблении сетки
	auto reseacrhConvergence = [&](bool isGridUniform, bool isTimeUniform) {
		for (int i = 0; i < u_s.size(); i++) {
			f_s[i] = calc_f_s(u_s[i], u_c[i], lambda, sigma, omega, hi);
			f_c[i] = calc_f_c(u_s[i], u_c[i], lambda, sigma, omega, hi);

			ofstream fout("report/file_u" + to_string(i) + "." + to_string(isGridUniform) + to_string(isTimeUniform) + ".txt");
			fout << scientific << setprecision(3);
			fout << "i\tnodes\titers\tnorm\n";
			for (int coefGrid = 0; coefGrid < 5; coefGrid++)
			{
				FEM fem;
				fem.init(u_s[i], u_c[i], f_s[i], f_c[i], lambda, sigma, omega, hi, isGridUniform, isTimeUniform, 1, coefGrid, 0);
				fem.inputGrid();
				fem.buildGrid();
				result = fem.solve(1);
				fout << coefGrid << "\t"
					<< fem.getNodesCount() << "\t"
					<< result.first << "\t"
					<< result.second << endl;
			}
			fout.close();
		}
	};


	reseacrhConvergence(true, true);
	cout << "Research 2.1: convergence with grid crushing" << endl;
	reseacrhConvergence(true, false);
	cout << "Research 2.2: convergence with grid crushing" << endl;
	reseacrhConvergence(false, true);
	cout << "Research 2.3: convergence with grid crushing" << endl;
	reseacrhConvergence(false, false);
	cout << "Research 2.4: convergence with grid crushing" << endl;

}