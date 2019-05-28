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

	vector <function1D> u_s(3), u_c(3), f_s(3), f_c(3);
	u_s[0] = { [](double x) -> double { return 3 * x; } };
	u_s[1] = { [](double x) -> double { return 2 * x * x; } };
	u_s[2] = { [](double x) -> double { return 2 * x * x * x; } };

	u_c[0] = { [](double x) -> double { return x; } };
	u_c[1] = { [](double x) -> double { return x * x; } };
	u_c[2] = { [](double x) -> double { return x * x* x; } };

	for (size_t i = 0; i < 3; i++)
	{
		f_s[i] = calc_f_s(u_s[i], u_c[i], lambda, sigma, omega, hi);
		f_c[i] = calc_f_c(u_s[i], u_c[i], lambda, sigma, omega, hi);
	}

	vector <string> u_names = {
		"$ 3x + t $",
		"$ 2x ^ 2 + t $",
		"$ x ^ 3 + t $",
		"$ x ^ 4 + t $",
		"$ e^x + t $",
		"$ 3x + t $",
		"$ 3x + t ^ 2 $",
		"$ 3x + t ^ 3 $",
		"$ 3x + e ^ t $",
		"$ 3x + sin(t) $",
		"$ e^x + t ^ 2 $",
		"$ e^x + t ^ 3 $",
		"$ e^x + e ^ t $",
		"$ e^x + sin(t) $"
	};






	FEM fem;
	fem.init(u_s[0], u_c[0], f_s[0], f_c[0], lambda, sigma, omega, hi, isGridUniform, isTimeUniform, 1, coefGrid, 0);
	fem.inputGrid();
	fem.buildGrid();
	fem.inputTime();
	fem.buildTimeGrid();
	fem.solve();


	/*
	ofstream foutTable("report/table.txt");
	foutTable << scientific << setprecision(2);
	foutTable << "a\t$1$\t$u+x+1$\t$u^2$\t$u^2+1$\t$u^3$\t$u^4$\t$e^u$\tsinu" << endl;
	cout << "Research 1: convergence with diferent u and lambda" << endl;
	for (size_t i = 0; i < u.size(); i++)
	{
		foutTable << u_names[i] << "\t";
		for (size_t j = 0; j < lambda.size(); j++)
		{
			std::cout << int(float(i*lambda.size() + j) * 100.0 / (u.size()*lambda.size())) << " %\r";
			f[i] = calcRightPart(lambda[j], u[i], sigma);
			FEM fem;
			fem.init(u[i], f[i], lambda[j], sigma, isGridUniform, isTimeUniform, 1, coefGrid, 0);
			fem.inputGrid();
			fem.buildGrid();
			fem.inputTime();
			fem.buildTimeGrid();
			foutTable << fem.solve().second;
			if (j + 1 < lambda.size())
				foutTable << "\t";
		}
		foutTable << endl;
	}
	foutTable.close();
	cout << endl;


	// Исследование точности при дроблении сетки
	auto reseacrhConvergence = [&](bool isGridUniform, bool isTimeUniform) {
		for (int i = 0; i < u.size(); i++) {
			f[i] = calcRightPart(lambda[1], u[i], sigma);

			string prefix = "";
			if (!isGridUniform)
				prefix = "Non";

			ofstream fout("report/file_u" + to_string(i) + "." + to_string(isGridUniform) + to_string(isTimeUniform) + ".txt");
			fout << scientific << fixed;
			fout << "i\tnodes\titers\tnorm\n";
			for (int coefGrid = 0; coefGrid < 5; coefGrid++)
			{
				string gridFile = "grids/" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";
				string gridBorderFile = "grids/Border" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";
				FEM fem;
				fem.init(u[i], f[i], lambda[1], sigma, isGridUniform, isTimeUniform, 1, coefGrid, 0);
				fem.inputGrid();
				fem.buildGrid();
				fem.inputTime();
				fem.buildTimeGrid();
				auto result = fem.solve();
				fout << coefGrid << "\t"
					<< fem.getNodesCount() << "\t"
					<< result.first << "\t"
					<< result.second << endl;
			}
			fout.close();
		}
	};


	thread R21(reseacrhConvergence, true, true);
	thread R22(reseacrhConvergence, true, false);
	thread R23(reseacrhConvergence, false, true);
	thread R24(reseacrhConvergence, false, false);

	R21.join();
	cout << "Research 2.1: convergence with grid crushing" << endl;
	R22.join();
	cout << "Research 2.2: convergence with grid crushing" << endl;
	R23.join();
	cout << "Research 2.3: convergence with grid crushing" << endl;
	R24.join();
	cout << "Research 2.4: convergence with grid crushing" << endl;*/
}