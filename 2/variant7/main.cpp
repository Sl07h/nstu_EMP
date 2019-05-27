#include "fem.h"
#include <thread>


function1D calcFirstDerivative(const function1D& f) {
	return [f](double x) -> double {
		const double h = 0.00001;
		return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
	};
}


function2D calcRightPart(const function2D& lambda, const function2D& u, double sigma) {
	return [=](double x, double t) -> double {
		using namespace std::placeholders;
		auto duBydt = calcFirstDerivative(std::bind(u, x, _1));
		auto duBydx = calcFirstDerivative(std::bind(u, _1, t));
		auto lambda_grad = [=](double x, double t) -> double {
			//return lambda(u(x, t)) * duBydx(x); // var 5
			return lambda(duBydx(x), x) * duBydx(x); // var 7
		};
		auto div = calcFirstDerivative(std::bind(lambda_grad, _1, t));
		return -div(x) + sigma * duBydt(t);
	};
}



void main() {

	vector <function2D> u(14), f(14);
	u[0] = { [](double x, double t) -> double { return 3 * x + t; } };
	u[1] = { [](double x, double t) -> double { return 2 * x*x + t; } };
	u[2] = { [](double x, double t) -> double { return x * x*x + t; } };
	u[3] = { [](double x, double t) -> double { return x * x*x*x + t; } };
	u[4] = { [](double x, double t) -> double { return exp(x) + t; } };
	u[5] = { [](double x, double t) -> double { return 3 * x + t; } };
	u[6] = { [](double x, double t) -> double { return 3 * x + t * t; } };
	u[7] = { [](double x, double t) -> double { return 3 * x + t * t*t; } };
	u[8] = { [](double x, double t) -> double { return 3 * x + exp(t); } };
	u[9] = { [](double x, double t) -> double { return 3 * x + sin(t); } };
	u[10] = { [](double x, double t) -> double { return exp(x) + t * t; } };
	u[11] = { [](double x, double t) -> double { return exp(x) + t * t*t; } };
	u[12] = { [](double x, double t) -> double { return exp(x) + exp(t); } };
	u[13] = { [](double x, double t) -> double { return exp(x) + sin(t); } };

	vector <function2D>  lambda(8);
	lambda[0] = { [](double u, double x) -> double {return 1; } };
	lambda[1] = { [](double u, double x) -> double {return u + x + 1; } };
	lambda[2] = { [](double u, double x) -> double {return u * u; } };
	lambda[3] = { [](double u, double x) -> double {return u * u + 1; } };
	lambda[4] = { [](double u, double x) -> double {return u * u*u; } };
	lambda[5] = { [](double u, double x) -> double {return u * u*u*u; } };
	lambda[6] = { [](double u, double x) -> double {return exp(u); } };
	lambda[7] = { [](double u, double x) -> double {return sin(u); } };

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

	double sigma = 1;
	int coefGrid = 0;
	bool isGridUniform = true;
	bool isTimeUniform = true;

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
	cout << "Research 2.4: convergence with grid crushing" << endl;
}