#include "fem.h"



function1D calcFirstDerivative(const function1D& f) {
	return [f](double x) -> double {
		const double h = 0.00001;
		return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
	};
}


function2D calcRightPart(const function1D& lambda, const function2D& u, double sigma) {
	return [=](double x, double t) -> double {
		using namespace std::placeholders;
		auto duBydt = calcFirstDerivative(std::bind(u, x, _1));
		auto duBydx = calcFirstDerivative(std::bind(u, _1, t));
		auto lambda_grad = [=](double x, double t) -> double {
			return lambda(u(x, t)) * duBydx(x);
		};
		auto div = calcFirstDerivative(std::bind(lambda_grad, _1, t));
		return -div(x) + sigma * duBydt(x);
	};
}



void main() {
	double sigma = 1;

	vector <function2D> function_u(6), function_f(6);
	function_u[0] = { [](double x, double t) -> double { return /*t + */x; } };
	function_u[1] = { [](double x, double t) -> double { return /*t + */10 * x; } };
	function_u[2] = { [](double x, double t) -> double { return /*t + */pow(x, 2); } };
	function_u[3] = { [](double x, double t) -> double { return /*t + */pow(x, 3); } };
	function_u[4] = { [](double x, double t) -> double { return /*t + */sin(x); } };
	function_u[5] = { [](double x, double t) -> double { return /*t + */exp(x); } };
	vector <function1D>  function_lambda(3);
	function_lambda[0] = { [](double u) -> double {return 1; } };
	/*function_lambda[1] = { [](double u, double x) -> double {return u + x; } };
	function_lambda[2] = { [](double u, double x) -> double {return u * u; } };*/
	cout << "Ready:" << endl;
	// Зависимость от функции u
	for (int i = 0; i < 6; i++) {
		function_f[i] = calcRightPart(function_lambda[0], function_u[i], sigma);

		string prefix = "";

		bool isGridUniform = true;
		bool isTimeUniform = true;

		if (!isGridUniform)
			prefix = "Non"; 

		ofstream fout("report/file_u" + to_string(i) + ".txt");
		fout << scientific << fixed;
		fout << "i\tnodes\titers\tnorm\n";
		for (int coefGrid = 0; coefGrid < 8; coefGrid++)
		{
			//cout << float(i * 5 + coefGrid) / 30 << "\r";
			string gridFile = "grids/" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";
			string gridBorderFile = "grids/Border" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";

			FEM fem;
			fem.init(function_u[i], function_f[i], function_lambda[0], sigma, isGridUniform, isTimeUniform, 1, coefGrid, 0);
			fem.inputGrid();
			fem.buildGrid();

			fem.inputTime();
			fem.buildTimeGrid();

			/*fem.saveGridAndBorder(gridFile, gridBorderFile);
			string runVisualisation = "python plot.py " + gridFile + " " + gridBorderFile;
			system(runVisualisation.c_str());*/
			fem.solve(fout);

		}
		fout.close();
	}

	//string runVisualisation = "python plot.py " + gridFile + " " + gridBorderFile;
	//system(runVisualisation.c_str());
}
