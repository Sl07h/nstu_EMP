#include "fem.h"


void main() {

	int testsCount = 7;
	vector <function1D> function_u(testsCount), function_f(testsCount);
	function_u[0] = { [](double x) -> double { return x; } };
	function_f[0] = { [](double x) -> double { return 0; } };

	function_u[1] = { [](double x) -> double { return 10 * x; } };
	function_f[1] = { [](double x) -> double { return 0; } };

	function_u[2] = { [](double x) -> double { return pow(x,2); } };
	function_f[2] = { [](double x) -> double { return 2; } };

	function_u[3] = { [](double x) -> double { return  pow(x,3); } };
	function_f[3] = { [](double x) -> double { return 6 * x; } };

	function_u[4] = { [](double x) -> double { return pow(x,4); } };
	function_f[4] = { [](double x) -> double { return 12 * pow(x,2); } };

	function_u[5] = { [](double x) -> double { return sin(x); } };
	function_f[5] = { [](double x) -> double { return -sin(x); } };

	function_u[6] = { [](double x) -> double { return exp(x); } };
	function_f[6] = { [](double x) -> double { return exp(x); } };

	string prefix = "";
	int coefGrid = 0;

	bool isGridUniform = true;
	bool isTimeUniform = false;

	if (!isGridUniform)
		prefix = "Non";

	string gridFile = "grids/" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";
	string gridBorderFile = "grids/Border" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";

	FEM fem;
	fem.init(function_u[0], function_f[1], isGridUniform, isTimeUniform, 1, 0, 0);
	fem.inputGrid();
	fem.buildGrid();

	fem.inputTime();
	fem.buildTimeGrid();
	
	fem.saveGridAndBorder(gridFile, gridBorderFile);
	fem.solve();

	string runVisualisation = "python plot.py " + gridFile + " " + gridBorderFile;
	system(runVisualisation.c_str());


}
