#include "fdm.h"



// Исследование работы МКР на функциях u,f
double performSingleTest(function2D &u, function2D&f, int index, bool isGridUniform, int condNumber, bool useJacopbNotGaussSeidel, int coef = 0) {

	string absResidualFile = "tables/Uniform_" + to_string(index) + "_" + to_string(condNumber) + "_" + to_string(coef) + "_AbsResidual.txt";

	FDM fdm;
	fdm.init(u, f, isGridUniform, condNumber, coef);
	fdm.inputEquationParameters();
	fdm.inputSLAEParameters();
	fdm.inputGrid();
	fdm.buildGrid();
	//fdm.showGrid();
	fdm.transformGridToSLAE();
	fdm.convMatrixToDense();
	//fdm.outputSLAE("tables/Uniform_" + to_string(index) + "_A.txt");
	fdm.generateInitualGuess();
	cout << "Count of steps: " << fdm.calcIterative(useJacopbNotGaussSeidel, 0.8) << endl << endl;
	return fdm.calcAbsResidual(absResidualFile);
}



// Исследование сходимости при дроблении сетки
void exploreConvergence(function2D &u, function2D&f, int index, bool isGridUniform, int condNumber, bool useJacopbNotGaussSeidel) {

	ofstream fout;
	string prefix = "";
	if (!isGridUniform)
		prefix = "Non";

	fout.open("tables/Convergence" + prefix + "Uniform" + to_string(index) + "_" + to_string(condNumber) + ".txt");
	fout << std::scientific;

	// Увеличиваем сетку в kx^coef раз
	for (size_t coef = 0; coef < 4; coef++)
		fout << performSingleTest(u, f, index, isGridUniform, condNumber, useJacopbNotGaussSeidel, coef) << endl;

	fout.close();
}



// Отрисовка сетки
void drawGrids(bool isGridUniform) {

	string prefix = "";
	if (!isGridUniform)
		prefix = "Non";

	function2D u = { [](double x, double y) -> double { return x + y; } };
	function2D f = { [](double x, double y) -> double { return 0; } };

	for (size_t coef = 0; coef < 4; coef++)
	{
		string gridFile = "grids/" + prefix + "Uniform_" + to_string(coef) + ".txt";
		string gridBorderFile = "grids/Border" + prefix + "Uniform_" + to_string(coef) + ".txt";

		FDM fdm;
		fdm.init(u, f, isGridUniform, 1, coef);
		//fdm.inputEquationParameters();
		//fdm.inputSLAEParameters();
		fdm.inputGrid();
		fdm.buildGrid();
		fdm.saveGridAndBorder(gridFile, gridBorderFile);

		string runVisualisation = "python plot.py " + gridFile + " " + gridBorderFile;
		system(runVisualisation.c_str());
	}
}




int main() {

	int testsCount = 7;
	bool isGridUniform;

	vector <function2D> function_u(testsCount), function_f(testsCount);
	function_u[0] = { [](double x, double y) -> double { return x + y; } };
	function_f[0] = { [](double x, double y) -> double { return 0; } };

	function_u[1] = { [](double x, double y) -> double { return 10 * x + y; } };
	function_f[1] = { [](double x, double y) -> double { return 0; } };

	function_u[2] = { [](double x, double y) -> double { return pow(x,2) + pow(y,2); } };
	function_f[2] = { [](double x, double y) -> double { return 4; } };

	function_u[3] = { [](double x, double y) -> double { return  pow(x,3) + pow(y,3); } };
	function_f[3] = { [](double x, double y) -> double { return 6 * (x + y); } };

	function_u[4] = { [](double x, double y) -> double { return pow(x,4) + pow(y,4); } };
	function_f[4] = { [](double x, double y) -> double { return 12 * (pow(x,2) + pow(y,2)); } };

	function_u[5] = { [](double x, double y) -> double { return sin(x) + cos(y); } };
	function_f[5] = { [](double x, double y) -> double { return -sin(x) - cos(y); } };

	function_u[6] = { [](double x, double y) -> double { return exp(x + y); } };
	function_f[6] = { [](double x, double y) -> double { return 2 * exp(x + y); } };


	// Таблица 0
	// Отрисовка сеток
	drawGrids(true);
	drawGrids(false);



	// Таблица 1
	// Точность на разных функциях
	/*std::ofstream fout1("tables/table1.txt");
	std::ofstream fout3("tables/table3.txt");
	fout1 << std::scientific;
	fout3 << std::scientific;
	for (size_t i = 0; i < testsCount; i++) {
		fout1 << performSingleTest(function_u[i], function_f[i], i, true, 1, true, 0) << endl;
		fout3 << performSingleTest(function_u[i], function_f[i], i, true, 3, true, 0) << endl;
	}
	fout1.close();
	fout3.close();*/


	// Таблица 2
	// Точность в зависимости от дробления сетки
	/*for (size_t i = 0; i < testsCount; i++) {
		exploreConvergence(function_u[i], function_f[i], i, true, 1, true);
		exploreConvergence(function_u[i], function_f[i], i, false, 1, true);
		exploreConvergence(function_u[i], function_f[i], i, true, 3, true);
		exploreConvergence(function_u[i], function_f[i], i, false, 3, true);
	}*/



	return 0;
}
