#include "fdm.h"


//#include <Eigen/Eigen>
//#include <Eigen/Dense>
// Включаемые каталоги
// C:\Users\User\Desktop\eigen-eigen-323c052e1731
//typedef Eigen::MatrixXd Matrix;
//typedef Eigen::VectorXd Vector;


// Исследование работы МКР на функциях u,f
void performSingleTest(function2D &u, function2D&f, int index, bool isGridUniform, int condNumber, int method, ofstream &fout) {

	FDM fdm;
	fdm.init(u, f, isGridUniform, condNumber, 1);
	fdm.inputEquationParameters();
	fdm.inputSLAEParameters();
	fdm.inputGrid();
	fdm.buildGrid();
	fdm.showGrid();

	fdm.transformGridToSLAE();
	fdm.convMatrixToDense();
	//fdm.outputSLAE("tables/SLAE_" + to_string(index) + "_" + to_string(condNumber) + "_AbsResidual.txt");
	double w = fdm.findOptimalW(method);
	cout << "Optimal w = " << w << endl;
	//fdm.writeDenseMatrixToFile("tables/A.txt");
	//fdm.writeSecondDenseMatrixToFile("tables/A2.txt");
	if (isGridUniform)
		fout << fdm.calcAbsResidual("tables/Uniform_" + to_string(index) + "_" + to_string(condNumber) + "_AbsResidual.txt") << endl;
	else
		fout << fdm.calcAbsResidual("tables/NonUniform_" + to_string(index) + "_" + to_string(condNumber) + "_AbsResidual.txt") << endl;

	cout << endl;
}


void exploreConvergence(function2D &u, function2D&f, int index, bool isGridUniform, int condNumber, bool useJacopbNotGaussSeidel, bool doVisualisation) {

	ofstream fout;
	string prefix = "";
	if (!isGridUniform)
		prefix = "Non";

	fout.open("tables/Convergence" + prefix + "Uniform" + to_string(index) + "_" + to_string(condNumber) + ".txt");


	for (size_t coef = 0; coef < 4; coef++)
	{
		string gridFile = "grids/" + prefix + "Uniform_" + to_string(index) + "_" + to_string(condNumber) + "_" + to_string(coef) + ".txt";
		string gridBorderFile = "grids/Border" + prefix + "Uniform_" + to_string(index) + "_" + to_string(condNumber) + "_" + to_string(coef) + ".txt";
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


		fdm.checkAnswer();

		fdm.saveGridAndBorder(gridFile, gridBorderFile);
		fout << fdm.calcAbsResidual(absResidualFile) << endl;
		if (doVisualisation) {
			string runVisualisation = "python plot.py " + gridFile + " " + gridBorderFile;
			system(runVisualisation.c_str());
		}
	}

	fout.close();
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


	//cout << setprecision(1) << fixed;
	//int i = 0;
	//exploreConvergence(function_u[i], function_f[i], i, true, 1, true, true);
	//exploreConvergence(function_u[i], function_f[i], i, false, 1, true, true);
	//exploreConvergence(function_u[i], function_f[i], i, true, 3, true, true);
	//exploreConvergence(function_u[i], function_f[i], i, false, 3, true, true);

	//ofstream fout("tables/table.txt");
	//fout << setprecision(numeric_limits<double>::digits10 + 1) << fixed;
	
	for (size_t i = 0; i < testsCount - 3; i++) {

		exploreConvergence(function_u[i], function_f[i], i, true, 1, true, true);
		exploreConvergence(function_u[i], function_f[i], i, false, 1, true, true);
		exploreConvergence(function_u[i], function_f[i], i, true, 3, true, true);
		exploreConvergence(function_u[i], function_f[i], i, false, 3, true, true);
	}

	//fout.close();

	return 0;
}
