#include "fdm.h"



//#include <Eigen/Eigen>
//#include <Eigen/Dense>
// Включаемые каталоги
// C:\Users\User\Desktop\eigen-eigen-323c052e1731
//typedef Eigen::MatrixXd Matrix;
//typedef Eigen::VectorXd Vector;






int main()
{
	int testsCount = 5;
	FDM fdm;


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



	for (size_t i = 0; i < testsCount; i++)
	{
		bool isGridUniform = true;
		fdm.init(function_u[i], function_f[i], i, isGridUniform);
		fdm.inputEquationParameters();
		fdm.inputSLAEParameters();
		fdm.inputGrid();
		fdm.buildGrid();
		fdm.transformGridToSLAE();
		//fdm.writeDenseMatrixToFile("report/A.txt");
		//fdm.convMatrixToDense();
		//fdm.writeDenseMatrixToFile("report/A2.txt");
		double w = fdm.findOptimalW(1);
		cout << w << endl;
		//fdm.generateInitualGuess();
		//fdm.Jacobi(w);
		//fdm.writeXToFile("report/x.txt");
		//fdm.writebToFile("report/b.txt");
		fdm.calcAbsResidual("report/" + to_string(i) + "_" + to_string(int(isGridUniform)) + "_" + "abs_residual.txt");
	}


	for (size_t i = 0; i < testsCount; i++)
	{
		bool isGridUniform = false;
		fdm.init(function_u[i], function_f[i], i, isGridUniform);
		fdm.inputEquationParameters();
		fdm.inputSLAEParameters();
		fdm.inputGrid();
		fdm.buildGrid();
		fdm.transformGridToSLAE();
		//fdm.writeDenseMatrixToFile("report/A.txt");
		//fdm.convMatrixToDense();
		//fdm.writeDenseMatrixToFile("report/A2.txt");
		double w = fdm.findOptimalW(1);
		cout << w << endl;
		//fdm.generateInitualGuess();
		//fdm.Jacobi(w);
		//fdm.writeXToFile("report/x.txt");
		//fdm.writebToFile("report/b.txt");
		fdm.calcAbsResidual("report/" + to_string(i) + "_" + to_string(int(isGridUniform)) + "_" + "abs_residual.txt");
	}


	return 0;
}
