#include "fem.h"


// Инициализируем модель, задавая функции u, f и тип сетки
void FEM::init(const function2D & _u, const function2D & _f, const function1D &_lambda, double _sigma, bool _isGridUniform, bool _isTimeUniform, int _condType, int _coefGrid, int _coefTime)
{
	ifstream fin("input/SLAE_parameters.txt");
	fin >> E >> delta >> maxiter;
	fin.close();
	u = _u;
	f = _f;
	lambda = _lambda;
	sigma = _sigma;
	isGridUniform = _isGridUniform;
	isTimeUniform = _isTimeUniform;
	condType = _condType;
	coefGrid = _coefGrid;
	coefTime = _coefTime;
}



// Решаем во всех узлах по времени
void FEM::solve(std::ofstream& _fout)
{
	//ofstream fout("output/solution.txt");
	// Задаём начальные условия
	q.resize(nodesCount, 0);
	qPrev.resize(nodesCount, 0);
	vector1D qExact(nodesCount);
	for (size_t i = 0; i < nodesCount; i++)
		qExact[i] = u(nodes[i].x, times[0]);

	qPrev = qExact;
	//cout << "Expected:" << endl << qPrev << endl;

	/*for (size_t i = 0; i < nodesCount; i++)
		qPrev[i] = 1;
	cout << "X_0:" << endl << qPrev << endl;*/

	//printTeXVector(fout, qPrev, coefGrid);
	//fout << endl;
	int count = 0;

	// Решаем в каждый момент временной сетки
	for (size_t i = 1; i < times.size(); i++)
	{
		dt = times[i] - times[i - 1];
		t = times[i];

		//fout << "Expected" << endl;
		//for (size_t i = 0; i < nodesCount; i++)
		//	fout << u(nodes[i].x, t) << " ";
		//fout << endl;

		do {
			qPrev = q;
			buildGlobalMatrixA(dt);
			buildGlobalVectorb();
			/*printGlobalMatrixA();
			printGlobalVectorb();*/
			calcWithLUDecomposition();
			count++;
		} while (shouldCalc(count));
		//fout << "Reality" << endl;
		//fout << q;
		//fout << endl << endl;
	}
	//qPrev = q;
	//buildGlobalMatrixA(dt);
	//buildGlobalVectorb();
	////calcWithLUDecomposition();
	//printGlobalMatrixA();
	//printGlobalVectorb();
	//cout << endl << "Result:" << endl << q << endl;
	//fout << count << "\t";
	//printTeXVector(fout, q, coefGrid);
	//fout << "\t" << calcNormE(q) << endl;
	_fout << coefGrid << "\t" 
		<< nodesCount << "\t"
		<< count << "\t"
		<< calcNormAtMainNodes(q) << endl;
	//fout.close();
}




//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


// Строим глобальную матрицу системы нелинейных уравнений
void FEM::buildGlobalMatrixA(double _dt)
{
	dt = _dt;
	A.clear();
	di.clear();
	au.clear();
	al.clear();
	A.resize(nodesCount);
	for (size_t i = 0; i < nodesCount; i++)
		A[i].resize(nodesCount, 0);

	di.resize(nodesCount, 0);
	al.resize(nodesCount - 1, 0);
	au.resize(nodesCount - 1, 0);
	for (size_t elemNumber = 0; elemNumber < finiteElementsCount; elemNumber++)
	{
		buildLocalmatrixA(elemNumber);

		//cout << ALocal << endl;
		di[elemNumber] += ALocal[0][0];		au[elemNumber] += ALocal[0][1];
		al[elemNumber] += ALocal[1][0];		di[elemNumber + 1] += ALocal[1][1];

		A[elemNumber][elemNumber] += ALocal[0][0];		A[elemNumber][elemNumber + 1] += ALocal[0][1];
		A[elemNumber + 1][elemNumber] += ALocal[1][0];	A[elemNumber + 1][elemNumber + 1] += ALocal[1][1];
	}

	// Первые краевые условия
	A[0][0] = 1; A[0][1] = 0;
	A[nodesCount - 1][nodesCount - 1] = 1; A[nodesCount - 1][nodesCount - 2] = 0;
	di[0] = 1;
	au[0] = 0;
	di[nodesCount - 1] = 1;
	al[al.size() - 1] = 0;
}



// Строим глобальный вектор правой части системы нелинейных уравнений
void FEM::buildGlobalVectorb()
{
	b.clear();
	b.resize(nodesCount, 0);

	for (size_t elemNumber = 0; elemNumber < finiteElementsCount; elemNumber++)
	{
		buildLocalVectorb(elemNumber);
		b[elemNumber] += bLocal[0];
		b[elemNumber + 1] += bLocal[1];
	}

	b[0] = u(nodes[0].x, t);
	b[nodesCount - 1] = u(nodes[nodesCount - 1].x, t);
}



// Вывод матрицы А в консоль
void FEM::printGlobalMatrixA()
{
	ofstream fout("output/A.txt");
	cout << fixed << setprecision(3);
	fout << "A = [ ";
	for (size_t i = 0; i < nodesCount; i++)
	{
		for (size_t j = 0; j < nodesCount; j++)
		{
			fout << A[i][j] << "\t";
			cout << A[i][j] << "\t";
		}
		fout << ";" << endl;
		cout << endl;
	}
	fout << "]";

	fout << endl;
	fout << "di = {" << di << "};" << endl;
	fout << "al = {" << al << "};" << endl;
	fout << "au = {" << au << "};" << endl;
	fout.close();
}



// Вывод матрицы А в консоль
void FEM::printGlobalVectorb()
{
	ofstream fout("output/b.txt");
	cout << endl << fixed << setprecision(3);
	fout << "b = [ ";
	for (size_t i = 0; i < nodesCount; i++)
	{
		fout << b[i] << ";" << endl;
		cout << b[i] << endl;
	}
	fout << "]";
	fout << endl;
	fout << "b = {" << b << "};" << endl;
	fout.close();
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


// Построение локальной матрицы жёсткости
void FEM::buildLocalMatrixG(int elemNumber)
{
	lambda0 = lambda(q[elemNumber]);
	//cout << lambda0 << endl;
	lambda1 = lambda(q[elemNumber + 1]);
	double numerator = (lambda0 + lambda1) / (2 * hx);
	GLocal[0][0] = GLocal[1][1] = numerator;
	GLocal[0][1] = GLocal[1][0] = -numerator;
}



// Построение локальной матрицы масс
void FEM::buildLocalMatrixM(int elemNumber)
{
	double numerator = (sigma * hx) / (6 * dt);
	MLocal[0][0] = MLocal[1][1] = 2 * numerator;
	MLocal[0][1] = MLocal[1][0] = numerator;
}



// Построение локальной матрицы А
void FEM::buildLocalmatrixA(int elemNumber)
{
	ALocal = GLocal = MLocal = { {0,0}, {0,0} };
	buildLocalMatrixG(elemNumber);
	buildLocalMatrixM(elemNumber);
	for (size_t i = 0; i < 2; i++)
	{
		for (size_t j = 0; j < 2; j++)
		{
			ALocal[i][j] = GLocal[i][j] + MLocal[i][j];
		}
	}
}



// Построение локального вектора b
void FEM::buildLocalVectorb(int elemNumber)
{
	bLocal = { 0, 0 };
	bLocal[0] = f(nodes[elemNumber].x, t) * hx / 2 + sigma * hx * (2 * qPrev[elemNumber] + qPrev[elemNumber + 1]) / (6 * dt);
	bLocal[1] = f(nodes[elemNumber + 1].x, t) * hx / 2 + sigma * hx * (qPrev[elemNumber] + 2 * qPrev[elemNumber + 1]) / (6 * dt);
}
