#include "fem.h"


// Инициализируем модель, задавая функции u, f и тип сетки
void FEM::init(function1D & _u, function1D & _f, bool _isGridUniform, bool _isTimeUniform, int _condType, int _coefGrid, int _coefTime)
{
	ifstream fin("input/SLAE_parameters.txt");
	fin >> E >> delta >> maxiter;
	fin.close();
	u = _u;
	f = _f;
	isGridUniform = _isGridUniform;
	isTimeUniform = _isTimeUniform;
	condType = _condType;
	coefGrid = _coefGrid;
	coefTime = _coefTime;
}



// Решаем во всех узлах по времени
void FEM::solve()
{
	ofstream fout("output/solution.txt");
	// Задаём начальные условия
	q.resize(nodesCount, 0);
	qPrev.resize(nodesCount, 0);
	for (size_t i = 0; i < nodesCount; i++)
		qPrev[i] = u(nodes[i].x);


	// Решаем в каждый момент временной сетки
	for (size_t i = 1; i < times.size(); i++)
	{
		dt = times[i] - times[i - 1];
		double t = times[i];
		
		bool doCalculation = true;
		while (doCalculation) {

			buildGlobalMatrixA();
			buildGlobalVectorb();
			//printGlobalMatrixA();
			//printGlobalVectorb();

			calcWithLUDecomposition();
			qPrev = q;

			if (shouldCalc(i) == false) {
				cout << "Iteration: " << i
					<< "Time: " << t << endl;
				fout << q << endl;
				doCalculation = false;
			}
		}
	}
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


// Строим глобальную матрицу системы нелинейных уравнений
void FEM::buildGlobalMatrixA()
{
	A.resize(nodesCount);
	for (size_t i = 0; i < nodesCount; i++)
		A[i].resize(nodesCount, 0);

	di.resize(nodesCount, 0);
	al.resize(nodesCount - 1, 0);
	au.resize(nodesCount - 1, 0);
	for (size_t elemNumber = 0; elemNumber < finiteElementsCount; elemNumber++)
	{
		buildLocalmatrixA(elemNumber);

		di[elemNumber] += ALocal[0][0];
		di[elemNumber + 1] += ALocal[0][1];
		al[elemNumber] += ALocal[1][0];
		au[elemNumber] += ALocal[0][1];

		A[elemNumber][elemNumber] += ALocal[0][0];
		A[elemNumber][elemNumber + 1] += ALocal[0][1];
		A[elemNumber + 1][elemNumber] += ALocal[1][0];
		A[elemNumber + 1][elemNumber + 1] += ALocal[1][1];
	}
}



// Строим глобальный вектор правой части системы нелинейных уравнений
void FEM::buildGlobalVectorb()
{
	b.resize(nodesCount, 0);

	for (size_t elemNumber = 0; elemNumber < finiteElementsCount; elemNumber++)
	{
		buildLocalVectorb(elemNumber);
		b[elemNumber] += bLocal[0];
		b[elemNumber + 1] += bLocal[1];
	}
}



// Вывод матрицы А в консоль
void FEM::printGlobalMatrixA()
{
	cout << fixed << setprecision(2);
	for (size_t i = 0; i < nodesCount; i++)
	{
		for (size_t j = 0; j < nodesCount; j++)
		{
			cout << A[i][j] << "\t";
		}
		cout << endl;
	}
}



// Вывод матрицы А в консоль
void FEM::printGlobalVectorb()
{
	cout << endl << fixed << setprecision(2);
	for (size_t i = 0; i < nodesCount; i++)
	{
		cout << b[i] << endl;
	}
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


// Построение локальной матрицы жёсткости
void FEM::buildLocalMatrixG(int elemNumber)
{
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
	bLocal = { 0,0 };
	bLocal[0] = f(nodes[elemNumber].x) * hx / 2 + sigma * hx / (6 * dt) * (2 * qPrev[0] + qPrev[1]);
	bLocal[1] = f(nodes[elemNumber + 1].x) * hx / 2 + sigma * hx / (6 * dt) * (qPrev[0] + 2 * qPrev[1]);
}
