#include "fem.h"


// Инициализируем модель, задавая функции u, f и тип сетки
void FEM::init(const function1D & _u, const function1D & _f, const function2D &_lambda, bool _isGridUniform, bool _isTimeUniform, int _condType, int _coefGrid, int _coefTime)
{
	ifstream fin("input/SLAE_parameters.txt");
	fin >> E >> delta >> maxiter;
	fin.close();
	u = _u;
	f = _f;
	lambda = _lambda;
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
	fout << qPrev << endl;


	// Решаем в каждый момент временной сетки
	for (size_t i = 1; i < times.size(); i++)
	{
		dt = times[i] - times[i - 1];
		double t = times[i];
		int count = 0;
		
		do {
			buildGlobalMatrixA(dt);
			buildGlobalVectorb();
			//printGlobalMatrixA();
			//printGlobalVectorb();
			calcWithLUDecomposition();
			qPrev = q;
			/*cout << q << endl;
			cout << "Time Iteration: " << i
				<< " Count: " << count << endl;
			fout << q << endl;*/
			count++;
		} while (shouldCalc(count));
		fout << q << endl;
	}
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

	b[0] = u(nodes[0].x);
	b[nodesCount - 1] = u(nodes[nodesCount - 1].x);
}



// Вывод матрицы А в консоль
void FEM::printGlobalMatrixA()
{
	ofstream fout("output/A.txt");
	cout << fixed << setprecision(2);
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
	cout << endl << fixed << setprecision(2);
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
	lambda0 = lambda(q[elemNumber], nodes[elemNumber].x);
	//cout << lambda0 << endl;
	lambda1 = lambda(q[elemNumber + 1], nodes[elemNumber + 1].x);
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
	bLocal[0] = f(nodes[elemNumber].x) * hx / 2 + sigma * hx * (2 * qPrev[elemNumber] + qPrev[elemNumber + 1]) / (6 * dt);
	bLocal[1] = f(nodes[elemNumber + 1].x) * hx / 2 + sigma * hx * (qPrev[elemNumber] + 2 * qPrev[elemNumber + 1]) / (6 * dt);
}
