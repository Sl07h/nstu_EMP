#include "fem.h"





void FEM::outputALocal() {
	cout << endl;
	for (int i = 0; i < ALocal.size(); ++i) {
		for (int j = 0; j < ALocal.size(); ++j) {
			cout << ALocal[i][j] << "  ";
		}
		cout << endl;
	}
}

void FEM::convAToDense() {

	A.clear();
	A.resize(n);
	for (int i = 0; i < n; ++i) {
		A[i].resize(n, 0);
	}

	for (int i = 0; i < n; ++i) {

		A[i][i] = di[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];

		for (int k = i0; k < i1; ++k) {
			A[i][ja[k]] = al[k];
			A[ja[k]][i] = au[k];
		}
	}
}




void FEM::outputA() {
	cout << endl;
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[i][j] << "\t";
		}
		cout << endl;
	}
}






// Инициализируем модель, задавая функции u, f и тип сетки
void FEM::init(
	function1D _u,
	function1D _f,
	double _lambda,
	double _gamma,
	double _sigma,
	double _chi,
	bool _isGridUniform,
	bool _isTimeUniform,
	int _coefGrid,
	int _coefTime
) {
	ifstream fin("input/SLAE_parameters.txt");
	fin >> E >> delta >> maxiter;
	fin.close();
	u = _u;
	f = _f;
	lambda = _lambda;
	gamma = _gamma;
	sigma = _sigma;
	chi = _chi;
	isGridUniform = _isGridUniform;
	isTimeUniform = _isTimeUniform;
	coefGrid = _coefGrid;
	coefTime = _coefTime;
}

void FEM::solve()
{
	n = 2 * nodesCount;
	buildGlobalMatrixA();
	buildGlobalVectorb();


	auto result = LOS();
	cout << "Iters count: " << result.first << endl;
	//BiCG();

	cout << x << endl;

	//double result = 0.0;
	for (size_t i = 0; i < nodesCount; i++)
	{
		//result += pow((x[2 * i] - u_c(nodes[i].x)), 2);

		cout << u(nodes[i].x) << "\t";
	}
	//cout << sqrt(result);
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Строим глобальную матрицу системы нелинейных уравнений
void FEM::buildGlobalMatrixA()
{
	A.clear();
	A.resize(nodesCount);
	for (size_t i = 0; i < 2 * nodesCount; i++)
		A[i].resize(2 * nodesCount, 0);

	for (size_t elemNumber = 0; elemNumber < finiteElementsCount; elemNumber++)
	{
		buildLocalmatrixA(elemNumber);

		int t = 0;
		for (size_t i = elemNumber; i < elemNumber + 4; i++, t++)
			di[i] += ALocal[t][t];

		t = 0;
		for (size_t i = elemNumber; i < elemNumber + 4; i++, t++)
		{
			int i0 = ia[i];
			int i1 = ia[i + 1];
			int tLocal = 0;
			if (t < 2)
				i0 = i1 - 1;
			for (size_t k = i0; k < i1; k++, tLocal++)
			{
				al[k] += ALocal[t][tLocal];
				au[k] += ALocal[tLocal][t];
			}
		}

		outputALocal();
		convAToDense();
		outputA();
	}

	// Первые краевые условия
	di[0] = 1; di[1] = 1; al[0] = 0;
	for (size_t i = 0; i < 5; i++)
		au[i] = 0;

	di[di.size() - 1] = 1; di[di.size() - 2] = 1; au[au.size() - 1] = 0;
	for (size_t i = 1; i < 6; i++)
		al[al.size() - i] = 0;

	convAToDense();
	cout << endl << "Added 1st boundary conditions:" << endl;
	outputA();
}




// Строим глобальный вектор правой части системы нелинейных уравнений
void FEM::buildGlobalVectorb()
{
	b.clear();
	b.resize(2 * nodesCount, 0);

	for (size_t elemNumber = 0; elemNumber < finiteElementsCount; elemNumber++)
	{
		buildLocalVectorb(elemNumber);
		int k = 0;
		for (size_t i = elemNumber; i < elemNumber + 4; i++, k++)
			b[i] = bLocal[k];
	}

	// Первые краевые условия
	b[0] = u(nodes[0].x);
	b[b.size() - 1] = u(nodes[b.size() - 1].x);
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


// Построение локальной матрицы
void FEM::buildLocalmatrixA()
{



}



// Построение локальной матрицы жёсткости
void FEM::buildLocalmatrixG()
{
	GLocal = { {0,0,0,0},
				{0,0,0,0},
				{0,0,0,0},
				{0,0,0,0} };

	tmpLocal = { {2, -2, 1, -1},
				{-2, 2, -1, 1},
				{1, -1, 2, -2},
				{-1, 1, -2, 2} };

	for (size_t i = 0; i < 4; i++)
		for (size_t j = 0; j < 4; j++)
			GLocal[i][j] += tmpLocal[i][j] * lambda * hy / (6 * hx);

	tmpLocal = { {2, 1, -2, -1},
				{1, 2, -1, -2},
				{-2, -1, 2, 1},
				{-1, -2, 1, 2} };


	for (size_t i = 0; i < 4; i++)
		for (size_t j = 0; j < 4; j++)
			GLocal[i][j] += tmpLocal[i][j] * lambda * hx / (6 * hy);
}


// Построение локальной матрицы масс
void FEM::buildLocalmatrixM()
{
	MLocal = { {4, 2, 2, 1},
				{2, 4, 1, 2},
				{2, 1, 4, 2},
				{1, 2, 2, 4} };

	for (size_t i = 0; i < 4; i++)
		for (size_t j = 0; j < 4; j++)
			MLocal[i][j] *= hx * hy / 36;
}


// Построение локального вектора b
void FEM::buildLocalVectorb(int elemNumber)
{
	bLocal = { 0, 0, 0, 0 };
	f1 = f(nodes[elemNumber].x);
	f2 = f(nodes[elemNumber + 1].x);
	f3 = f(nodes[elemNumber + 2].x);
	f4 = f(nodes[elemNumber + 3].x);

	bLocal[0] = 4 * f1 + 2 * f2 + 2 * f3 + 1 * f4;
	bLocal[1] = 2 * f1 + 4 * f2 + 1 * f3 + 2 * f4;
	bLocal[2] = 2 * f1 + 1 * f2 + 4 * f3 + 2 * f4;
	bLocal[3] = 1 * f1 + 2 * f2 + 2 * f3 + 4 * f4;
}