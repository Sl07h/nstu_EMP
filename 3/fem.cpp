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
	function1D _u_s,
	function1D _u_c,
	function1D _f_s,
	function1D _f_c,
	double _lambda,
	double _sigma,
	double _omega,
	double _hi,
	bool _isGridUniform,
	bool _isTimeUniform,
	int _condType,
	int _coefGrid,
	int _coefTime
) {
	ifstream fin("input/SLAE_parameters.txt");
	fin >> E >> delta >> maxiter;
	fin.close();
	u_s = _u_s;
	u_c = _u_c;
	f_s = _f_s;
	f_c = _f_c;
	lambda = _lambda;
	sigma = _sigma;
	omega = _omega;
	hi = _hi;
	isGridUniform = _isGridUniform;
	isTimeUniform = _isTimeUniform;
	condType = _condType;
	coefGrid = _coefGrid;
	coefTime = _coefTime;
}

pair<int, double> FEM::solve(int solver)
{
	pair<int, double> result;
	n = 2 * nodesCount;
	buildGlobalMatrixA();
	buildGlobalVectorb();

	switch (solver)
	{
	case 1:
		result = LOSfactLUsq();
		break;
	case 2:
		result = BiCG();
		break;
	case 3:
		//result = BiCG();
		break;
	default:
		break;
	}
	return result;
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Строим глобальную матрицу системы нелинейных уравнений
void FEM::buildGlobalMatrixA()
{
	di.clear();
	au.clear();
	al.clear();
	/*A.clear();
	A.resize(2 * nodesCount, 0);
	for (size_t i = 0; i < 2 * nodesCount; i++)
	{
		A[i].resize(2 * nodesCount, 0);
	}*/

	di.resize(2 * nodesCount, 0);
	ia.resize(2 * nodesCount + 1);
	ia[0] = 0; ia[1] = 0; ia[2] = 1;
	double iaLast = 1;
	for (size_t i = 3; i < ia.size(); i++)
	{
		if (i % 2 != 0)
			iaLast += 2;
		else
			iaLast += 3;
		ia[i] = iaLast;
	}
	ja.resize(5 * finiteElementsCount + 1, 0);
	ja[0] = 0;
	for (size_t i = 1; i < ja.size(); i += 5)
	{
		ja[i] = ja[i - 1];
		ja[i + 1] = ja[i - 1] + 1;
		ja[i + 2] = ja[i - 1];
		ja[i + 3] = ja[i - 1] + 1;
		ja[i + 4] = ja[i - 1] + 2;
	}
	al.resize(5 * finiteElementsCount + 1, 0);
	au.resize(5 * finiteElementsCount + 1, 0);
	for (size_t elemNumber = 0; elemNumber < 2 * finiteElementsCount; elemNumber += 2)
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
			if (t == 0)
				i0 = i1;
			else if (t == 1)
				i0 = i1 - 1;
			for (size_t k = i0; k < i1; k++, tLocal++)
			{
				al[k] += ALocal[t][tLocal];
				au[k] += ALocal[tLocal][t];
			}
		}

		/*outputALocal();
		convAToDense();
		outputA();*/
	}

	// Первые краевые условия
	di[0] = 1; di[1] = 1; al[0] = 0;
	for (size_t i = 0; i < 5; i++)
		au[i] = 0;

	di[di.size() - 1] = 1; di[di.size() - 2] = 1; au[au.size() - 1] = 0;
	for (size_t i = 1; i < 6; i++)
		al[al.size() - i] = 0;

	/*convAToDense();
	cout << endl << "Added 1st boundary conditions:" << endl;
	outputA();*/
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
	b[0] = u_c(nodes[0].x);
	b[1] = u_s(nodes[0].x);
	b[b.size() - 2] = u_s(nodes[nodes.size() - 1].x);
	b[b.size() - 1] = u_c(nodes[nodes.size() - 1].x);
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



// Построение локальной матрицы А
void FEM::buildLocalmatrixA(int elemNumber)
{
	p00 = lambda / (hx*hx) - (omega*omega)*hi / 3;
	p11 = p00;
	p01 = -lambda / (hx*hx) - (omega*omega)*hi / 6;
	p10 = p01;

	c00 = omega * sigma / 3;
	c11 = c00;
	c01 = omega * sigma / 6;
	c10 = c01;

	ALocal = { {p00,	-c00,	p01,	-c01},
				{c00,	p00,	c01,	p01},
				{p10,	-c10,	p11,	-c11},
				{c10,	p10,	c11,	p11} };
}


// Построение локального вектора b
void FEM::buildLocalVectorb(int elemNumber)
{
	bLocal = { 0, 0, 0, 0 };
	double f0_s = f_s(nodes[elemNumber].x);
	double f0_c = f_c(nodes[elemNumber].x);
	double f1_s = f_s(nodes[elemNumber + 1].x);
	double f1_c = f_c(nodes[elemNumber + 1].x);
	bLocal[0] = hx * (2 * f0_s + f1_s) / 6;
	bLocal[1] = hx * (2 * f0_c + f1_c) / 6;
	bLocal[2] = hx * (f0_s + 2 * f1_s) / 6;
	bLocal[3] = hx * (f0_c + 2 * f1_c) / 6;
}