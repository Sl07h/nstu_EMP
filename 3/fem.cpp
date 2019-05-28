#include "fem.h"


// Инициализируем модель, задавая функции u, f и тип сетки
void FEM::init(
	const function1D & _u_s,
	const function1D & _u_c,
	const function1D & _f_s,
	const function1D & _f_c,
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

void FEM::solve()
{
	buildGlobalMatrixA();
	buildGlobalVectorb();

	LOS();
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


	di.resize(2 * nodesCount, 0);
	ia.resize(5 * nodesCount + 1);
	ja.resize(5 * nodesCount + 1, 0);
	al.resize(5 * nodesCount + 1, 0);
	au.resize(5 * nodesCount + 1, 0);
	for (size_t elemNumber = 0; elemNumber < finiteElementsCount; elemNumber++)
	{
		buildLocalmatrixA(elemNumber);

		int k = 0;
		for (size_t i = elemNumber; i < elemNumber + 4; i++, k++)
			di[i] += ALocal[k][k];

		for (size_t i = elemNumber; i < elemNumber + 4; i++)
		{
			int i0 = ia[i];
			int i1 = ia[i + 1];
			int j = ja[i0];
			for (size_t k = i0; k < i1; k++, j++)
			{
				al[k] += ALocal[i][j];
				au[k] += ALocal[j][i];
			}
		}
	}

	// Первые краевые условия
	di[0] = 1;
	au[0] = 0;
	au[1] = 0;
	au[3] = 0;

	di[di.size() - 1] = 1;
	al[al.size() - 1] = 0;
	al[al.size() - 2] = 0;
	al[al.size() - 3] = 0;
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

	/*b[0] = u(nodes[0].x, t);
	b[nodesCount - 1] = u(nodes[nodesCount - 1].x, t);*/
}

/*
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
*/

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



// Построение локальной матрицы А
void FEM::buildLocalmatrixA(int elemNumber)
{
	p00 = lambda / (hx*hx) - (omega*omega)*hi / 3;
	p01 = -lambda / (hx*hx) - (omega*omega)*hi / 6;
	c00 = omega * sigma / 6;
	c01 = omega * sigma / 6;
	ALocal = { {p00,	-c00,	p01,	-c01},
				{c00,	p00,	c01,	p01},
				{p01,	-c01,	p00,	-c00},
				{c01,	p01,	c00,	p00} };
}


// Построение локального вектора b
void FEM::buildLocalVectorb(int elemNumber)
{
	bLocal = { 0,0,0,0 };
	bLocal[0] = f_s(nodes[elemNumber].x);
	bLocal[1] = f_c(nodes[elemNumber].x);
	bLocal[2] = f_s(nodes[elemNumber + 1].x);
	bLocal[3] = f_c(nodes[elemNumber + 1].x);
}