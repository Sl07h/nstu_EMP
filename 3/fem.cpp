#include "fem.h"





void FEM::outputALocal() {
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
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A.size(); ++j) {
			cout << A[i][j] << "  ";
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

void FEM::solve()
{
	n = 2 * nodesCount;
	buildGlobalMatrixA();
	buildGlobalVectorb();


	//LOS();
	BiCG();

	cout << x << endl;

	for (size_t i = 0; i < nodesCount; i++)
	{
		cout << u_s(nodes[i].x) << "\t";
		cout << u_c(nodes[i].x) << "\t";
	}
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
	ia.resize(5 * nodesCount + 1);
	ja.resize(5 * nodesCount + 1, 0);
	al.resize(5 * nodesCount + 1, 0);
	au.resize(5 * nodesCount + 1, 0);
	for (size_t elemNumber = 0; elemNumber < nodesCount; elemNumber+=2)
	{
		buildLocalmatrixA(elemNumber);
		
		/*int k = 0;
		for (size_t i = elemNumber * 2; i < elemNumber * 2 + 2; i++, k++) {
			int kInternal = i * 2;
			for (size_t j = elemNumber*2; j < elemNumber *2+ 2; j++, kInternal++)
			{
				A[i][kInternal]
			}
		}*/


		int t = 0;
		for (size_t i = elemNumber; i < elemNumber + 4; i++, t++)
			di[i] += ALocal[t][t];

		t = 0;
		for (size_t i = elemNumber; i < elemNumber + 4; i++, t++)
		{
			int i0 = ia[i];
			int i1 = ia[i + 1];
			int j = ja[i0];
			int tLocal = 0;
			for (size_t k = i0; k < i1; k++, j++)
			{
				al[k] += ALocal[i][j];
				au[k] += ALocal[j][i];
			}
		}

		outputALocal();
		convAToDense();
		outputA();
	}

	// Первые краевые условия
	di[0] = 1; di[1] = 1;
	for (size_t i = 0; i < 5; i++)
		au[i] = 0;


	di[di.size() - 1] = 1; di[di.size() - 2] = 1;
	for (size_t i = 1; i < 6; i++)
		al[al.size() - i];
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