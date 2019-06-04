#include "fem.h"





//void FEM::outputALocal() {
//	cout << endl;
//	for (int i = 0; i < ALocal.size(); ++i) {
//		for (int j = 0; j < ALocal.size(); ++j) {
//			cout << ALocal[i][j] << "  ";
//		}
//		cout << endl;
//	}
//}

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


void FEM::outputG() {
	cout << endl;
	for (int i = 0; i < G.size(); ++i) {
		for (int j = 0; j < G.size(); ++j) {
			cout << G[i][j] << "\t";
		}
		cout << endl;
	}
}

void FEM::outputM() {
	cout << endl;
	for (int i = 0; i < M.size(); ++i) {
		for (int j = 0; j < M.size(); ++j) {
			cout << M[i][j] << "\t";
		}
		cout << endl;
	}
}




// Инициализируем модель, задавая функции u, f и тип сетки
void FEM::init(
	function3D _u,
	function3D _f,
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
	n = nodesCount;

	q1.resize(nodesCount);
	q2.resize(nodesCount);
	q3.resize(nodesCount);

	for (size_t i = 0; i < heigth; i++)
		for (size_t j = 0; j < width; j++)
		{
			int k = i * width + j;
			q3[k] = u(nodes[k].x, nodes[k].y, times[0]);
			q2[k] = u(nodes[k].x, nodes[k].y, times[1]);
			q1[k] = u(nodes[k].x, nodes[k].y, times[2]);
		}

	b3 = buildGlobalVectorb(0);
	b2 = buildGlobalVectorb(1);
	b1 = buildGlobalVectorb(2);

	for (size_t timeLayer = 3; timeLayer < times.size(); timeLayer++) {
		implicitCheme4(timeLayer);

		for (int j = 0; j < A.size(); ++j) {
			for (int i = j + 1; i < A.size(); ++i) {

				double toMult = A[i][j] / A[j][j];

				for (int k = 0; k < A.size(); ++k)
					A[i][k] -= toMult * A[j][k];

				b[i] -= toMult * b[j];
			}
		}


		vector <double> x;
		x.resize(A.size(), 0);

		for (int i = n - 1; i >= 0; --i) {

			double tmp = 0.0;
			for (int j = i + 1; j < A.size(); ++j) {
				tmp += A[i][j] * x[j];
			}
			x[i] = (b[i] - tmp) / A[i][i];
		}
		q = x;

		b3 = b2;
		b2 = b1;
		b1 = b;

		q3 = q2;
		q2 = q1;
		q1 = q;

		//cout << endl << "q:" << endl << q << endl;
		cout << endl << "Residual: " << endl << calcNormAtMainNodes(q, t0);
	}
}


void FEM::solveParabolic()
{
	n = nodesCount;

	q1.resize(nodesCount);
	q2.resize(nodesCount);

	for (size_t i = 0; i < heigth; i++)
		for (size_t j = 0; j < width; j++)
		{
			int k = i * width + j;
			q2[k] = u(nodes[k].x, nodes[k].y, times[0]);
			q1[k] = u(nodes[k].x, nodes[k].y, times[1]);
		}

	b2 = buildGlobalVectorb(0);
	b1 = buildGlobalVectorb(1);

	for (size_t timeLayer = 2; timeLayer < times.size(); timeLayer++) {
		implicitCheme3(timeLayer);
		//implicitCheme3Parabolic(timeLayer);
		//CranckNicolson(timeLayer);

		for (int j = 0; j < A.size(); ++j) {
			for (int i = j + 1; i < A.size(); ++i) {

				double toMult = A[i][j] / A[j][j];

				for (int k = 0; k < A.size(); ++k)
					A[i][k] -= toMult * A[j][k];

				b[i] -= toMult * b[j];
			}
		}


		vector <double> x;
		x.resize(A.size(), 0);

		for (int i = n - 1; i >= 0; --i) {

			double tmp = 0.0;
			for (int j = i + 1; j < A.size(); ++j) {
				tmp += A[i][j] * x[j];
			}
			x[i] = (b[i] - tmp) / A[i][i];
		}
		q = x;

		b2 = b1;
		b1 = b;

		q2 = q1;
		q1 = q;

		//cout << endl << "q:" << endl << q << endl;
		cout << endl << "Residual: " << endl << calcNormAtMainNodes(q, t0);
	}
}




//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Строим глобальную матрицу системы нелинейных уравнений (см. с. 239)
void FEM::CranckNicolson(int timeLayer)
{
	A.clear();
	A.resize(nodesCount);
	for (size_t i = 0; i < nodesCount; i++)
		A[i].resize(nodesCount, 0);

	buildGlobalMatrixG();
	buildGlobalMatrixM();
	b = buildGlobalVectorb(timeLayer);

	t0 = times[timeLayer];
	t1 = times[timeLayer - 1];
	t2 = times[timeLayer - 2];

	d1 = t0 - t2;
	d2 = (t0 * (t0 - t1 - t2) + t1 * t2) / 2;
	m1 = (t0 - t2) / (t1 - t2);
	m2 = (t0 - t1) / (t1 - t2);


	// Собираем левую часть
	for (size_t i = 0; i < nodesCount; i++)
		for (size_t j = 0; j < nodesCount; j++)
			A[i][j] = G[i][j] / 2;

	double tmp = gamma / 2 + sigma / d1 + chi / d2;
	for (size_t i = 0; i < nodesCount; i++)
		for (size_t j = 0; j < nodesCount; j++)
			A[i][j] += M[i][j] * tmp;


	// Собираем правую часть
	b = (b + b2) / 2.0
		- G * q2 / 2.0
		+ M * (q1 * m1 * chi / d2 + q2 * ((-gamma / 2.0) + (sigma / d1) - (m2 * chi / d2)));

	//outputA();
	//cout << endl << "b:" << endl << b << endl;

	// Добавляем краевые условия
	for (size_t i = 0; i < nodesCount; i++)
	{
		if (nodes[i].type == 1) {
			A[i].clear();
			A[i].resize(nodesCount, 0);
			A[i][i] = 1;
			b[i] = u(nodes[i].x, nodes[i].y, t0);
		}
	}


	/*cout << endl << "Added 1st boundary conditions:" << endl;
	cout << endl << "A:" << endl;
	outputA();*/

	//cout << endl << "b:" << endl << b << endl;

}




void FEM::implicitCheme3(int timeLayer)
{
	A.clear();
	A.resize(nodesCount);
	for (size_t i = 0; i < nodesCount; i++)
		A[i].resize(nodesCount, 0);

	buildGlobalMatrixG();
	buildGlobalMatrixM();
	b = buildGlobalVectorb(timeLayer);

	t0 = times[timeLayer];
	t1 = times[timeLayer - 1];
	t2 = times[timeLayer - 2];


	t01 = t0 - t1;
	t02 = t0 - t2;
	t20 = t2 - t0;
	t21 = t2 - t1;
	t10 = t1 - t0;
	t12 = t1 - t2;



	// Собираем левую часть
	double tmp = 2.0 * chi / (t01 * t02)
		+ sigma * (t02 + t01) / (t02 * t01) + gamma;
	for (size_t i = 0; i < nodesCount; i++)
		for (size_t j = 0; j < nodesCount; j++)
			A[i][j] += M[i][j] * tmp + G[i][j];


	// Собираем правую часть
	b = b - (2 * chi / (t01 * t02) + sigma * t01 / (t02 * t12)) * M * q2
		- (2 * chi / (t01 * t12) + sigma * t02 / (t01 * t12)) * M * q1;

	//outputA();
	//cout << endl << "b:" << endl << b << endl;


	// Добавляем краевые условия
	for (size_t i = 0; i < nodesCount; i++)
	{
		if (nodes[i].type == 1) {
			A[i].clear();
			A[i].resize(nodesCount, 0);
			A[i][i] = 1;
			b[i] = u(nodes[i].x, nodes[i].y, t0);
		}
	}


	//cout << endl << "Added 1st boundary conditions:" << endl;
	//cout << endl << "A:" << endl;
	//outputA();

	//cout << endl << "b:" << endl << b << endl;

}


void FEM::implicitCheme3Parabolic(int timeLayer)
{
	A.clear();
	A.resize(nodesCount);
	for (size_t i = 0; i < nodesCount; i++)
		A[i].resize(nodesCount, 0);

	buildGlobalMatrixG();
	buildGlobalMatrixM();
	b = buildGlobalVectorb(timeLayer);

	t0 = times[timeLayer];
	t1 = times[timeLayer - 1];
	t2 = times[timeLayer - 2];


	t01 = t0 - t1;
	t02 = t0 - t2;
	t20 = t2 - t0;
	t21 = t2 - t1;
	t10 = t1 - t0;
	t12 = t1 - t2;



	// Собираем левую часть
	double tmp = sigma * (t02 + t01) / (t02 * t01) + gamma;
	for (size_t i = 0; i < nodesCount; i++)
		for (size_t j = 0; j < nodesCount; j++)
			A[i][j] += M[i][j] * tmp + G[i][j];


	// Собираем правую часть
	b = b
		- (sigma * t01 / (t02 * t12)) * M * q2
		+ (sigma * t02 / (t01 * t12)) * M * q1;

	//outputA();
	//cout << endl << "b:" << endl << b << endl;

	// Добавляем краевые условия
	for (size_t i = 0; i < nodesCount; i++)
	{
		if (nodes[i].type == 1) {
			A[i].clear();
			A[i].resize(nodesCount, 0);
			A[i][i] = 1;
			b[i] = u(nodes[i].x, nodes[i].y, t0);
		}
	}

	//cout << endl << "Added 1st boundary conditions:" << endl;
	//cout << endl << "A:" << endl;
	//outputA();
	//cout << endl << "b:" << endl << b << endl;
}

// Строим глобальную матрицу системы нелинейных уравнений (см. с. 239)
void FEM::implicitCheme4(int timeLayer)
{
	A.clear();
	A.resize(nodesCount);
	for (size_t i = 0; i < nodesCount; i++)
		A[i].resize(nodesCount, 0);

	buildGlobalMatrixG();
	buildGlobalMatrixM();
	b = buildGlobalVectorb(timeLayer);

	t0 = times[timeLayer];
	t1 = times[timeLayer - 1];
	t2 = times[timeLayer - 2];
	t3 = times[timeLayer - 3];

	t01 = t0 - t1;
	t02 = t0 - t2;
	t03 = t0 - t3;
	t20 = t2 - t0;
	t21 = t2 - t1;
	t10 = t1 - t0;
	t12 = t1 - t2;
	t23 = t2 - t3;
	t30 = t3 - t0;
	t31 = t3 - t1;
	t32 = t3 - t2;
	t13 = t1 - t3;



	// Собираем левую часть
	double tmp = 2 * chi * (t01 + t02 + t03) / (t01 * t02 * t03)
				+ sigma * (1.0 / t01 + 1.0 / t02 + 1.0 / t03) + gamma;
	for (size_t i = 0; i < nodesCount; i++)
		for (size_t j = 0; j < nodesCount; j++)
			A[i][j] += M[i][j] * tmp + G[i][j];


	// Собираем правую часть
	b = b
		+ (2 * chi * (t01 + t02) / (t03 * t13 * t23) + sigma * (t01 * t02) / (t03 * t13 * t23)) * M * q3
		- (2 * chi * (t01 + t03) / (t02 * t12 * t23) + sigma * (t01 * t03) / (t02 * t12 * t23)) * M * q2
		+ (2 * chi * (t02 + t03) / (t01 * t12 * t13) + sigma * (t02 * t03) / (t01 * t12 * t13)) * M * q1;

	//outputA();
	//cout << endl << "b:" << endl << b << endl;

	// Добавляем краевые условия
	for (size_t i = 0; i < nodesCount; i++)
	{
		if (nodes[i].type == 1) {
			A[i].clear();
			A[i].resize(nodesCount, 0);
			A[i][i] = 1;
			b[i] = u(nodes[i].x, nodes[i].y, t0);
		}
	}

	//cout << endl << "Added 1st boundary conditions:" << endl;
	//cout << endl << "A:" << endl;
	//outputA();
	//cout << endl << "b:" << endl << b << endl;
}



void FEM::buildGlobalMatrixG()
{
	G.clear();
	G.resize(nodesCount);
	for (size_t i = 0; i < nodesCount; i++)
		G[i].resize(nodesCount, 0);


	for (size_t i = 0; i < heigth - 1; i++)
	{
		for (size_t j = 0; j < width - 1; j++)
		{
			int k = i * width + j;
			buildLocalmatrixG(k);

			vector <int> nearestNodesIndexes = { k, k + 1, k + width, k + width + 1 };
			for (size_t i1 = 0; i1 < 4; i1++)
				for (size_t j1 = 0; j1 < 4; j1++)
					G[nearestNodesIndexes[i1]][nearestNodesIndexes[j1]] += GLocal[i1][j1];
		}
	}
}



void FEM::buildGlobalMatrixM()
{
	M.clear();
	M.resize(nodesCount);
	for (size_t i = 0; i < nodesCount; i++)
		M[i].resize(nodesCount, 0);

	for (size_t i = 0; i < heigth - 1; i++)
	{
		for (size_t j = 0; j < width - 1; j++)
		{
			int k = i * width + j;
			buildLocalmatrixM(k);

			vector <int> nearestNodesIndexes = { k, k + 1, k + width, k + width + 1 };
			for (size_t i1 = 0; i1 < 4; i1++)
				for (size_t j1 = 0; j1 < 4; j1++)
					M[nearestNodesIndexes[i1]][nearestNodesIndexes[j1]] += MLocal[i1][j1];
		}
	}
}




// Строим глобальный вектор правой части системы нелинейных уравнений
vector1D FEM::buildGlobalVectorb(int timeLayer)
{
	t = times[timeLayer];
	tmpVector.clear();
	tmpVector.resize(nodesCount, 0);


	for (size_t i = 0; i < heigth - 1; i++)
	{
		for (size_t j = 0; j < width - 1; j++)
		{
			int k = i * width + j;
			buildLocalVectorb(k);

			vector <int> nearestNodesIndexes = { k, k + 1, k + width, k + width + 1 };
			for (size_t i1 = 0; i1 < 4; i1++)
				tmpVector[nearestNodesIndexes[i1]] += bLocal[i1];
		}
	}

	return tmpVector;
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



// Построение локальной матрицы жёсткости
void FEM::buildLocalmatrixG(int elemNumber)
{
	hx = nodes[elemNumber + 1].x - nodes[elemNumber].x;
	hy = nodes[elemNumber + width].y - nodes[elemNumber].y;

	GLocal = { {0, 0, 0, 0},
				{0, 0, 0, 0},
				{0, 0, 0, 0},
				{0, 0, 0, 0} };

	tmpMatrix = { {2, -2, 1, -1},
				{-2, 2, -1, 1},
				{1, -1, 2, -2},
				{-1, 1, -2, 2} };

	for (size_t i = 0; i < 4; i++)
		for (size_t j = 0; j < 4; j++)
			GLocal[i][j] += tmpMatrix[i][j] * lambda * hy / (6 * hx);

	tmpMatrix = { {2, 1, -2, -1},
				{1, 2, -1, -2},
				{-2, -1, 2, 1},
				{-1, -2, 1, 2} };


	for (size_t i = 0; i < 4; i++)
		for (size_t j = 0; j < 4; j++)
			GLocal[i][j] += tmpMatrix[i][j] * lambda * hx / (6 * hy);
}


// Построение локальной матрицы масс
void FEM::buildLocalmatrixM(int elemNumber)
{
	hx = nodes[elemNumber + 1].x - nodes[elemNumber].x;
	hy = nodes[elemNumber + width].y - nodes[elemNumber].y;

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
	hx = nodes[elemNumber + 1].x - nodes[elemNumber].x;
	hy = nodes[elemNumber + width].y - nodes[elemNumber].y;

	bLocal = { 0, 0, 0, 0 };
	f1 = f(nodes[elemNumber].x, nodes[elemNumber].y, t);
	f2 = f(nodes[elemNumber + 1].x, nodes[elemNumber + 1].y, t);
	f3 = f(nodes[elemNumber + width].x, nodes[elemNumber + width].y, t);
	f4 = f(nodes[elemNumber + width + 1].x, nodes[elemNumber + width + 1].y, t);

	bLocal[0] = (hx*hy / 36) * (4 * f1 + 2 * f2 + 2 * f3 + 1 * f4);
	bLocal[1] = (hx*hy / 36) * (2 * f1 + 4 * f2 + 1 * f3 + 2 * f4);
	bLocal[2] = (hx*hy / 36) * (2 * f1 + 1 * f2 + 4 * f3 + 2 * f4);
	bLocal[3] = (hx*hy / 36) * (1 * f1 + 2 * f2 + 2 * f3 + 4 * f4);
}