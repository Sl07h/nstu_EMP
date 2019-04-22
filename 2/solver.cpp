#include "solver.h"



// Решение СЛАУ
void SOLVER::calcWithLUDecomposition()
{
	// LU разложение
	LUdecomposition();

	// Прямой ход
	executeDirectTraversal();

	// Обратный ход
	executeReverseTraversal();
}



// Проверяем условие выхода
bool SOLVER::shouldCalc(int i)
{
	// Выход по числу итераций
	if (i > maxiter)
		return false;

	//// Выход шагу
	//if (calcNormE(q - qPrevTime) / calcNormE(q) < delta)
	//	return false;

	//// Выход по относительной невязке
	//if (calcNormE(multAonQ() - b) / calcNormE(b) < E)
	//	return false;
	return true;
}



// LU разложение
void SOLVER::LUdecomposition()
{
	// 1 0 0 0   1 2 0 0   1 2 0 0
	// 3 2 0 0 * 0 1 2 0 = 3 8 4 0
	// 0 3 3 0   0 0 1 2   0 3 9 6
	// 0 0 3 4   0 0 0 1   0 0 3 10
	int lIndex = di.size();
	for (size_t i = 1; i < lIndex; i++)
	{
		au[i - 1] = au[i - 1] / di[i - 1];
		di[i] = di[i] - al[i - 1] * au[i - 1];
	}
}



// Прямой ход
// LUq=b, y=Uq
// Ly=b
void SOLVER::executeDirectTraversal()
{
	q[0] = b[0] / di[0];

	for (size_t i = 1; i < di.size(); i++)
		q[i] = (b[i] - al[i - 1] * q[i - 1]) / di[i];

	b = q;
}



// Обратный ход
// LUq=b, y=Uq
// Uq = y
void SOLVER::executeReverseTraversal()
{
	int lIndex = di.size() - 1;
	q[lIndex] = b[lIndex];

	for (int i = lIndex - 1; i >= 0; i--)
		q[i] = (b[i] - q[i + 1] * au[i]);
}


// Проверка решения СЛАУ
void SOLVER::testSLAE()
{
	di = { 1, 8, 9, 10 };
	al = { 3, 3, 3 };
	au = { 2, 4, 6 };
	b = { 5, 31, 57, 49 };
	q.resize(4, 0);
	calcWithLUDecomposition();
	cout << q << endl;
}



// Умножение матрицы A на вектор q
vector1D SOLVER::multAonQ()
{
	vector1D tmp;
	tmp.resize(di.size());

	if (di.size() >= 2)
		tmp[0] = di[0] * q[0] + au[0] * q[1];

	if (di.size() >= 3)
		for (size_t i = 1; i < di.size() - 1; i++)
			tmp[i] = al[i - 1] * q[i - 1] + di[i] * q[i] + au[i] * q[i + 1];

	int lIndex = di.size() - 1;
	tmp[lIndex] = al[lIndex - 1] * q[lIndex - 1] + di[lIndex] * q[lIndex];
	return tmp;
}
