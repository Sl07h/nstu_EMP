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

	// Выход шагу
	if (calcNormE(q - qPrev) / calcNormE(q) < delta)
		return false;

	// Выход по относительной невязке
	if (calcNormE(multAonQ() - b) / calcNormE(b) < E)
		return false;
}


// LU разложение
void SOLVER::LUdecomposition()
{
	for (size_t i = 1; i < di.size(); i++)
	{
		di[i] = di[i] - al[i - 1] * au[i - 1];
		au[i - 1] = au[i - 1] / di[i];
	}
}



// Прямой ход
// LUq=b, y=Uq
// Ly=b
void SOLVER::executeDirectTraversal()
{
	q[0] = b[0] / al[0];

	for (size_t i = 1; i < di.size(); i++)
		q[i] = (b[i] - al[i - 1] * q[i - 1]) / di[i];
}



// Обратный ход
// LUq=b, y=Uq
// Uq = y
void SOLVER::executeReverseTraversal()
{
	for (size_t i = di.size() - 1; i > 0; i--)
		q[i] -= q[i - 1] * au[i - 1];
}

vector1D SOLVER::multAonQ()
{
	if (di.size() >= 2)
		b[0] = di[0] * q[0] + au[0] * q[1];

	if (di.size() >= 3)
		for (size_t i = 1; i < di.size() - 1; i++)
			b[i] = al[i - 1] * q[i - 1] + di[i] * q[i] + au[i] * q[i + 1];

	int lIndex = di.size() - 1;
	b[lIndex] = al[lIndex - 1] * q[lIndex - 1] + di[lIndex] * q[lIndex];
	return vector1D();
}
