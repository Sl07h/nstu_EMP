#include "fem.h"


// Строим глобальную матрицу системы нелинейных уравнений
void FEM::buildGlobalMatrixA()
{
	for (size_t elemNumber = 0; elemNumber < elemCount; elemNumber++)
	{
		buildLocalmatrixA(elemNumber);
		di[elemNumber] = ALocal[0][0];
		di[elemNumber + 1] = ALocal[0][1];
		al[elemNumber] = ALocal[1][0];
		au[elemNumber] = ALocal[0][1];
	}
}


// Строим глобальный вектор правой части системы нелинейных уравнений
void FEM::buildGlobalVectorb()
{
	for (size_t elemNumber = 0; elemNumber < elemCount; elemNumber++)
	{
		buildLocalVectorb(elemNumber);
		b[elemNumber] = bLocal[elemNumber];
		b[elemNumber + 1] = bLocal[elemNumber + 1];
	}
}






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
	bLocal[0] = f(elemNumber) * hx / 2 +
				sigma * hx / (6 * dt) * (2 * qPrev[0] + qPrev[1]);
	bLocal[1] = f(elemNumber + 1) * hx / 2 +
				sigma * hx / (6 * dt) * (qPrev[0] + 2 * qPrev[1]);
}
