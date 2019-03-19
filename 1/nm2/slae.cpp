#include "slae.h"

void SLAE::convMatrixToSparse()
{

}


// Преобразование 5-ми диагональной матрицы в плотный формат
void SLAE::convMatrixToDense() {

	A.clear();
	A.resize(n);

	for (int i = 0; i < n; ++i) {
		A[i].resize(n, 0);
		A[i][i] = di[i];
	}


	int j = 1;
	for (int i = 0; i < al1.size(); ++i, ++j) {

		A[i][j] = au1[i];
		A[j][i] = al1[i];
	}

	j = m + 2;
	for (int i = 0; i < al2.size(); ++i, ++j) {

		A[i][j] = au2[i];
		A[j][i] = al2[i];
	}

	j = m + 3;
	for (int i = 0; i < al3.size(); ++i, ++j) {

		A[i][j] = au3[i];
		A[j][i] = al3[i];
	}
}


// Вывод плотной матрицы из файла
int SLAE::readDenseMatrixFromFile(const char * filename)
{
	std::ifstream fin;
	fin.open(filename);
	fin >> n;
	F.resize(n);

	A.resize(n);
	for (size_t i = 0; i < n; i++)
	{
		A[i].resize(n);
	}

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			fin >> A[i][j];
		}
	}

	return 0;
}


// Вывод плотной матрицы в файл
void SLAE::writeDenseMatrixToFile(const char * filename) {

	std::ofstream fout;
	fout.open(filename);


	for (int i = 0; i < n;++i) {
		for (int j = 0; j < n; ++j)
			fout << A[i][j] << "\t";
		fout << endl;
	}

	fout.close();
}


// Умножение i-й строки матрицы на вектор
real SLAE::multLine(vector <real> &line, int i, int mode) {

	real sum = 0;
	if (mode == 1 || mode == 3) {	// Нижний треугольник

		if (i > 0) {

			sum += al1[i - 1] * line[i - 1];
			if (i > m + 1) {
				sum += al2[i - m - 2] * line[i - m - 2];
				if (i > m + 2)
					sum += al3[i - m - 3] * line[i - m - 3];
			}
		}
	}


	if (mode == 2 || mode == 3) {	// Главная диагональ
									// и верхний треугольник
		sum += di[i] * line[i];
		if (i < n - 1) {

			sum += au1[i] * line[i + 1];

			if (i < n - m - 2) {
				sum += au2[i] * line[i + m + 2];
				if (i < n - m - 3)
					sum += au3[i] * line[i + m + 3];
			}
		}
	}

	return sum;
}


// Умножение матрицы на вектор
void SLAE::mult() {
	
	int index;
	F.clear();
	F.resize(n, 0);
	// Нижний треугольник
	index = 1;
	for (int i = 0; i < al1.size(); ++i, ++index)
		F[index] += al1[i] * xk[i];
	index = m + 2;
	for (int i = 0; i < al2.size(); ++i, ++index)
		F[index] += al2[i] * xk[i];
	index = m + 3;
	for (int i = 0; i < al3.size(); ++i, ++index)
		F[index] += al3[i] * xk[i];


	// Главная диагональ
	for (int i = 0; i < di.size(); ++i)
		F[i] += di[i] * xk[i];


	// Верхний треугольник
	index = 1;
	for (int i = 0; i < au1.size(); ++i, ++index)
		F[i] += au1[i] * xk[index];
	index = m + 2;
	for (int i = 0; i < au2.size(); ++i, ++index)
		F[i] += au2[i] * xk[index];
	index = m + 3;
	for (int i = 0; i < au3.size(); ++i, ++index)
		F[i] += au3[i] * xk[index];
}


// Метод Якоби. 0 < w < 1
// Используется общая память для xk и xk1
void SLAE::Jacobi(real w) {

	real sum;
	vector <real> xk1;
	xk1.resize(n);

	for (int i = 0; i < n; ++i) {
		sum = multLine(xk, i, 3);
		//xk[i] += w * (F[i] - sum) / di[i];
		xk1[i] = xk[i] + w * (F[i] - sum) / di[i];
	}
	xk = xk1;
}


// Метод Гаусса-Зейделя. 0 < w < 2
void SLAE::GaussSeildel(real w) {

	real sum;
	vector <real> xk1 = xk;

	for (int i = 0; i < n; ++i) {

		sum = multLine(xk1, i, 1);
		sum += multLine(xk, i, 2);

		xk1[i] = xk[i] + w * (F[i] - sum) / di[i];
	}
	xk = xk1;
}


// Решение СЛАУ итерационным методом
// 1  Метод Якоби
// 2  Метод Гаусса-Зейделя
int SLAE::calcIterative(int mode, real w) {
	int i = 0;
	while (i < maxiter && calcRelativeDiscrepancy() >= E) {

		if (mode == 1)
			Jacobi(w);
		else
			GaussSeildel(w);

		++i;
	}


	return i;
}


// Поиск оптимального веса
real SLAE::findOptimalW(int mode, std::ofstream& fout) {

	real optimalW, tmpW;
	int max_i, min_i = maxiter, tmp_i;
	if (mode == 1) max_i = 101;
	else max_i = 200;

	for (int i = 0; i < max_i; ++i) {

		generateInitualGuess(getDimention());
		tmpW = real(i) / 100;
		tmp_i = calcIterative(mode, tmpW);
		if (tmp_i < min_i) {
			min_i = tmp_i;
			optimalW = tmpW;
		}


		if (i % 10 == 0) // Выводим таблицу с точность 0.1
			writeTableToFile(fout, 1, tmpW, tmp_i, calcCondNumber());
	}
	generateInitualGuess(getDimention());
	min_i = calcIterative(mode, optimalW);
	writeTableToFile(fout, 2, optimalW, min_i, calcCondNumber());

	return optimalW;
}


// Вычисление нормы в Евклидовом пространстве
real SLAE::calcNormE(vector <real> &x) {

	real normE = 0;
	for (int i = 0; i < n; i++)
		normE += x[i] * x[i];

	return sqrt(normE);
}


// Рассчёт относительной невязки
real SLAE::calcRelativeDiscrepancy() {

	vector <real> numerator, denominator = F;
	numerator.resize(n);

	mult(); // F = A*xk

	for (int i = 0; i < n; ++i)
		numerator[i] = denominator[i] - F[i]; // F - A*xk

	// || F - A*xk || / || F ||
	real res = calcNormE(numerator) / calcNormE(denominator);
	F = denominator;
	return res;
}


// Рассчёт числа обусловленности
int SLAE::calcCondNumber() {

	vector <real> numerator, denominator;
	numerator.resize(n);
	denominator.resize(n);
	for (int i = 0; i < n; ++i) {
		numerator[i] = xk[i] - (i + 1);	// x - x*
		denominator[i] = i + 1;			// x
	}

	return calcNormE(numerator) / (calcNormE(denominator) * calcRelativeDiscrepancy());
}

