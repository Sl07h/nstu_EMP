#pragma once
#include "slae.h"

/////////////////////////////////////////////////////////////////////////////////////////////
// matrix

// Рассчёт суммы элементов строки при генерации матрицы A(k)
double SLAE::calcAii(int i) {

	double sum = 0;

	if (i >= 1) { // Нижний треугольник
		sum += al1[i - 1];
		if (i >= m) {
			sum += al2[i - m];
		}
	}

	sum += di[i];

	if (i < n - 1) { // Верхний треугольник
		sum += au1[i];
		if (i < n - m) {
			sum += au2[i];
		}
	}

	return sum;
}


/////////////////////////////////////////////////////////////////////////////////////////////
// vect


void SLAE::inputSLAEParameters()
{
	ifstream fin("input/SLAE_parameters.txt");

	fin >> E >> maxiter;

	fin.close();
}

// Вывод вектора X в файл 
void SLAE::writeXToFile(const string &filename) {

	std::ofstream fout;
	fout.open(filename);

	for (int i = 0; i < x.size(); ++i)
		fout << x[i] << endl;
	fout.close();
}

void SLAE::writebToFile(const string & filename) {

	std::ofstream fout;
	fout.open(filename);

	for (int i = 0; i < x.size(); ++i)
		fout << b[i] << endl;
	fout.close();
}

/////////////////////////////////////////////////////////////////////////////////////////////
// slae



// Преобразование 5-ми диагональной матрицы в плотный формат
void SLAE::convMatrixToDense() {

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

	j = m;
	for (int i = 0; i < al2.size(); ++i, ++j) {

		A[i][j] = au2[i];
		A[j][i] = al2[i];
	}
}


// Вывод плотной матрицы в файл
void SLAE::writeDenseMatrixToFile(const string & filename) {

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
double SLAE::multLine(vector <double> &line, int i, int mode) {

	double sum = 0;
	if (mode == 1 || mode == 3) {	// Нижний треугольник

		if (i > 0) {

			sum += al1[i - 1] * line[i - 1];
			if (i > m) {
				sum += al2[i - m] * line[i - m];
			}
		}
	}


	if (mode == 2 || mode == 3) {	// Главная диагональ
									// и верхний треугольник
		sum += di[i] * line[i];
		if (i < n - 1) {

			sum += au1[i] * line[i + 1];

			if (i < n - m) {
				sum += au2[i] * line[i + m];
			}
		}
	}

	return sum;
}


// Умножение матрицы на вектор
void SLAE::mult() {
	
	int index;
	b.clear();
	b.resize(n, 0);
	// Нижний треугольник
	index = 1;
	for (int i = 0; i < al1.size(); ++i, ++index)
		b[index] += al1[i] * x[i];
	index = m;
	for (int i = 0; i < al2.size(); ++i, ++index)
		b[index] += al2[i] * x[i];


	// Главная диагональ
	for (int i = 0; i < di.size(); ++i)
		b[i] += di[i] * x[i];


	// Верхний треугольник
	index = 1;
	for (int i = 0; i < au1.size(); ++i, ++index)
		b[i] += au1[i] * x[index];
	index = m;
	for (int i = 0; i < au2.size(); ++i, ++index)
		b[i] += au2[i] * x[index];
}


// Метод Якоби. 0 < w < 1
// Используется общая память для x и x1
void SLAE::Jacobi(double w) {

	double sum;
	vector <double> x1;
	x1.resize(n);

	for (int i = 0; i < n; ++i) {
		sum = multLine(x, i, 3);
		//x[i] += w * (b[i] - sum) / di[i];
		x1[i] = x[i] + w * (b[i] - sum) / di[i];
	}
	x = x1;
}


// Метод Гаусса-Зейделя. 0 < w < 2
void SLAE::GaussSeildel(double w) {

	double sum;
	vector <double> x1 = x;

	for (int i = 0; i < n; ++i) {

		sum = multLine(x1, i, 1);
		sum += multLine(x, i, 2);

		x1[i] = x[i] + w * (b[i] - sum) / di[i];
	}
	x = x1;
}


// Решение СЛАУ итерационным методом
// 1  Метод Якоби
// 2  Метод Гаусса-Зейделя
int SLAE::calcIterative(int mode, double w) {
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
double SLAE::findOptimalW(int mode) {

	double optimalW = 0.0, tmpW;
	int max_i, min_i = maxiter, tmp_i;
	if (mode == 1) max_i = 101;
	else max_i = 200;

	for (int i = 0; i < max_i; ++i) {

		generateInitualGuess();
		tmpW = double(i) / 100;
		tmp_i = calcIterative(mode, tmpW);
		if (tmp_i < min_i) {
			min_i = tmp_i;
			optimalW = tmpW;
		}


	}
	generateInitualGuess();
	min_i = calcIterative(mode, optimalW);
	
	return optimalW;
}


// Вычисление нормы в Евклидовом пространстве
double SLAE::calcNormE(vector <double> &x) {

	double normE = 0;
	for (int i = 0; i < n; i++)
		normE += x[i] * x[i];

	return sqrt(normE);
}


// Рассчёт относительной невязки
double SLAE::calcRelativeDiscrepancy() {

	vector <double> numerator, denominator = b;
	numerator.resize(n);

	mult(); // b = A*x

	for (int i = 0; i < n; ++i)
		numerator[i] = denominator[i] - b[i]; // b - A*x

	// || b - A*x || / || b ||
	double res = calcNormE(numerator) / calcNormE(denominator);
	b = denominator;
	return res;
}