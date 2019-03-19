#include "matrix.h"




// Ввод 5-диагональной матрицы из файла
int matrix::readMatrixFromFile(const const char *filename) {

	std::ifstream fin;
	fin.open(filename);

	fin >> n >> m;
	fin >> E >> maxiter;

	di.resize(n);
	for (int i = 0; i < n; ++i)
		fin >> di[i];

	al1.resize(n - 1);
	for (int i = 0; i < al1.size(); ++i)
		fin >> al1[i];

	al2.resize(n - m - 2);
	for (int i = 0; i < al2.size(); ++i)
		fin >> al2[i];

	al3.resize(n - m - 3);
	for (int i = 0; i < al3.size(); ++i)
		fin >> al3[i];


	au1.resize(n - 1);
	for (int i = 0; i < au1.size(); ++i)
		fin >> au1[i];

	au2.resize(n - m - 2);
	for (int i = 0; i < au2.size(); ++i)
		fin >> au2[i];

	au3.resize(n - m - 3);
	for (int i = 0; i < au3.size(); ++i)
		fin >> au3[i];




	fin.close();
	return 0;
}


// Вывод 5-ми диагональной матрицы в файл
void matrix::writeMatrixToFile(const char * filename) {

	std::ofstream fout;
	fout.open(filename);

	fout << n << " " << m << endl;
	fout << E << " " << maxiter << endl;

	for (int i = 0; i < n; ++i)
		fout << di[i] << " ";
	fout << endl;


	for (int i = 0; i < al1.size(); ++i)
		fout << al1[i] << " ";
	fout << endl;

	for (int i = 0; i < al2.size(); ++i)
		fout << al2[i] << " ";
	fout << endl;

	for (int i = 0; i < al3.size(); ++i)
		fout << al3[i] << " ";
	fout << endl;



	for (int i = 0; i < au1.size(); ++i)
		fout << au1[i] << " ";
	fout << endl;

	for (int i = 0; i < au2.size(); ++i)
		fout << au2[i] << " ";
	fout << endl;

	for (int i = 0; i < au3.size(); ++i)
		fout << au3[i] << " ";
	fout << endl;


	fout.close();
}


// Создаём матрицу A(k)
void matrix::generateMatrixWith7Diagonals(int new_n, int new_m) {

	n = new_n;
	m = new_m;
	di.clear();
	di.resize(n, 0);
	au1.resize(n - 1);
	al1.resize(n - 1);

	au2.resize(n - m - 2);
	al2.resize(n - m - 2);

	au3.resize(n - m - 3);
	al3.resize(n - m - 3);

	for (int i = 0; i < al1.size(); ++i) {
		au1[i] = -rand() % 5;
		al1[i] = -rand() % 5;
	}

	for (int i = 0; i < al2.size(); ++i) {
		au2[i] = -rand() % 5;
		al2[i] = -rand() % 5;
	}

	for (int i = 0; i < al3.size(); ++i) {
		au3[i] = -rand() % 5;
		al3[i] = -rand() % 5;
	}


	for (int i = 0; i < di.size(); ++i)
		di[i] = -calcAii(i);

	di[0]++;
}


// Рассчёт суммы элементов строки при генерации матрицы A(k)
real matrix::calcAii(int i) {

	real sum = 0;

	if (i >= 1) { // Нижний треугольник
		sum += al1[i - 1];
		if (i >= m + 2) {
			sum += al2[i - m - 2];
			if (i >= m + 3)
				sum += al3[i - m - 3];
		}
	}

	sum += di[i];

	if (i < n - 1) { // Верхний треугольник
		sum += au1[i];
		if (i < n - m - 2) {
			sum += au2[i];
			if (i < n - m - 3)
				sum += au3[i];
		}
	}

	return sum;
}


// Меняем знак внедиагональных элементов на противоположный
void matrix::invertSigns() {

	for (int i = 0; i < al1.size(); ++i) {
		au1[i] = abs(au1[i]);
		al1[i] = abs(al1[i]);
	}

	for (int i = 0; i < al2.size(); ++i) {
		au2[i] = abs(au2[i]);
		al2[i] = abs(al2[i]);
	}

	for (int i = 0; i < al3.size(); ++i) {
		au3[i] = abs(au3[i]);
		al3[i] = abs(al3[i]);
	}
}