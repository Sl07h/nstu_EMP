#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "nm2/slae.h"

using namespace std;

typedef double real;


// Класс МКР. L-образная область.
class FDM
{
public:
	void inputData(const char *filename);
	void outputSLAE(const char *fileA, const char *fileB);
	void transformToSLAE();

private:
	int width, heigth, widthLeft, widthRight, heigthLower, heigthUpper, elemCount;

	real lambda, gamma;                     // коэффициенты диффуров
	real xLeft, xRight, yLower, yUpper;
	real hx, hy;                            // Приращение по осям ox и oy

	real E = 1e-15, w = 1.0;                // параметры СЛАУ
	int maxiter;

	vector <vector <real>> A;				// матрица в плотном формате
	vector<real> di, au1, au2, al1, al2;    // 5-ти диагональная матрица
	vector<real> b;                         // вектор правой части СЛАУ
	vector<real> x;                         // вектор решения


	real f(real x, real y);                 // функция правой части
	real u(real x, real y);                 // решение нашего ДУ
	real calcNormE(vector<real> &v);        // вычисление Евклидовой нормы вектора
};


// Функция правой части
real FDM::f(real x, real y)
{
	real f;
	// Test 1,2,3
	//f = 2.0 * x + 2.0 * y;
	// Test 4
	//f = 2.0 * pow(x,2) + 2.0 * pow(y,2) - 4.0;
	// Test 5
	//f = 2.0 * pow(x,3) + 2.0 * pow(y,3) - 6.0 * x - 6.0 * y;
	// Test 6
	//f = 2.0 * pow(x,4) + 2.0 * pow(y,4) - 12.0 * pow(x,2) - 12.0 * pow(y,2);
	// Test 7
	f = 9.0 * cos(2.0 * x + 2.0 * y);
	return f;
}


// Задаём искомую функцию  u
real FDM::u(real i, real j)
{
	real u;
	real x = xLeft + (xRight - xLeft) * i / width;
	real y = yLower + (yUpper - yLower) * j / heigth;

	u = x + y;

	return u;
}


// Евклидова норма
real FDM::calcNormE(vector<real> &v)
{
	real res = 0;
	for (int i = 0; i < v.size(); i++)
		res += v[i] * v[i];

	return sqrt(res);
}


// Ввод данных из файла
void FDM::inputData(const char *filename) {

	std::ifstream fin(filename);

	fin >> E >> w;
	fin >> lambda >> gamma;

	fin >> xLeft >> xRight >> yLower >> yUpper;
	fin >> width >> heigth >>
		widthLeft >> widthRight >>
		heigthLower >> heigthUpper;

	hx = (xRight - xLeft) / real(width);
	hy = (yUpper - yLower) / real(heigth);

	elemCount = width * heigth;

	A.resize(elemCount);
	for (size_t i = 0; i < elemCount; i++)
	{
		A[i].resize(elemCount);
	}
	di.resize(elemCount);
	al1.resize(elemCount);
	al2.resize(elemCount);
	au1.resize(elemCount);
	au2.resize(elemCount);
	b.resize(elemCount);
	x.resize(elemCount);
}


// Вывод слау в файлы
void FDM::outputSLAE(const char * fileA, const char * fileB)
{
	ofstream foutA(fileA), foutB(fileB);

	//foutA << fixed << setprecision(2);
	foutA << elemCount << " " << width << endl;
	foutA << E << " " << maxiter << endl;

	for (size_t i = 0; i < elemCount; i++)
	{
		for (size_t j = 0; j < elemCount; j++)
			foutA << A[i][j] << "\t";

		foutA << endl;
	}

	for (size_t i = 0; i < elemCount; i++)
	{
		foutB << b[i] << endl;
	}
}


//      c    d          Where:
//    -----====         a - height
//  | 66665xxxx !       b - width
//  | 10005xxxx ! e     c - widthLeft
// a| 10005xxxx !       d - widthRight
//  | 100004443 |       e - heightUpper
//  | 100000003 | f     f - heightLower
//  | 100000003 |
//  | 122222222 |
//    ---------
//        b
//
// 66665xxxx
// 10005xxxx
// 10005xxxx
// 100004443
// 100000003
// 100000003
// 122222222
//
// Перевод конечно-разностной схемы к СЛАУ
void FDM::transformToSLAE()
{
	size_t i, j, elem;
	// Обходим все внутренние элементы нижней части "L" по пяти точкам
	for (j = 1; j < heigthLower - 1; j++)
	{
		i = 1;
		for (elem = j * width + 1; elem < (j + 1) * width - 1; elem++, i++)
		{
			A[elem][elem] = -2 / (hx * hx) - 2 / (hy * hy);
			A[elem][elem - 1] = 1 / (hx * hx);
			A[elem][elem + 1] = 1 / (hx * hx);
			A[elem][elem - width] = 1 / (hy * hy);
			A[elem][elem + width] = 1 / (hy * hy);
			b[elem] = u(i, j);
			cout << elem << "\t";
		}
		cout << endl;
	}

	// Обходим все внутренние элементы средней части части "L" по пяти точкам
	i = 1;
	j = heigthLower - 1;
	for (elem = j * width + 1; elem < j * width + widthLeft; elem++)
	{
		A[elem][elem] = -2 / (hx * hx) - 2 / (hy * hy);
		A[elem][elem - 1] = 1 / (hx * hx);
		A[elem][elem + 1] = 1 / (hx * hx);
		A[elem][elem - width] = 1 / (hy * hy);
		A[elem][elem + width] = 1 / (hy * hy);
		b[elem] = u(i, j);
		cout << elem << "\t";
	}
	cout << endl;


	// Обходим все внутренние элементы верхней части "L" по пяти точкам
	for (j = heigthLower; j < heigth - 1; j++)
	{
		i = 1;
		for (elem = j * width + 1; elem < j * width + widthLeft - 1; elem++, i++)
		{
			A[elem][elem] = -2 / (hx * hx) - 2 / (hy * hy);
			A[elem][elem - 1] = 1 / (hx * hx);
			A[elem][elem + 1] = 1 / (hx * hx);
			A[elem][elem - width] = 1 / (hy * hy);
			A[elem][elem + width] = 1 / (hy * hy);
			b[elem] = u(i, j);
			cout << elem << "\t";
		}
		cout << endl;
	}

	// 1
	i = 0;
	j = 0;
	for (elem = 0; elem < heigth - 1; elem += width, j++)
	{
		A[elem][elem] = 1;
		b[elem] = u(i, j);
	}

	// 2
	i = 1;
	j = 0;
	for (elem = 1; elem < width; elem++, i++)
	{
		A[elem][elem] = 1;
		b[elem] = u(i, j);
	}

	// 3
	i = width - 1;
	j = 1;
	for (elem = 2 * width - 1; elem < (width * heigthLower); elem += width, j++)
	{
		A[elem][elem] = 1;
		b[elem] = u(i, j);
	}

	// 4
	i = widthLeft;
	j = heigthLower - 1;
	for (elem = width * heigthLower - widthRight; elem < width * heigthLower - 1; elem++, i++)
	{
		A[elem][elem] = 1;
		b[elem] = u(i, j);
	}

	// 5
	i = widthLeft - 1;
	j = heigthLower;
	for (elem = width * j + widthLeft - 1; elem < elemCount; elem += width, j++)
	{
		A[elem][elem] = 1;
		b[elem] = u(i, j);
	}

	// 6
	i = 0;
	j = heigth - 1;
	for (elem = width * j; elem < width * j + heigthLower - 1; elem++, i++)
	{
		A[elem][elem] = 1;
		b[elem] = u(i, j);
	}

	// fictitious nodes
	for (j = heigthLower; j < heigth; j++)
	{
		i = widthLeft;
		for (elem = widthLeft; elem < width; elem++)
		{
			A[elem][elem] = 1;
			b[elem] = 0;
		}
	}
}



int main()
{
	SLAE slae;
	FDM fdm;
	const char *fileA = "report/A.txt", *fileF = "report/F.txt";
	ofstream fout_res("report/table.txt");

	fdm.inputData("input/input.txt");
	fdm.transformToSLAE();
	fdm.outputSLAE(fileA, fileF);

	slae.readMatrixFromFile(fileA);
	slae.readFFromFile(fileF);
	//slae.convMatrixToSparse();
	//cout << slae.findOptimalW(2, fout_res) << endl;

	return 0;
}