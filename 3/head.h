#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <functional>
#include <cmath>

using namespace std;


typedef std::function<double(double)> function1D;
typedef std::function<double(double, double)> function2D;

typedef vector <double> vector1D;
typedef vector <vector <double>> matrix2D;


// Сравнение векторов
inline bool operator==(const vector1D& a, const vector1D& b) {
#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
#endif
	for (int i = 0; i < a.size(); ++i)
		if (a[i] != b[i])
			return false;

	return true;
}

// Сложение векторов
inline vector1D operator+(const vector1D& a, const vector1D& b) {
#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
#endif
	vector1D result = a;
	for (int i = 0; i < b.size(); i++)
		result[i] += b[i];
	return result;
}
// Сложение матриц
inline matrix2D operator+(const matrix2D& a, const matrix2D& b) {
#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
#endif
	matrix2D result = a;
	for (int i = 0; i < b.size(); i++)
		for (int j = 0; j < b.size(); j++)
			result[i][j] += b[i][j];
	return result;
}


// Деление матрицы на число
inline matrix2D operator/(const matrix2D& a, const double& b) {

	matrix2D result = a;
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.size(); j++)
			result[i][j] /= b;
	return result;
}


// Вычитание векторов
inline vector1D operator-(const vector1D& a, const vector1D& b) {
#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
#endif
	vector1D result = a;
	for (int i = 0; i < b.size(); i++)
		result[i] -= b[i];
	return result;
}
// Обратный знак вектора
inline vector1D operator-(const vector1D& a) {
	vector1D result = a;
	for (int i = 0; i < a.size(); i++)
		result[i] = -result[i];
	return result;
}



// Умножение матрицы на вектор
inline vector1D operator*(const matrix2D& a, const vector1D& b) {
	vector1D result = { 0.0, 0.0 };
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.size(); j++)
			result[i] += a[i][j] * b[j];
	return result;
}



// Умножение на число
inline vector1D operator*(const vector1D& a, double b) {
	vector1D result = a;
	for (int i = 0; i < result.size(); i++)
		result[i] *= b;
	return result;
}
// Умножение на число
inline vector1D operator*(double b, const vector1D& a) {
	return operator*(a, b);
}



// Деление на число
inline vector1D operator/(const vector1D& a, double b) {
	vector1D result = a;
	for (int i = 0; i < result.size(); i++)
		result[i] /= b;
	return result;
}
// Деление на число
inline vector1D operator/(double b, const vector1D& a) {
	return operator/(a, b);
}



// Скалярное произведение
inline double operator*(const vector1D& a, const vector1D& b) {
#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
#endif
	double sum = 0;
	for (int i = 0; i < a.size(); i++)
		sum += a[i] * b[i];
	return sum;
}


// Потоковый вывод вектора
inline std::ostream& operator<<(std::ostream& out, const vector1D& v) {
	for (int i = 0; i < v.size() - 1; ++i)
		out << v[i] << ", ";
	out << v.back();
	return out;
}
// Потоковый вывод матрицы
inline std::ostream& operator<<(std::ostream& out, const matrix2D& v) {
	for (int i = 0; i < v.size() - 1; ++i)
		out << v[i] << "  ";
	out << v.back();
	return out;
}


// Потоковый вывод вектора для TeX
inline void printTeXVector(std::ofstream &fout, const vector1D &v, int coefGrid) {
	fout << "$(";
	for (int i = 0; i < v.size() - 1; ++i)
		if (i % int(pow(2, coefGrid)) == 0)
			fout << v[i] << ", ";
	fout << v.back() << ")^T$";
}

