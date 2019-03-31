#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <functional>

using namespace std;


typedef std::function<double(double, double)> function2D;


// https://en.wikipedia.org/wiki/Numerical_differentiation#Higher_Derivatives

//
//function2D calcSecondDerivativeX(const function2D& f) {
//	return [f](double x, double y) -> double {
//		const double h = 0.0001;
//		return (-f(x + 2 * h, y) + 16 * f(x + h, y) - 30 * f(x, y) + 16 * f(x - h, y) - f(x - 2 * h, y)) / (12 * h*h);
//	};
//}
//
//function2D calcSecondDerivativeY(const function2D& f) {
//	return [f](double x, double y) -> double {
//		const double h = 0.0001;
//		return (-f(x, y + 2 * h) + 16 * f(x, y + h) - 30 * f(x, y) + 16 * f(x, y - h) - f(x, y - 2 * h)) / (12 * h*h);
//	};
//}
//
//function2D calcLaplacian(const function2D& f) {
//	auto fxx = calcSecondDerivativeX(f);
//	auto fyy = calcSecondDerivativeY(f);
//	return [fxx, fyy](double x, double y) -> double {
//		return fxx(x, y) + fyy(x, y);
//	};
//}