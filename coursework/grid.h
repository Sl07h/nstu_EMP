#pragma once
#include "head.h"


struct NODE {

	bool isFirstNode = false;
	int i;
	double x;
	int type = -9000;		// -9000	начение при инициализации
							// -1		фиктивный узел
							// 0		внутренний узел
							// n		номер границы
	int border;				// номер границы

	void setNodesData(double _x, double _y, int _i, int _j, int _type, double _coef) {
		x = _x;
		y = _y;
		i = _i;
		j = _j;
		type = _type;
		if (i % int(pow(2, _coef)) == 0 && j % int(pow(2.0, _coef)) == 0)
			isFirstNode = true;
	}
};


class GRID
{
public:
	void inputGrid();
	void inputTime();
	void buildGrid();
	void buildTimeGrid();
	void showGrid();
	void saveGridAndBorder(const string &filepathGrid, const string &filepathGridBorder);

protected:

	int coefGrid, // Сколько раз дробили сетку по пространству
		coefTime; // Сколько раз дробили сетку по времени

	// Пространство
	int width;
	double xLeft, xRight;
	double hx, nx, kx;
	double dx;

	int heigth;
	double yLower, yUpper;
	double hy, ny, ky;
	double dy;

	bool isGridUniform;
	int nodesCount, finiteElementsCount;

	// Время
	bool isTimeUniform;
	int tCount;
	double tFirst, tLast;
	double ht, nt, kt;
	double dt;

	// Узлы
	vector <NODE> nodes;
	vector1D times;
};