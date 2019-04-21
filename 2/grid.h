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
	double time;

	void setNodesData(double _x, int _i, int _type, double _coef) {
		x = _x;
		i = _i;
		type = _type;
		if (i % int(pow(2, _coef)) == 0)
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
	bool isGridUniform;
	int width;
	double xLeft, xRight;
	double hx, nx, kx, hxPrev;
	double dx;
	int nodesCount, finiteElementsCount;
	int condType;

	// Время
	bool isTimeUniform;
	int tCount;
	double tFirst, tLast;
	double ht, nt, kt, htPrev;
	double dt = 1;

	// Узлы
	vector <NODE> nodes;
	vector1D times;
};