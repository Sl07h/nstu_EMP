#pragma once
#include "head.h"


struct NODE {

	bool isFirstNode = false;
	int i, j;
	double x, y;
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
	void buildGrid();
	void showGrid();
	void saveGridAndBorder(const string &filepathGrid, const string &filepathGridBorder);

protected:
	bool isGridUniform;
	int condType, coef;
	int width, heigth, widthLeft, widthRight, heigthLower, heigthUpper, elemCount;
	double xLeft, xRight, yLower, yUpper;
	double hx, hy, nx, ny, kx, ky, hxPrev, hyPrev;
	double dx, dy;	// hx_n = hx_1 + kx^(n-1) геометрическая прогрессия
	vector <NODE> nodes;
};