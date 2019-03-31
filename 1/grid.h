#pragma once
#include "head.h"




struct NODE {

	int i, j;
	double x, y;
	double value;			// значение u
	int type = -9000;		// -9000	начение при инициализации
							// -1		фиктивный узел
							// 0		внутренний узел
							// n		номер границы
	int border;				// номер границы
	void setNodesData(double _x, double _y, int _i, int _j, int _type) {
		x = _x;
		y = _y;
		i = _i;
		j = _j;
		type = _type;
	}
};



class GRID
{
public:
	void inputGrid();
	void buildGrid();

protected:
	bool isGridUniform;
	int condType;
	int width, heigth, widthLeft, widthRight, heigthLower, heigthUpper, elemCount;
	double xLeft, xRight, yLower, yUpper;
	double hx, hy;                            // Приращение по осям ox и oy у равномерной сетки
	double stepX0, stepY0;					// x_0, y_0
	double dx, dy;							// x_n = x_0 + dx*(n-1), y_n = y_0 + dy*(n-1) арифметическая прогрессия
	vector <NODE> nodes;
};