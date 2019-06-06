#include "grid.h"


void GRID::inputGrid()
{
	string filepath;
	if (isGridUniform)
		filepath = "input/uniform_grid.txt";
	else
		filepath = "input/nonuniform_grid.txt";

	std::ifstream fin(filepath);
	fin >> xLeft >> xRight
		>> yLower >> yUpper;

	fin >> width >> heigth;
	if (!isGridUniform) {
		fin >> kx >> ky;
		nx = width - 1;
		ny = heigth - 1;
	}
	fin.close();
}


void GRID::inputTime()
{
	string filepath;
	if (isTimeUniform)
		filepath = "input/uniform_time.txt";
	else
		filepath = "input/nonuniform_time.txt";

	std::ifstream fin(filepath);
	fin >> tFirst >> tLast;
	fin >> tCount;
	if (!isTimeUniform) {
		fin >> kt;
		nt = tCount - 1;
	}
	fin.close();
}


void GRID::buildGrid()
{
	//  xLeft         xRight
	//	  *-------------*
	//    0    width    1
	if (isGridUniform) {
		hx = ((xRight - xLeft) / double(width - 1)) / pow(2, coefGrid);
		hy = ((yUpper - yLower) / double(heigth - 1)) / pow(2, coefGrid);
		if (coefGrid != 0) {
			width = (width - 1) * pow(2, coefGrid) + 1;
			heigth = (heigth - 1) * pow(2, coefGrid) + 1;
		}
	}
	else {
		if (coefGrid != 0) {
			width = (width - 1) * pow(2, coefGrid) + 1;
			heigth = (heigth - 1) * pow(2, coefGrid) + 1;
			nx *= pow(2, coefGrid);
			ny *= pow(2, coefGrid);
			kx *= pow(kx, 1.0 / coefGrid);
			ky *= pow(ky, 1.0 / coefGrid);
		}
		hx = (xRight - xLeft) * (1 - kx) / (1 - pow(kx, nx));
		hy = (yUpper - yLower) * (1 - ky) / (1 - pow(ky, ny));
	}


	nodesCount = width * heigth;
	finiteElementsCount = (width - 1) * (heigth - 1);
	nodes.resize(nodesCount);

	if (isGridUniform) {

		size_t i, j, elem;
		double x, y;

		// Внутренние узлы
		for (size_t j = 1; j < heigth - 1; j++)
		{
			i = 1;
			for (size_t elem = j * width + 1; elem < (j + 1) * width - 1; elem++, i++)
			{
				x = xLeft + hx * i;
				y = yLower + hy * j;
				nodes[elem].setNodesData(x, y, i, j, 0, coefGrid);
				nodes[elem].border = 0;
			}
		}


		// 1 Нижняя граница
		i = 0;
		j = 0;
		y = yLower;
		for (size_t elem = 0; elem < width; elem++, i++)
		{
			x = xLeft + hx * i;
			nodes[elem].setNodesData(x, y, i, j, 1, coefGrid);
			nodes[elem].border = 1;
		}


		// 2 Правая граница
		i = width - 1;
		j = 1;
		x = xRight;
		for (size_t elem = 2 * width - 1; elem < width * heigth; elem += width, j++)
		{
			y = yLower + hy * j;
			nodes[elem].setNodesData(x, y, i, j, 1, coefGrid);
			nodes[elem].border = 2;
		}


		// 3 Верхняя граница
		i = 0;
		j = heigth - 1;
		y = yUpper;
		for (size_t elem = width * j; elem < (j + 1) * width - 1; elem++, i++)
		{
			x = xLeft + hx * i;
			nodes[elem].setNodesData(x, y, i, j, 1, coefGrid);
			nodes[elem].border = 3;
		}


		// 4 Левая граница
		i = 0;
		j = 1;
		x = xLeft;
		for (size_t elem = width; elem < (heigth - 1) * width; elem += width, j++)
		{
			y = yLower + hy * j;
			nodes[elem].setNodesData(x, y, i, j, 1, coefGrid);
			nodes[elem].border = 4;
		}

	}
	else {

		size_t i, j, elem;
		double x, y;

		// Внутренние узлы
		dy = hy * ky;
		y = yLower + hy;
		for (size_t j = 1; j < heigth - 1; j++, dy *= ky)
		{
			i = 1;
			dx = hx * kx;
			x = xLeft + hx;
			for (size_t elem = j * width + 1; elem < (j + 1) * width - 1; elem++, i++, dx *= kx)
			{
				nodes[elem].setNodesData(x, y, i, j, 0, coefGrid);
				nodes[elem].border = 0;
				x += dx;
			}
			y += dy;
		}


		// 1 Нижняя граница
		i = 0;
		dx = hx;
		x = xLeft;

		j = 0;
		y = yLower;
		for (size_t elem = 0; elem < width; elem++, i++, dx *= kx)
		{
			nodes[elem].setNodesData(x, y, i, j, 1, coefGrid);
			nodes[elem].border = 1;
			x += dx;
		}


		// 2 Правая граница
		i = width - 1;
		x = xRight;

		j = 1;
		dy = hy * ky;
		y = yLower + hy;

		for (size_t elem = 2 * width - 1; elem < width * heigth; elem += width, j++, dy *= ky)
		{
			nodes[elem].setNodesData(x, y, i, j, 1, coefGrid);
			nodes[elem].border = 2;
			y += dy;
		}


		// 3 Верхняя граница
		i = 0;
		dx = hx;
		x = xLeft;

		j = heigth - 1;
		y = yUpper;
		for (size_t elem = width * j; elem < (j + 1) * width - 1; elem++, i++, dx *= kx)
		{
			nodes[elem].setNodesData(x, y, i, j, 1, coefGrid);
			nodes[elem].border = 3;
			x += dx;
		}


		// 4 Левая граница
		i = 0;
		x = xLeft;

		j = 1;
		dy = hy * ky;
		y = yLower + hy;

		for (size_t elem = width; elem < (heigth - 1) * width; elem += width, j++, dy *= ky)
		{
			nodes[elem].setNodesData(x, y, i, j, 1, coefGrid);
			nodes[elem].border = 4;
			y += dy;
		}
	}
}


void GRID::buildTimeGrid()
{
	//  tFirst         tLast
	//	  *-------------*
	//    0    tCount   1


	if (isTimeUniform) {

		ht = ((tLast - tFirst) / double(tCount - 1)) / pow(2, coefTime);
		if (coefTime != 0)
			tCount = (tCount - 1) * pow(2, coefTime) + 1;

		times.resize(tCount);
		size_t i, elem;
		double t;
		// Первый элемент
		times[0] = tFirst;
		i = 1;
		for (elem = 1; elem < tCount; elem++, i++)
			times[elem] = tFirst + ht * i;

		// Последний элемент
		times[tCount - 1] = tLast;
	}

	else {

		if (coefTime != 0) {
			tCount = (tCount - 1) * pow(2, coefTime) + 1;
			nt *= pow(2, coefTime);
			kt *= pow(kt, 1.0 / coefTime);
		}
		times.resize(tCount);
		ht = (tLast - tFirst) * (1 - kt) / (1 - pow(kt, nt));
		double t;
		size_t i, elem;
		i = 1;
		dt = ht * kt;
		t = tFirst + ht;
		// Первый элемент
		times[0] = tFirst;
		for (elem = 1; elem < tCount; elem++, i++, dt *= kt)
		{
			times[elem] = t;
			t += dt;
		}
		// Последний элемент
		times[tCount - 1] = tLast;
	}
}


// Отображние сетки на экран
void GRID::showGrid() {

	for (size_t i = 0; i < width; i++)
		cout << nodes[i].x << " ";
}


// Сохранение внутренних и внешних узлов в 2 файлах
void GRID::saveGridAndBorder(const string &filepathGrid, const string &filepathGridBorder) {

	ofstream grid(filepathGrid);
	ofstream border(filepathGridBorder);
	for (size_t i = 0; i < nodesCount; i++)
	{
		if (nodes[i].type > 0)
			border << nodes[i].x << " " << nodes[i].y << endl;
		else
			grid << nodes[i].x << " " << nodes[i].y << endl;
	}

	border.close();
	grid.close();
}