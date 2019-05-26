#include "grid.h"


void GRID::inputGrid()
{
	string filepath;
	if (isGridUniform)
		filepath = "input/uniform_grid.txt";
	else
		filepath = "input/nonuniform_grid.txt";

	std::ifstream fin(filepath);
	fin >> xLeft >> xRight;
	fin >> width;
	if (!isGridUniform) {
		fin >> kx;
		nx = width - 1;
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
		if (coefGrid != 0)
			width = (width - 1) * pow(2, coefGrid) + 1;
	}
	else {
		if (coefGrid != 0) {
			width = (width - 1) * pow(2, coefGrid) + 1;
			nx *= pow(2, coefGrid);
			kx *= pow(kx, 1.0 / coefGrid);
		}
		hx = (xRight - xLeft) * (1 - kx) / (1 - pow(kx, nx));
	}


	nodesCount = width;
	finiteElementsCount = nodesCount - 1;
	nodes.resize(width);

	if (isGridUniform) {

		size_t i, elem;
		double x;

		// Первый элемент
		nodes[0].setNodesData(xLeft, 0, 1, coefGrid);
		i = 1;
		for (elem = 1; elem < nodesCount - 1; elem++, i++)
		{
			x = xLeft + hx * i;
			nodes[elem].setNodesData(x, i, 0, coefGrid);
			nodes[elem].border = 0;
		}
		// Последний элемент
		nodes[nodesCount - 1].setNodesData(xRight, width, 1, coefGrid);

	}
	else {

		double x;
		size_t i, elem;

		i = 1;
		dx = hx * kx;
		x = xLeft + hx;
		// Первый элемент
		nodes[0].setNodesData(xLeft, 0, 1, coefGrid);
		for (elem = 1; elem < width; elem++, i++, dx *= kx)
		{
			nodes[elem].setNodesData(x, i, 0, coefGrid);
			nodes[elem].border = 0;
			x += dx;
		}
		// Последний элемент
		nodes[nodesCount - 1].setNodesData(xRight, width, 1, coefGrid);
	}
}


void GRID::buildTimeGrid()
{
	//  tFirst         tLast
	//	  *-------------*
	//    0    tCount   10
	times.resize(tCount);

	if (isTimeUniform) {

		ht = ((tLast - tFirst) / double(tCount - 1)) / pow(2, coefTime);
		if (coefTime != 0)
			width = (width - 1) * pow(2, coefTime) + 1;

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
			width = (width - 1) * pow(2, coefTime) + 1;
			nt *= pow(2, coefTime);
			kt *= pow(kt, 1.0 / coefTime);
		}

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
		if (nodes[i].type > 0)
			border << nodes[i].x << endl;
		else
			grid << nodes[i].x << endl;

	border.close();
	grid.close();
}