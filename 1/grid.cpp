#include "grid.h"
#include "fdm.h"



void GRID::inputGrid()
{
	string filepath;
	if (isGridUniform)
		filepath = "input/uniform_grid.txt";
	else
		filepath = "input/nonuniform_grid.txt";

	std::ifstream fin(filepath);

	fin >> xLeft >> xRight >>
		yLower >> yUpper;

	fin >> width >> heigth >>
		widthLeft >> widthRight >>
		heigthLower >> heigthUpper;

	if (!isGridUniform) {
		fin >> kx >> ky;
		nx = width - 1;
		ny = heigth - 1;
	}

	fin.close();
}




void GRID::buildGrid()
{
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

	if (isGridUniform) {
		hx = ((xRight - xLeft) / double(width - 1)) / pow(2, coef);
		hy = ((yUpper - yLower) / double(heigth - 1)) / pow(2, coef);
		cout << hx << "\t" << hy << endl;
		if (coef != 0) {
			width = (width - 1) * pow(2, coef) + 1;
			heigth = (heigth - 1) * pow(2, coef) + 1;
			widthLeft = (widthLeft)* pow(2, coef);
			widthRight = (widthRight - 1) * pow(2, coef) + 1;
			heigthLower = (heigthLower)* pow(2, coef);
			heigthUpper = (heigthUpper - 1) * pow(2, coef) + 1;
		}
	}
	else {
		if (coef != 0) {
			width = (width - 1) * pow(2, coef) + 1;
			heigth = (heigth - 1) * pow(2, coef) + 1;
			widthLeft = (widthLeft)* pow(2, coef);
			widthRight = (widthRight - 1) * pow(2, coef) + 1;
			heigthLower = (heigthLower)* pow(2, coef);
			heigthUpper = (heigthUpper - 1) * pow(2, coef) + 1;
			nx *= pow(2, coef);
			ny *= pow(2, coef);
			kx *= pow(kx, 1.0 / coef);
			ky *= pow(ky, 1.0 / coef);
		}

		hx = (xRight - xLeft) * (1 - kx) / (1 - pow(kx, nx));
		hy = (yUpper - yLower) * (1 - ky) / (1 - pow(ky, ny));
	}




	elemCount = width * heigth;
	cout << "Grid is uniform" << endl
		<< "Coef: " << coef << endl
		<< "Width: " << width << endl
		<< "Height: " << heigth << endl
		<< "Count of elements: " << elemCount << endl
		<< "hx: " << hx << endl
		<< "hy: " << hy << endl;
	nodes.resize(elemCount);


	if (isGridUniform) {


		size_t i, j, elem;
		double x, y, xPrev, yPrev;
		// Обходим все внутренние элементы нижней части "L" по пяти точкам
		for (j = 1; j < heigthLower - 1; j++)
		{
			i = 1;
			for (elem = j * width + 1; elem < (j + 1) * width - 1; elem++, i++)
			{
				x = xLeft + hx * i;
				y = yLower + hy * j;
				nodes[elem].setNodesData(x, y, i, j, 0, coef);
			}
		}

		// Обходим все внутренние элементы средней части части "L" по пяти точкам
		i = 1;
		j = heigthLower - 1;
		y = yLower + hy * j;
		for (elem = j * width + 1; elem < j * width + widthLeft; elem++, i++)
		{
			x = xLeft + hx * i;
			nodes[elem].setNodesData(x, y, i, j, 0, coef);
		}


		// Обходим все внутренние элементы верхней части "L" по пяти точкам
		for (j = heigthLower; j < heigth - 1; j++)
		{
			i = 1;
			for (elem = j * width + 1; elem < j * width + widthLeft - 1; elem++, i++)
			{
				x = xLeft + hx * i;
				y = yLower + hy * j;
				nodes[elem].setNodesData(x, y, i, j, 0, coef);
			}
		}

		// 1
		i = 0;
		j = 0;
		x = xLeft + hx * i;
		for (elem = 0; elem < (heigth - 1) * width; elem += width, j++)
		{
			y = yLower + hy * j;
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 1;
		}

		// 2
		i = 1;
		j = 0;
		y = yLower + hy * j;
		for (elem = 1; elem < width; elem++, i++)
		{
			x = xLeft + hx * i;
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 2;
		}

		// 3
		i = width - 1;
		j = 1;
		x = xLeft + hx * i;
		for (elem = 2 * width - 1; elem < (width * heigthLower); elem += width, j++)
		{
			y = yLower + hy * j;
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 3;
		}

		// 4
		i = widthLeft;
		j = heigthLower - 1;
		y = yLower + hy * j;
		for (elem = width * heigthLower - widthRight; elem < width * heigthLower - 1; elem++, i++)
		{
			x = xLeft + hx * i;
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 4;
		}

		// 5
		i = widthLeft - 1;
		j = heigthLower;
		x = xLeft + hx * i;
		for (elem = width * j + widthLeft - 1; elem < elemCount; elem += width, j++)
		{
			y = yLower + hy * j;
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 5;
		}

		// 6
		i = 0;
		j = heigth - 1;
		y = yLower + hy * j;
		for (elem = width * j; elem < width * j + widthLeft - 1; elem++, i++)
		{
			x = xLeft + hx * i;
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 6;
		}

		// fictitious nodes
		for (j = heigthLower; j < heigth; j++)
		{
			i = widthLeft;
			for (elem = width * j + widthLeft; elem < width * (j + 1); elem++, i++)
			{
				x = xLeft + hx * i;
				y = yLower + hy * j;
				nodes[elem].setNodesData(x, y, i, j, -1, coef);
			}
		}
	}

	else {

		double x, y;
		size_t i, j, elem;


		// Обходим все внутренние элементы нижней части "L" по пяти точкам
		dy = hy * ky;
		y = yLower + hy;
		for (j = 1; j < heigthLower - 1; j++, dy *= ky)
		{
			i = 1;
			dx = hx * kx;
			x = xLeft + hx;
			for (elem = j * width + 1; elem < (j + 1) * width - 1; elem++, i++, dx *= kx)
			{
				nodes[elem].setNodesData(x, y, i, j, 0, coef);
				x += dx;
			}
			y += dy;
		}



		// Обходим все внутренние элементы средней части части "L" по пяти точкам
		i = 1;
		dx = hx * kx;
		x = xLeft + hx;

		j = heigthLower - 1;
		dy = hy * pow(ky, j);
		y = yLower + hy * (1 - pow(ky, j)) / (1 - ky);

		for (elem = j * width + 1; elem < j * width + widthLeft; elem++, i++, dx *= kx)
		{
			nodes[elem].setNodesData(x, y, i, j, 0, coef);
			x += dx;
		}


		// Обходим все внутренние элементы верхней части "L" по пяти точкам
		dy = hy * pow(ky, heigthLower);
		y = yLower + hy * (1 - pow(ky, heigthLower)) / (1 - ky);
		for (j = heigthLower; j < heigth - 1; j++, dy *= ky)
		{
			i = 1;
			dx = hx * kx;
			x = xLeft + hx;

			for (elem = j * width + 1; elem < j * width + widthLeft - 1; elem++, i++, dx *= kx)
			{
				nodes[elem].setNodesData(x, y, i, j, 0, coef);
				x += dx;
			}
			y += dy;
		}



		// 1
		i = 0;
		dx = hx;
		x = xLeft;

		j = 0;
		dy = hy;
		y = yLower;

		for (elem = 0; elem < (heigth - 1) * width; elem += width, j++, dy *= ky)
		{
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 1;
			y += dy;
		}



		// 2
		i = 1;
		dx = hx * kx;
		x = xLeft + hx;

		j = 0;
		dy = hy;
		y = yLower;

		for (elem = 1; elem < width; elem++, i++, dx *= kx)
		{
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 2;
			x += dx;
		}



		// 3
		i = width - 1;
		dx = hx * pow(kx, i - 1);
		x = xLeft + hx * (1 - pow(kx, i)) / (1 - kx);

		j = 1;
		dy = hy * ky;
		y = yLower + hy;

		for (elem = 2 * width - 1; elem < (width * heigthLower); elem += width, j++, dy *= ky)
		{
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 3;
			y += dy;
		}



		// 4
		i = widthLeft;
		dx = hx * pow(kx, i);
		x = xLeft + hx * (1 - pow(kx, i)) / (1 - kx);

		j = heigthLower - 1;
		dy = hy * pow(ky, j);
		y = yLower + hy * (1 - pow(ky, j)) / (1 - ky);

		for (elem = width * heigthLower - widthRight; elem < width * heigthLower - 1; elem++, i++, dx *= kx)
		{
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 4;
			x += dx;
		}



		// 5
		i = widthLeft - 1;
		x = xLeft + hx * (1 - pow(kx, i)) / (1 - kx);

		j = heigthLower;
		dy = hy * pow(ky, j);
		y = yLower + hy * (1 - pow(ky, j)) / (1 - ky);

		for (elem = width * j + widthLeft - 1; elem < elemCount; elem += width, j++, dy *= ky)
		{
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 5;
			y += dy;
		}



		// 6
		i = 0;
		dx = hx;
		x = xLeft;

		j = heigth - 1;
		y = yLower + hy * (1 - pow(ky, j)) / (1 - ky);

		for (elem = width * j; elem < width * j + widthLeft - 1; elem++, i++, dx *= kx)
		{
			nodes[elem].setNodesData(x, y, i, j, condType, coef);
			nodes[elem].border = 6;
			x += dx;
		}



		// fictitious nodes
		j = heigthLower;
		dy = hy * pow(ky, j);
		y = yLower + hy * (1 - pow(ky, j)) / (1 - ky);

		for (j = heigthLower; j < heigth; j++, dy *= ky)
		{
			i = widthLeft;
			dx = hx * pow(kx, i);
			x = xLeft + hx * (1 - pow(kx, i)) / (1 - kx);
			for (elem = width * j + widthLeft; elem < width * (j + 1); elem++, i++, dx *= kx)
			{
				nodes[elem].setNodesData(x, y, i, j, -1, coef);
				x += dx;
			}
			y += dy;
		}
	}
}



// Отображние сетки на экран
void GRID::showGrid() {

	cout << "X:" << endl;
	for (size_t j = 0; j < heigth; j++)
	{
		size_t elem = j * width;
		for (size_t i = 0; i < width; i++, elem++)
		{
			cout << nodes[elem].x << " ";
		}
		cout << endl;
	}

	cout << endl << endl << "Y:" << endl;

	for (size_t j = 0; j < heigth; j++)
	{
		size_t elem = j * width;
		for (size_t i = 0; i < width; i++, elem++)
		{
			cout << nodes[elem].y << " ";
		}
		cout << endl;
	}
}



// Сохранение внутренних и внешних узлов в 2 файлах
void GRID::saveGridAndBorder(const string &filepathGrid, const string &filepathGridBorder) {

	ofstream grid(filepathGrid);
	ofstream border(filepathGridBorder);
	for (size_t i = 0; i < elemCount; i++)
	{
		if (nodes[i].type > 0)
			border << nodes[i].x << " " << nodes[i].y << endl;
		else
			grid << nodes[i].x << " " << nodes[i].y << endl;
	}

	border.close();
	grid.close();
}