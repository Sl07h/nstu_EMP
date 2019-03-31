#include "grid.h"



void GRID::inputGrid()
{
	string filename;
	if (isGridUniform)
		filename = "input/uniform_grid.txt";
	else
		filename = "input/nonuniform_grid.txt";

	std::ifstream fin(filename);

	fin >> xLeft >> xRight >>
		yLower >> yUpper;

	fin >> width >> heigth >>
		widthLeft >> widthRight >>
		heigthLower >> heigthUpper;

	if (!isGridUniform) {
		fin >> stepX0 >> stepY0;
		fin >> dx >> dy;
	}

	fin >> condType;


	elemCount = width * heigth;


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
	nodes.resize(elemCount);

	if (isGridUniform) {
		hx = (xRight - xLeft) / double(width);
		hy = (yUpper - yLower) / double(heigth);


		size_t i, j, elem;
		double x, y;
		// Обходим все внутренние элементы нижней части "L" по пяти точкам
		for (j = 1; j < heigthLower - 1; j++)
		{
			i = 1;
			for (elem = j * width + 1; elem < (j + 1) * width - 1; elem++, i++)
			{
				x = xLeft + (xRight - xLeft) * i / width;
				y = yLower + (yUpper - yLower) * j / heigth;
				
				nodes[elem].setNodesData(x, y, i, j, 0);
			}
		}

		// Обходим все внутренние элементы средней части части "L" по пяти точкам
		i = 1;
		j = heigthLower - 1;
		y = yLower + (yUpper - yLower) * j / heigth;
		for (elem = j * width + 1; elem < j * width + widthLeft; elem++, i++)
		{
			x = xLeft + (xRight - xLeft) * i / width;

			nodes[elem].setNodesData(x, y, i, j, 0);
		}


		// Обходим все внутренние элементы верхней части "L" по пяти точкам
		for (j = heigthLower; j < heigth - 1; j++)
		{
			i = 1;
			for (elem = j * width + 1; elem < j * width + widthLeft - 1; elem++, i++)
			{
				x = xLeft + (xRight - xLeft) * i / width;
				y = yLower + (yUpper - yLower) * j / heigth;
				
				nodes[elem].setNodesData(x, y, i, j, 0);
			}
		}

		// 1
		i = 0;
		j = 0;
		x = xLeft + (xRight - xLeft) * i / width;
		for (elem = 0; elem < (heigth - 1) * width; elem += width, j++)
		{
			y = yLower + (yUpper - yLower) * j / heigth;
			
			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 1;
		}

		// 2
		i = 1;
		j = 0;
		y = yLower + (yUpper - yLower) * j / heigth;
		for (elem = 1; elem < width; elem++, i++)
		{
			x = xLeft + (xRight - xLeft) * i / width;
			
			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 2;
		}

		// 3
		i = width - 1;
		j = 1;
		x = xLeft + (xRight - xLeft) * i / width;
		for (elem = 2 * width - 1; elem < (width * heigthLower); elem += width, j++)
		{

			y = yLower + (yUpper - yLower) * j / heigth;
			
			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 3;
		}

		// 4
		i = widthLeft;
		j = heigthLower - 1;
		y = yLower + (yUpper - yLower) * j / heigth;
		for (elem = width * heigthLower - widthRight; elem < width * heigthLower - 1; elem++, i++)
		{
			x = xLeft + (xRight - xLeft) * i / width;
			
			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 4;
		}

		// 5
		i = widthLeft - 1;
		j = heigthLower;
		x = xLeft + (xRight - xLeft) * i / width;
		for (elem = width * j + widthLeft - 1; elem < elemCount; elem += width, j++)
		{
			y = yLower + (yUpper - yLower) * j / heigth;
			
			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 5;
		}

		// 6
		i = 0;
		j = heigth - 1;
		y = yLower + (yUpper - yLower) * j / heigth;
		for (elem = width * j; elem < width * j + widthLeft - 1; elem++, i++)
		{
			x = xLeft + (xRight - xLeft) * i / width;
			
			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 6;
		}

		// fictitious nodes
		for (j = heigthLower; j < heigth; j++)
		{
			i = widthLeft;
			for (elem = width * j + widthLeft; elem < width * (j + 1); elem++, i++)
			{
				x = xLeft + (xRight - xLeft) * i / width;
				y = yLower + (yUpper - yLower) * j / heigth;
				
				nodes[elem].setNodesData(x, y, i, j, -1);
			}
		}

	}
	else {

		double x, y, stepX = 0, stepY = 0;
		size_t i, j, elem;


		// Обходим все внутренние элементы нижней части "L" по пяти точкам
		stepY = stepY0;
		for (j = 1; j < heigthLower - 1; j++)
		{
			stepY = stepY0 + (j - 1) * dy;
			y = yLower + (stepY0 + stepY) * j / 2;

			i = 1;
			for (elem = j * width + 1; elem < (j + 1) * width - 1; elem++, i++)
			{
				stepX = stepX0 + (i - 1) * dx;
				x = xLeft + (stepX0 + stepX) * i / 2;

				nodes[elem].setNodesData(x, y, i, j, 0);
			}
		}



		// Обходим все внутренние элементы средней части части "L" по пяти точкам
		i = 1;
		j = heigthLower - 1;
		stepY = stepY0 + (j - 1) * dy;
		y = yLower + (stepY0 + stepY) * j / 2;
		for (elem = j * width + 1; elem < j * width + widthLeft; elem++, i++)
		{
			stepX = stepX0 + (i - 1) * dx;
			x = xLeft + (stepX0 + stepX) * i / 2;

			nodes[elem].setNodesData(x, y, i, j, 0);
		}


		// Обходим все внутренние элементы верхней части "L" по пяти точкам
		for (j = heigthLower; j < heigth - 1; j++)
		{
			stepY = stepY0 + (j - 1) * dy;
			y = yLower + (stepY0 + stepY) * j / 2;

			i = 1;
			for (elem = j * width + 1; elem < j * width + widthLeft - 1; elem++, i++)
			{
				stepX = stepX0 + (i - 1) * dx;
				x = xLeft + (stepX0 + stepX) * i / 2;

				nodes[elem].setNodesData(x, y, i, j, 0);
			}
		}

		// 1
		i = 0;
		j = 0;
		stepX = stepX0 + (i - 1) * dx;
		x = xLeft + (stepX0 + stepX) * i / 2;
		for (elem = 0; elem < (heigth - 1) * width; elem += width, j++)
		{
			stepY = stepY0 + (j - 1) * dy;
			y = yLower + (stepY0 + stepY) * j / 2;

			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 1;
		}

		// 2
		i = 1;
		j = 0;
		stepY = stepY0 + (j - 1) * dy;
		y = yLower + (stepY0 + stepY) * j / 2;
		for (elem = 1; elem < width; elem++, i++)
		{
			stepX = stepX0 + (i - 1) * dx;
			x = xLeft + (stepX0 + stepX) * i / 2;

			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 2;
		}

		// 3
		i = width - 1;
		j = 1;
		stepX = stepX0 + (i - 1) * dx;
		x = xLeft + (stepX0 + stepX) * i / 2;
		for (elem = 2 * width - 1; elem < (width * heigthLower); elem += width, j++)
		{
			stepY = stepY0 + (j - 1) * dy;
			y = yLower + (stepY0 + stepY) * j / 2;

			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 3;
		}

		// 4
		i = widthLeft;
		j = heigthLower - 1;
		stepY = stepY0 + (j - 1) * dy;
		y = yLower + (stepY0 + stepY) * j / 2;
		for (elem = width * heigthLower - widthRight; elem < width * heigthLower - 1; elem++, i++)
		{
			stepX = stepX0 + (i - 1) * dx;
			x = xLeft + (stepX0 + stepX) * i / 2;

			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 4;
		}

		// 5
		i = widthLeft - 1;
		j = heigthLower;
		stepX = stepX0 + (i - 1) * dx;
		x = xLeft + (stepX0 + stepX) * i / 2;
		for (elem = width * j + widthLeft - 1; elem < elemCount; elem += width, j++)
		{
			stepY = stepY0 + (j - 1) * dy;
			y = yLower + (stepY0 + stepY) * j / 2;

			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 5;
		}

		// 6
		i = 0;
		j = heigth - 1;
		stepY = stepY0 + (j - 1) * dy;
		y = yLower + (stepY0 + stepY) * j / 2;
		for (elem = width * j; elem < width * j + widthLeft - 1; elem++, i++)
		{
			stepX = stepX0 + (i - 1) * dx;
			x = xLeft + (stepX0 + stepX) * i / 2;

			nodes[elem].setNodesData(x, y, i, j, condType);
			nodes[elem].border = 6;
		}

		// fictitious nodes
		for (j = heigthLower; j < heigth; j++)
		{
			stepY = stepY0 + (j - 1) * dy;
			y = yLower + (stepY0 + stepY) * j / 2;

			i = widthLeft;
			for (elem = width * j + widthLeft; elem < width * (j + 1); elem++, i++)
			{
				stepX = stepX0 + (i - 1) * dx;
				x = xLeft + (stepX0 + stepX) * i / 2;
				
				nodes[elem].setNodesData(x, y, i, j, -1);
			}
		}
	}
}
