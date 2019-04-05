#include "head.h"
#include "fdm.h"





// Инициализируем модель, задавая функции u, f и тип сетки
void FDM::init(function2D & _u, function2D & _f, bool _isGridUniform, int _condType, int _coef)
{
	u = _u;
	f = _f;
	isGridUniform = _isGridUniform;
	condType = _condType;
	coef = _coef;
}


// Ввод параметров точности и коэффициентов уравнения из файла
void FDM::inputEquationParameters() {

	std::ifstream fin("equation_parameters.txt");

	fin >> lambda >> gamma;

	fin.close();
}


// Вывод слау в файлы
void FDM::outputSLAE(const string & fileA)
{
	ofstream foutA(fileA);

	//foutA << "A = [ ";
	for (size_t i = 0; i < elemCount; i++)
	{
		for (size_t j = 0; j < elemCount; j++)
			foutA << A[i][j] << "\t";

		foutA << ";" << endl;
	}

	//foutA << "b = [ ";
	for (size_t i = 0; i < elemCount; i++)
	{
		foutA << b[i] << ";" << endl;
	}
}


// Преобразуем сетку в СЛАУ
void FDM::transformGridToSLAE()
{
	n = elemCount;
	m = width;
	di.resize(elemCount);
	al1.resize(elemCount - 1);
	al2.resize(elemCount - width);
	au1.resize(elemCount - 1);
	au2.resize(elemCount - width);
	A2.resize(elemCount);
	for (size_t i = 0; i < elemCount; i++)
	{
		A2[i].resize(elemCount);
	}

	b.resize(elemCount);
	x.resize(elemCount);



	if (isGridUniform) {

		cout << "Grid is uniform" << endl;
		cout << hx << "\t" << hy << endl;

		for (size_t elem = 0; elem < elemCount; elem++)
		{
			int i = nodes[elem].i;
			int j = nodes[elem].j;
			double x = nodes[elem].x;
			double y = nodes[elem].y;

			switch (nodes[elem].type)
			{

			case -1: {	// фиктивные узлы
				di[elem] = 1;
				A2[elem][elem] = 1;
				b[elem] = 0;
			} break;


			case 0: {	// внутренние узлы
				di[elem] = -2 / (hx * hx) - 2 / (hy * hy);
				al1[elem - 1] = 1 / (hx * hx);
				au1[elem] = 1 / (hx * hx);
				al2[elem - width] = 1 / (hy * hy);
				au2[elem] = 1 / (hy * hy);

				A2[elem][elem] = -2 / (hx * hx) - 2 / (hy * hy);
				A2[elem][elem - 1] = 1 / (hx * hx);
				A2[elem][elem + 1] = 1 / (hx * hx);
				A2[elem][elem - width] = 1 / (hy * hy);
				A2[elem][elem + width] = 1 / (hy * hy);
				b[elem] = f(x, y);
			} break;


			case 1: { // 1-е краевые условия
				di[elem] = 1;
				A2[elem][elem] = 1;
				b[elem] = u(x, y);
			} break;


			case 3: {
				switch (nodes[elem].border)
				{
				case 1: {
					au1[elem] = 1.0 / hx;
					di[elem] = -1.0 / hx + C1;
					 
					A2[elem][elem + 1] = 1.0 / hx;
					A2[elem][elem] = -1.0 / hx + C1;
					b[elem] = calcFirstDerivativeX(x, y) + C1 * u(x, y) + C2;
					//
					// <---[elem][elem+1]
					//
				}break;

				default: { // 1-е краевые условия
					di[elem] = 1;
					A2[elem][elem] = 1;
					b[elem] = u(x, y);
				} break;
				}

			}break;
			}
		}

	}
	else {

		cout << "Grid is nonuniform" << endl;
		for (size_t elem = 0; elem < elemCount; elem++)
		{
			int i = nodes[elem].i;
			int j = nodes[elem].j;
			double x = nodes[elem].x;
			double y = nodes[elem].y;


			switch (nodes[elem].type)
			{
			case -1: {	// фиктивные узлы
				A2[elem][elem] = 1;
				di[elem] = 1;
				b[elem] = 0;
			} break;

			case 0: {	// внутренние узлы

				/* // Я памятник се6е воздвиг нерукотворный
				// Металлов твёрже он и выше пирамид...
				hx = xLeft + i * dx;
				hy = yLower + j * dy;
				double hxPrev = hx - dx;
				double hyPrev = hy - dy;
				*/

				hx = nodes[elem + 1].x - nodes[elem].x;
				hy = nodes[elem + width].y - nodes[elem].y;
				hxPrev = nodes[elem].x - nodes[elem - 1].x;
				hyPrev = nodes[elem].y - nodes[elem - width].y;

				di[elem] = -2.0 / (hxPrev * hx) - 2.0 / (hyPrev * hy);
				al1[elem - 1] = 2.0 / (hxPrev*(hx + hxPrev));
				au1[elem] = 2.0 / (hx*(hx + hxPrev));
				al2[elem - width] = 2.0 / (hyPrev*(hy + hyPrev));
				au2[elem] = 2.0 / (hy*(hy + hyPrev));

				A2[elem][elem] = -2.0 / (hxPrev * hx) - 2.0 / (hyPrev * hy);
				A2[elem][elem - 1] = 2.0 / (hxPrev*(hx + hxPrev));
				A2[elem][elem + 1] = 2.0 / (hx*(hx + hxPrev));
				A2[elem][elem - width] = 2.0 / (hyPrev*(hy + hyPrev));
				A2[elem][elem + width] = 2.0 / (hy*(hy + hyPrev));
				b[elem] = f(x, y);
			} break;

			case 1: { // 1-е краевые условия
				A2[elem][elem] = 1;
				di[elem] = 1;
				b[elem] = u(x, y);
			} break;
			case 3: {
				switch (nodes[elem].border)
				{
				case 1: {
					hx = nodes[elem + 1].x - nodes[elem].x;

					au1[elem] = 1.0 / hx;
					di[elem] = -1.0 / hx + C1;

					A2[elem][elem + 1] = 1.0 / hx;
					A2[elem][elem] = -1.0 / hx + C1;
					b[elem] = calcFirstDerivativeX(x, y) + C1 * u(x, y) + C2;
					//
					// <---[elem][elem+1]
					//
				}break;

				default: { // 1-е краевые условия
					di[elem] = 1;
					A2[elem][elem] = 1;
					b[elem] = u(x, y);
				} break;
				}

			}break;
			}
		}

	}
}


//
//// Преобразуем сетку в СЛАУ
//void FDM::transformGridToSLAE()
//{
//	n = elemCount;
//	m = width;
//	di.resize(elemCount);
//	al1.resize(elemCount - 1);
//	al2.resize(elemCount - width);
//	au1.resize(elemCount - 1);
//	au2.resize(elemCount - width);
//	A2.resize(elemCount);
//	for (size_t i = 0; i < elemCount; i++)
//	{
//		A2[i].resize(elemCount);
//	}
//
//	b.resize(elemCount);
//	x.resize(elemCount);
//
//
//
//	if (isGridUniform) {
//
//		cout << "Grid is uniform" << endl;
//		cout << hx << "\t" << hy << endl;
//
//		for (size_t elem = 0; elem < elemCount; elem++)
//		{
//			int i = nodes[elem].i;
//			int j = nodes[elem].j;
//			double x = nodes[elem].x;
//			double y = nodes[elem].y;
//
//			switch (nodes[elem].type)
//			{
//
//			case -1: {	// фиктивные узлы
//				di[elem] = 1;
//				A2[elem][elem] = 1;
//				b[elem] = 0;
//			} break;
//
//
//			case 0: {	// внутренние узлы
//				di[elem] = -2 / (hx * hx) - 2 / (hy * hy);
//				al1[elem - 1] = 1 / (hx * hx);
//				au1[elem] = 1 / (hx * hx);
//				al2[elem - width] = 1 / (hy * hy);
//				au2[elem] = 1 / (hy * hy);
//
//				A2[elem][elem] = -2 / (hx * hx) - 2 / (hy * hy);
//				A2[elem][elem - 1] = 1 / (hx * hx);
//				A2[elem][elem + 1] = 1 / (hx * hx);
//				A2[elem][elem - width] = 1 / (hy * hy);
//				A2[elem][elem + width] = 1 / (hy * hy);
//				b[elem] = f(x, y);
//			} break;
//
//
//			case 1: { // 1-е краевые условия
//				di[elem] = 1;
//				A2[elem][elem] = 1;
//				b[elem] = u(x, y);
//			} break;
//
//
//			case 3: {
//				switch (nodes[elem].border)
//				{
//				case 1: {
//					au1[elem] = 1.0 / hx;
//					di[elem] = -1.0 / hx + C1;
//
//					A2[elem][elem + 1] = 1.0 / hx;
//					A2[elem][elem] = -1.0 / hx + C1;
//					b[elem] = calcFirstDerivativeX(x, y) + C1 * u(x, y) + C2;
//					//
//					// <---[elem][elem+1]
//					//
//				}break;
//
//				case 2: {
//					au2[elem] = 1.0 / hy;
//					di[elem] = -1.0 / hy + C1;
//
//					A2[elem][elem + width] = 1.0 / hy;
//					A2[elem][elem] = -1.0 / hy + C1;
//					b[elem] = calcFirstDerivativeY(x, y) + C1 * u(x, y) + C2;
//					// [elem+width]
//					// [elem]
//					//    |
//				}break;
//
//				case 3: {
//					al1[elem - 1] = 1.0 / hx;
//					di[elem] = -1.0 / hx + C1;
//
//					A2[elem][elem - 1] = 1.0 / hx;
//					A2[elem][elem] = -1.0 / hx + C1;
//					b[elem] = calcFirstDerivativeX(x, y) + C1 * u(x, y) + C2;
//					//
//					// [elem-1][elem]--->
//					//
//				}break;
//
//				case 4: {
//					al2[elem - width] = 1.0 / hy;
//					di[elem] = -1.0 / hy + C1;
//
//					A2[elem][elem - width] = 1.0 / hy;
//					A2[elem][elem] = -1.0 / hy + C1;
//					b[elem] = calcFirstDerivativeY(x, y) + C1 * u(x, y) + C2;
//					//   |
//					// [elem]
//					// [elem-width]
//				}break;
//
//				case 5: {
//					al1[elem - 1] = 1.0 / hx;
//					di[elem] = -1.0 / hx + C1;
//
//					A2[elem][elem - 1] = 1.0 / hx;
//					A2[elem][elem] = -1.0 / hx + C1;
//					b[elem] = calcFirstDerivativeX(x, y) + C1 * u(x, y) + C2;
//					//
//					// [elem-1][elem]--->
//					//
//				}break;
//
//				case 6: {
//					al2[elem - width] = 1.0 / hy;
//					di[elem] = -1.0 / hy + C1;
//
//					A2[elem][elem - width] = 1.0 / hy;
//					A2[elem][elem] = -1.0 / hy + C1;
//					b[elem] = calcFirstDerivativeY(x, y) + C1 * u(x, y) + C2;
//					//   |
//					// [elem]
//					// [elem-width]
//				}break;
//				}
//
//			}break;
//			}
//		}
//
//	}
//	else {
//
//		cout << "Grid is nonuniform" << endl;
//		for (size_t elem = 0; elem < elemCount; elem++)
//		{
//			int i = nodes[elem].i;
//			int j = nodes[elem].j;
//			double x = nodes[elem].x;
//			double y = nodes[elem].y;
//
//
//			switch (nodes[elem].type)
//			{
//			case -1: {	// фиктивные узлы
//				A2[elem][elem] = 1;
//				di[elem] = 1;
//				b[elem] = 0;
//			} break;
//
//			case 0: {	// внутренние узлы
//
//				/* // Я памятник се6е воздвиг нерукотворный
//				// Металлов твёрже он и выше пирамид...
//				hx = xLeft + i * dx;
//				hy = yLower + j * dy;
//				double hxPrev = hx - dx;
//				double hyPrev = hy - dy;
//				*/
//
//				hx = nodes[elem + 1].x - nodes[elem].x;
//				hy = nodes[elem + width].y - nodes[elem].y;
//				hxPrev = nodes[elem].x - nodes[elem - 1].x;
//				hyPrev = nodes[elem].y - nodes[elem - width].y;
//
//				di[elem] = -2.0 / (hxPrev * hx) - 2.0 / (hyPrev * hy);
//				al1[elem - 1] = 2.0 / (hxPrev*(hx + hxPrev));
//				au1[elem] = 2.0 / (hx*(hx + hxPrev));
//				al2[elem - width] = 2.0 / (hyPrev*(hy + hyPrev));
//				au2[elem] = 2.0 / (hy*(hy + hyPrev));
//
//				A2[elem][elem] = -2.0 / (hxPrev * hx) - 2.0 / (hyPrev * hy);
//				A2[elem][elem - 1] = 2.0 / (hxPrev*(hx + hxPrev));
//				A2[elem][elem + 1] = 2.0 / (hx*(hx + hxPrev));
//				A2[elem][elem - width] = 2.0 / (hyPrev*(hy + hyPrev));
//				A2[elem][elem + width] = 2.0 / (hy*(hy + hyPrev));
//				b[elem] = f(x, y);
//			} break;
//
//			case 1: { // 1-е краевые условия
//				A2[elem][elem] = 1;
//				di[elem] = 1;
//				b[elem] = u(x, y);
//			} break;
//
//			case 3: {
//				switch (nodes[elem].border)
//				{
//				case 1: {
//					hx = nodes[elem + 1].x - nodes[elem].x;
//					au1[elem] = 1.0 / hx;
//					di[elem] = -1.0 / hx + C1;
//
//					A2[elem][elem + 1] = 1.0 / hx;
//					A2[elem][elem] = -1.0 / hx + C1;
//					b[elem] = calcFirstDerivativeX(x, y) + C1 * u(x, y) + C2;
//					//
//					// <---[elem][elem+1]
//					//
//				}break;
//
//				case 2: {
//					hy = nodes[elem + width].y - nodes[elem].y;
//					au2[elem] = 1.0 / hy;
//					di[elem] = -1.0 / hy + C1;
//
//					A2[elem][elem + width] = 1.0 / hy;
//					A2[elem][elem] = -1.0 / hy + C1;
//					b[elem] = calcFirstDerivativeY(x, y) + C1 * u(x, y) + C2;
//					// [elem+width]
//					// [elem]
//					//    |
//				}break;
//
//				case 3: {
//					hx = nodes[elem - 1].x - nodes[elem].x;
//					al1[elem - 1] = 1.0 / hx;
//					di[elem] = -1.0 / hx + C1;
//
//					A2[elem][elem - 1] = 1.0 / hx;
//					A2[elem][elem] = -1.0 / hx + C1;
//					b[elem] = calcFirstDerivativeX(x, y) + C1 * u(x, y) + C2;
//					//
//					// [elem-1][elem]--->
//					//
//				}break;
//
//				case 4: {
//					hy = nodes[elem - width].y - nodes[elem].y;
//					al2[elem - width] = 1.0 / hy;
//					di[elem] = -1.0 / hy + C1;
//
//					A2[elem][elem - width] = 1.0 / hy;
//					A2[elem][elem] = -1.0 / hy + C1;
//					b[elem] = calcFirstDerivativeY(x, y) + C1 * u(x, y) + C2;
//					//   |
//					// [elem]
//					// [elem-width]
//				}break;
//
//				case 5: {
//					hx = nodes[elem - 1].x - nodes[elem].x;
//					al1[elem - 1] = 1.0 / hx;
//					di[elem] = -1.0 / hx + C1;
//
//					A2[elem][elem - 1] = 1.0 / hx;
//					A2[elem][elem] = -1.0 / hx + C1;
//					b[elem] = calcFirstDerivativeX(x, y) + C1 * u(x, y) + C2;
//					//
//					// [elem-1][elem]--->
//					//
//				}break;
//
//				case 6: {
//					hy = nodes[elem - width].y - nodes[elem].y;
//					al2[elem - width] = 1.0 / hy;
//					di[elem] = -1.0 / hy + C1;
//
//					A2[elem][elem - width] = 1.0 / hy;
//					A2[elem][elem] = -1.0 / hy + C1;
//					b[elem] = calcFirstDerivativeY(x, y) + C1 * u(x, y) + C2;
//					//   |
//					// [elem]
//					// [elem-width]
//				}break;
//				}
//
//			}break;
//			}
//		}
//
//	}
//}
//




// Абсолютная невязка
double FDM::calcAbsResidual(const string &filename) {

	ofstream fout(filename);
	double tmp = 0.0, normAbsRes = 0.0;

	for (size_t elem = 0; elem < elemCount; elem++)
	{
		if (nodes[elem].type != -1) {
			tmp = abs(u(nodes[elem].x, nodes[elem].y) - x[elem]);
			fout << tmp << endl;
			normAbsRes += tmp;
		}
		else
			fout << 0 << endl;
	}

	fout.close();
	return sqrt(normAbsRes);
}
