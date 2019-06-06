#include "solver.h"



// A*x = b, где x - произвольный вектор
vector1D SOLVER::multA(const vector1D& x) {

	vector1D result(n);
	for (int i = 0; i < n; ++i) {

		result[i] = di[i] * x[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];

		for (int k = i0; k < i1; ++k) {

			result[i] += al[k] * x[ja[k]];
			result[ja[k]] += au[k] * x[i];
		}
	}

	return result;
}

// A*x = b, где x - произвольный вектор
vector1D SOLVER::multD(const vector1D&x) {

	vector1D result(n);

	for (int i = 0; i < n; ++i)
		result[i] = di_f[i] * x[i];

	return result;
}


// Создаём вектор x* = (1,2,...n)'
void SOLVER::generateVectX(int size) {

	x.resize(size);
	for (int i = 0; i < size; ++i) {
		x[i] = i + 1;
	}
}


// Вывод вектора b в файл 
void SOLVER::writeFToFile(const char *fileName) {

	std::ofstream fout;
	fout.open(fileName);

	for (int i = 0; i < b.size(); ++i)
		fout << b[i] << endl;
	fout.close();
}


// Вывод вектора x в файл
void SOLVER::writeXToFile(const char * fileName) {

	std::ofstream fout;
	fout.open(fileName);
	for (int i = 0; i < x.size(); ++i)
		fout << x[i] << " ";
	fout << " \t";
	fout.close();
}


// Вывод вектора x в поток
void SOLVER::writeXToStream(std::ofstream& fout) {

	for (int i = 0; i < x.size(); ++i)
		fout << x[i] << "\n";
	fout << "\n";
}


// Диагональное предобуславливание M = D
void SOLVER::decomposionD() {

	di_f.clear();
	di_f.resize(n);
	for (int i = 0; i < n; ++i)
		di_f[i] = 1.0 / sqrt(di[i]);
}


// LU_sq разложение матрицы А
void SOLVER::decomposionLUsq() {

	double sum_u, sum_l, sum_d;
	di_f = di;
	al_f = al;
	au_f = au;

	// Идём построчно в верхнем треугольнике, что экививалентно
	// Обходу нижнего треугольника по столбцам вниз, начиная с первого
	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];

		// Рассчёт элементов нижнего треугольника
		for (int k = i0; k < i1; ++k) {

			int j = ja[k]; // текущий j
			int j0 = ia[j]; // i0 строки j
			int j1 = ia[j + 1]; // i1 строки j
			sum_l = 0;
			sum_u = 0;
			int ki = i0; // Индекс l_ik
			int kj = j0; // Индекс u_kj

			while (ki < k && kj < j1) {

				if (ja[ki] == ja[kj]) { // l_ik * u_kj
					sum_l += al_f[ki] * au_f[kj];
					sum_u += au_f[ki] * al_f[kj];
					ki++;
					kj++;
				}
				else { // Ищем следующие элементы i и j строки, которые можем перемножить
					if (ja[ki] > ja[kj]) kj++;
					else ki++;
				}
			}

			al_f[k] = (al_f[k] - sum_l) / di_f[j];
			au_f[k] = (au_f[k] - sum_u) / di_f[j];
		}


		// Рассчёт диагонального элемента
		sum_d = 0.0;
		for (int k = i0; k < i1; ++k)
			sum_d += al_f[k] * au_f[k];
		di_f[i] = sqrt(di_f[i] - sum_d);
	}
}


// Прямой ход    L y = b    ==>    y = L^-1 b
vector1D SOLVER::execDirectTraversal(const vector1D &_F) {

	vector1D y;
	y.resize(n, 0);

	for (int i = 0; i < n; ++i) {
		double sum = 0;
		int i0 = ia[i];
		int i1 = ia[i + 1];
		for (int k = i0; k < i1; ++k)
			sum += al_f[k] * y[ja[k]];

		y[i] = (_F[i] - sum) / di_f[i];
	}
	return y;
}


// Обратный ход    U(sq) x = y    ==>    x = U(sq)^-1 y
vector1D SOLVER::execReverseTraversal(const vector1D &_y) {

	vector1D x, y = _y;
	x.resize(n);
	for (int i = n - 1; i >= 0; --i) {

		x[i] = y[i] / di_f[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];
		for (int k = i0; k < i1; ++k)
			y[ja[k]] -= au_f[k] * x[i];
	}

	return x;
}



// Полная очистка СЛАУ
void SOLVER::clearAll() {

	n = 0;
	E = 0.0;
	maxiter = 0;

	di.clear();
	ia.clear();
	ja.clear();
	al.clear();
	au.clear();

	di_f.clear();
	al_f.clear();
	au_f.clear();

	x.clear();
	r.clear();
	z.clear();
	p.clear();
	b.clear();
	bTmp.clear();
}



// Рассчёт относительной невязки
double SOLVER::calcRelativeDiscrepancy() {
	//return calcNormE(r) / calcNormE(b);
	return (r*r);
}


// Локально - оптимальная схема
int SOLVER::LOS() {

	x.clear();			// Задаём начальное приближение
	x.resize(n, 0);		// x_0 = (0, 0, ...)
	r.resize(n);

	vector1D xprev = x;
	r = b - multA(x);	// r_0 = b - A*x_0
	z = r;				// z_0 = r_0
	p = multA(z);		// p_0 = A*z_0


	for (int i = 0; i < maxiter; ++i) {

		double pp = (p * p);
		double alpha = (p * r) / pp;

		x = x + alpha * z;
		r = r - alpha * p;

		bTmp = multA(r);
		double beta = -(p * bTmp) / pp;

		z = r + beta * z;
		p = bTmp + beta * p;


		double relativeDiscrepancy = calcRelativeDiscrepancy();
		if (x == xprev || relativeDiscrepancy < E) {
			return i;
		}
		xprev = x;
	}
}


// Локально - оптимальная схема c неполной диагональной факторизацией
int SOLVER::LOSfactD() {

	x.clear();			// Задаём начальное приближение
	x.resize(n, 0);		// x_0 = (0, 0, ...)
	vector1D xprev = x;
	decomposionD();

	r = b - multA(x);	// r_0 = b - A*x_0
	r = multD(r);
	z = multD(r);		// z = U^-1 r
	p = multA(z);		// p = A*z
	p = multD(p);		// p = L^-1 A*z


	for (int i = 0; i < maxiter; ++i) {

		double pp = p * p;
		double alpha = (p*r) / pp;
		x = x + alpha * z;
		r = r - alpha * p;

		vector1D tmp = multD(r);
		tmp = multA(tmp);
		tmp = multD(tmp);
		double beta = -(p * tmp) / pp;
		p = tmp + beta * p;

		tmp = multD(r);
		z = tmp + beta * z;


		double relativeDiscrepancy = calcRelativeDiscrepancy();
		if (x == xprev || relativeDiscrepancy < E) {
			return i;
		}
		xprev = x;
	}
}





// Локально - оптимальная схема с неполной факторизацией LU(sq)
int SOLVER::LOSfactLUsq() {

	x.clear();						// Задаём начальное приближение
	x.resize(n, 0);					// x_0 = (0, 0, ...)
	vector1D xprev = x;
	decomposionLUsq();

	r = b - multA(x);				// r_0 = b - A*x_0
	r = execDirectTraversal(r);		// r = L^-1 (b - A*x_0)
	z = execReverseTraversal(r);	// z = U^-1 r
	p = multA(z);					// p = A*z
	p = execDirectTraversal(p);		// p = L^-1 A*z


	for (int i = 0; i < maxiter; ++i) {

		double pp = p * p;
		double alpha = (p*r) / pp;
		x = x + alpha * z;
		r = r - alpha * p;

		vector1D tmp = execReverseTraversal(r);
		tmp = multA(tmp);
		tmp = execDirectTraversal(tmp);
		double beta = -(p * tmp) / pp;
		p = tmp + beta * p;

		tmp = execReverseTraversal(r);
		z = tmp + beta * z;


		double relativeDiscrepancy = calcRelativeDiscrepancy();
		if (x == xprev || relativeDiscrepancy < E) {
			return i;
		}
		xprev = x;
	}
}




// Считываем все СЛАУ и её параметры из файлов
void SOLVER::initSLAE()
{
	ifstream finMatrix("input/matrix.txt");
	ifstream finVector("input/vector.txt");
	ifstream finParams("input/SLAE_parameters.txt");

	finParams >> maxiter >> E >> delta;
	finMatrix >> n;
	di.resize(n);
	ia.resize(n + 1);

	for (size_t i = 0; i < n; i++)
		finMatrix >> di[i];
	for (size_t i = 0; i < n + 1; i++)
		finMatrix >> ia[i];

	ja.resize(ia[n]);
	al.resize(ia[n]);
	au.resize(ia[n]);

	for (size_t i = 0; i < ja.size(); i++)
		finMatrix >> ja[i];
	for (size_t i = 0; i < al.size(); i++)
		finMatrix >> al[i];
	for (size_t i = 0; i < au.size(); i++)
		finMatrix >> au[i];

	b.resize(n);
	for (size_t i = 0; i < n; i++)
		finVector >> b[i];

	finMatrix.close();
	finVector.close();
	finParams.close();
}



// Метод бисопряжённых градиентов
int SOLVER::BiCG()
{
	ofstream fout("output/result.txt");
	int i = 0;
	double prPrev, pr;
	vector1D Az, Ats;
	generateInitialGuess();
	r = b - multAOn(x);
	p = z = s = r;
	do {
		xPrev = x;
		prPrev = (p*r);
		Az = multAOn(z);
		Ats = multAtOn(s);

		alpha = prPrev / (s*Az);

		x = x + alpha * z;
		r = r - alpha * Az;
		p = p - alpha * Ats;

		beta = (p*r) / prPrev;

		z = r + beta * z;
		s = p + beta * s;

		i++;
	} while (!doStop(i));

	fout.close();
	return i;
}



// Создание начального приближения x0 = (0,...,0)'
void SOLVER::generateInitialGuess()
{
	x.resize(n, 1);
}



// Умножение матрицы А на вектор
vector1D SOLVER::multAOn(const vector1D &v)
{
	vector1D result(n);
	for (size_t i = 0; i < n; i++)
	{
		result[i] = di[i] * v[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];

		for (size_t k = i0; k < i1; k++)
		{
			int j = ja[k];
			result[i] += al[k] * v[j];
			result[j] += au[k] * v[i];
		}
	}
	return result;
}



// Умножение транспонированной матрицы А на вектор
vector1D SOLVER::multAtOn(const vector1D &v)
{
	vector1D result(n);
	for (size_t i = 0; i < n; i++)
	{
		result[i] = di[i] * v[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];
		for (size_t k = i0; k < i1; k++)
		{
			int j = ja[k];
			result[j] += al[k] * v[i];
			result[i] += au[k] * v[j];
		}
	}
	return result;
}



// Проверяем условие выхода
bool SOLVER::doStop(int i)
{
	// Выход по числу итераций
	if (i > maxiter)
		return true;

	// Выход шагу
	if (calcNormE(x - xPrev) / calcNormE(x) < delta)
		return true;

	// Выход по относительной невязке
	if (calcNormE(multAOn(x) - b) / calcNormE(b) < E)
		return true;

	return false;
}
