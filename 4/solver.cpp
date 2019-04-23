#include "solver.h"



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
void SOLVER::BiCG()
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

	fout << x << endl
		<< "Iterations: " << i;

	fout.close();
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
