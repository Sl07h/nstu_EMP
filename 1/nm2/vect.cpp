#include "head.h"
#include "vect.h"


// Ввод вектора правой части из файла
int vect::readFFromFile(const char * filename)
{
	std::ifstream fin(filename);

	for (size_t i = 0; i < F.size(); i++)
		fin >> F[i];

	return 0;
}

// Создаём вектор xk* = (1,2,...n)'
void vect::generateVectX(int size) {

	xk.resize(size);
	for (int i = 0; i < size; ++i) {
		xk[i] = i + 1;
	}
}


// Вывод вектора F в файл 
void vect::writeFToFile(const char *filename) {

	std::ofstream fout;
	fout.open(filename);

	for (int i = 0; i < F.size();++i)
		fout << int(F[i]) << endl;
	fout.close();
}


// generates 1/3 of table in research
void vect::writeTableToFile(std::ofstream& fout, int presision, real w, int iterations, int condNumber) {

	fout << std::fixed << std::setprecision(presision) << w << "\t";
	fout << std::fixed << std::setprecision(std::numeric_limits<real>::digits10 + 1);
	for (int i = 0; i < xk.size();++i)
		fout << xk[i] << " ";
	fout << " \t";

	fout << std::scientific;
	for (int i = 0; i < xk.size();++i)
		fout << xk[i] - real(i + 1) << " ";
	fout << " \t";

	fout << iterations << "\t";
	fout << condNumber << endl;
}




// Проверяем насколько xk  близко к вектору xk* = (1,2,...,n)'
bool vect::isXcorrect() {

	for (int i = 0; i < xk.size();++i) {

		if (abs(xk[i] - (real)(i + 1)) > std::numeric_limits<real>::digits10 + 2)
			return false;
	}

	return true;
}
