#include "head.h"


class matrix {

public:

	int readMatrixFromFile(const char* filename);
	void writeMatrixToFile(const char* filename);
	int getDimention() { return n; }
	void setE(real new_E) { E = new_E; }
	void setMaxiter(real new_maxiter) { maxiter = new_maxiter; }
	void generateMatrixWith7Diagonals(int new_n, int new_m);
	void invertSigns();

protected:

	real calcAii(int i);

	vector <real> di, au1, au2, au3, al1, al2, al3;
	int n, m, maxiter;
	real E;
};