#include "head.h"
#include"matrix.h"
#include "vect.h"


class SLAE : public matrix, public vect {

public:
	void convMatrixToSparse();
	void convMatrixToDense();
	int readDenseMatrixFromFile(const char *filename);
	void writeDenseMatrixToFile(const char *filename);

	real multLine(vector<real>& line, int i, int mode);
	void mult();

	void Jacobi(real w);
	void GaussSeildel(real w);
	int calcIterative(int mode, real w);
	real findOptimalW(int mode, std::ofstream& fout);


protected:
	real calcNormE(vector <real> &x);
	real calcRelativeDiscrepancy();
	int calcCondNumber();
	vector <vector <real>> A;

};