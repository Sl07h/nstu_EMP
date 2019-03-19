#include "head.h"


// Векторы xk and F
class vect {

public:

	int readFFromFile(const char *filename);
	void generateInitualGuess(int size) { xk.clear(); xk.resize(size, real(0)); }

	void getVectX(vector <real> &x) { x = F; };

	void generateVectX(int size);
	
	void writeTableToFile(std::ofstream& fout, int presision, real w, int iterations, int condNumber);
	void writeFToFile(const char *filename);
	bool isXcorrect();

protected:
	vector <real> xk, F;
};