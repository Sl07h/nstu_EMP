#pragma once
#include "head.h"

class SNE
{
public:
	void calcByNewton();
	void calcBySimpleIteration();



protected:
	vector <vector <double>> A;
	
	vector <double> di, al, au;
	vector <double> b;

	vector <double> q, qPrev;

};

