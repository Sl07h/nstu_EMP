#include "fem.h"


void main() {

	string prefix = "";
	int coefGrid = 0;
	bool isGridUniform = false;
	if (!isGridUniform)
		prefix = "Non";

	string gridFile = "grids/" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";
	string gridBorderFile = "grids/Border" + prefix + "Uniform_" + to_string(coefGrid) + ".txt";

	FEM fem;
	fem.inputGrid();
	fem.buildGrid();
	fem.saveGridAndBorder(gridFile, gridBorderFile);
	
	
	string runVisualisation = "python plot.py " + gridFile + " " + gridBorderFile;
	system(runVisualisation.c_str());


}
