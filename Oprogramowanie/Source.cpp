#include <iostream>
#include <cstdio>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>

using namespace std;
double nKsi[4], nEta[4];
static double w1for3[9] = { 5.0/9.0, 8.0/9.0,
5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
static double w2for3[9] = { 5.0 / 9.0,
5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 };
static double w1for4[16] = { (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0,
							(18.0 - sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0,
							(18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };
static double w2for4[16] = { (18.0 - sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0,(18.0 + sqrt(30.0)) / 36.0,
							(18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0,(18.0 - sqrt(30.0)) / 36.0,
							(18.0 - sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };

//coordinates of single node
struct Node {
	float x, y;
	bool bc; 
};

//the id number and it's value
struct Element {
	int ID[4];
};

//grid with numbers of nodes and elements
struct Grid {
	int nN = 1; //number of nodes
	int nE = 1; //number of elements

	vector <Node> nodes; //all of the nodes on grid
	vector <Element> elements; ////all of the elements on grid
};

//global data which come from the file
struct GlobalData {
	int t;		//SimulationTime
	int st;		//SimulationStepTime 
	int lambda;	//Conductivity
	int alpha;	//Alfa 
	int tot;	//Tot
	int it;		//InitialTime
	int d;		//Density
	int cp;		//SpecificHeat
};

//scheme of integration, with points of integration and weights
struct Integration {
	double nodes2[2] = { -1 / sqrt(3), 1 / sqrt(3) };
	double weights2[2] = { 1,1 };

	double nodes3[3] = { -sqrt(0.6), 0 , sqrt(0.6) };
	double weights3[3] = { 5.0/9.0, 8.0/9.0, 5.0/9.0 };
};

//structure for derivatives with 2points
struct Elem4 {
	double** tabKsi = new double*[4];
	double** tabEta = new double*[4];

	double pcKsi[4] = { -1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3) };
	double pcEta[4] = { -1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3) };

	Elem4() {
		for (int i = 0; i < 4; i++) {
			tabKsi[i] = new double[4];
			tabEta[i] = new double[4];
		}
	}

	void printElem() {
		cout << "\ndN/dKsi\n";
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << tabKsi[i][j] << " ";
			}
			cout << endl;
		}

		cout << "\ndN/dEta\n";
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << tabEta[i][j] << " ";
			}
			cout << endl;
		}
	}
};

//structure for derivatives with 3points
struct Elem9 {
	double** tabKsi = new double* [9];
	double** tabEta = new double* [9];

	double pcKsi[9] = { -sqrt(0.6), 0, sqrt(0.6), -sqrt(0.6), 0, sqrt(0.6), -sqrt(0.6), 0 ,sqrt(0.6) };
	double pcEta[9] = { -sqrt(0.6), -sqrt(0.6), -sqrt(0.6), 0, 0, 0, sqrt(0.6), sqrt(0.6) , sqrt(0.6) };
	Elem9() {
		for (int i = 0; i < 9; i++) {
			tabKsi[i] = new double[4];
			tabEta[i] = new double[4];
		}
	}

	void printElem() {
		cout << "\ndN/dKsi\n";
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++) {
				cout << tabKsi[i][j] << " ";
			}
			cout << endl;
		}

		cout << "\ndN/dEta\n";
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++) {
				cout << tabEta[i][j] << " ";
			}
			cout << endl;
		}
	}

};

//structure for derivatives with 4points
struct Elem16 {
	double** tabKsi = new double* [16];
	double** tabEta = new double* [16];

	double pcKsi[16] = { - sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), -sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), -sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)),
		-sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), -sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),  -sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),
		-sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)) , -sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)) , sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)) ,
		sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)) , sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)) , sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),
		sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)) , sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)) ,
		sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)) };
	double pcEta[16] = { -sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), -sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)), sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),
		sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), -sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), -sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),
		sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)), sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)),-sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)),
		-sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)), sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)), sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)),
		-sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), -sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)), sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),
		sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)) };

	Elem16() {
		for (int i = 0; i < 16; i++) {
			tabKsi[i] = new double[4];
			tabEta[i] = new double[4];
		}
	}

	void printElem() {
		cout << "\ndN/dKsi\n";
		for (int i = 0; i < 16; i++) {
			for (int j = 0; j < 4; j++) {
				cout << tabKsi[i][j] << " ";
			}
			cout << endl;
		}

		cout << "\ndN/dEta\n";
		for (int i = 0; i < 16; i++) {
			for (int j = 0; j < 4; j++) {
				cout << tabEta[i][j] << " ";
			}
			cout << endl;
		}
	}
};

//scheme for caculating derivatives with 2,3,4 points
void derivativesScheme(int i, Elem4* elem4, Elem9* elem9, Elem16* elem16, int numOfPoints) {

	if (numOfPoints == 2) {
		nKsi[0] = -0.25 * (1 - elem4->pcEta[i]);
		nKsi[1] = 0.25 * (1 - elem4->pcEta[i]);
		nKsi[2] = 0.25 * (1 + elem4->pcEta[i]);
		nKsi[3] = -0.25 * (1 + elem4->pcEta[i]);

		nEta[0] = -0.25 * (1 - elem4->pcKsi[i]);
		nEta[1] = -0.25 * (1 + elem4->pcKsi[i]);
		nEta[2] = 0.25 * (1 + elem4->pcKsi[i]);
		nEta[3] = 0.25 * (1 - elem4->pcKsi[i]);
	}
	else if (numOfPoints == 3) {
		nKsi[0] = -0.25 * (1 - elem9->pcEta[i]);
		nKsi[1] = 0.25 * (1 - elem9->pcEta[i]);
		nKsi[2] = 0.25 * (1 + elem9->pcEta[i]);
		nKsi[3] = -0.25 * (1 + elem9->pcEta[i]);

		nEta[0] = -0.25 * (1 - elem9->pcKsi[i]);
		nEta[1] = -0.25 * (1 + elem9->pcKsi[i]);
		nEta[2] = 0.25 * (1 + elem9->pcKsi[i]);
		nEta[3] = 0.25 * (1 - elem9->pcKsi[i]);
	}
	else if (numOfPoints == 4) {
		nKsi[0] = -0.25 * (1 - elem16->pcEta[i]);
		nKsi[1] = 0.25 * (1 - elem16->pcEta[i]);
		nKsi[2] = 0.25 * (1 + elem16->pcEta[i]);
		nKsi[3] = -0.25 * (1 + elem16->pcEta[i]);

		nEta[0] = -0.25 * (1 - elem16->pcKsi[i]);
		nEta[1] = -0.25 * (1 + elem16->pcKsi[i]);
		nEta[2] = 0.25 * (1 + elem16->pcKsi[i]);
		nEta[3] = 0.25 * (1 - elem16->pcKsi[i]);
	}
}

//reading the file
void readFile(GlobalData* globaldata, Grid* grid) {
	string line;
	//int n, number;
	fstream dataFile;
	dataFile.open("Test1_4_4.txt", ios::in);

	if (dataFile.is_open()) {
		string name;

		//reading single values
		dataFile >> name >> globaldata->t >> name >> globaldata->st >> name >> globaldata->lambda >> name >> globaldata->alpha >> name >>
			globaldata->tot >> name >> globaldata->it >> name >> globaldata->d >> name >> globaldata->cp >> name >> name >> grid->nN >>
			name >> name >> grid->nE >> name;
		
		/*cout << "Simulation Time: " << globaldata->t << endl;
		cout << "Simulation Step Time: " << globaldata->st << endl;
		cout << "Conductivity: " << globaldata->lambda << endl;
		cout << "Alpha: " << globaldata->alpha << endl;
		cout << "Tot: " << globaldata->tot << endl;
		cout << "Initial Temp: " << globaldata->it << endl;
		cout << "Density: " << globaldata->d << endl;
		cout << "Specific Heat: " << globaldata->cp << endl;
		cout << "Nodes Numer: " << grid->nN << endl;
		cout << "Elements Number: " << grid->nE << endl;
		cout << "\nNodes: " << endl;*/
		
		//reading nodes into vector
		Node node;
		for (int i = 0; i < grid->nN; i++) {
			
			float x, y, n;
			dataFile >> n >> name >> x >> name >> y;
			node.x = x;
			node.y = y;
			
			//printf("%d %.10f  %.10f\n",i, x, y);
			grid->nodes.push_back(node);
		}
		getline(dataFile, line);
		getline(dataFile, line);
		
		//reading elements into vector
		Element element;
		for (int i = 0; i < grid->nE; i++) {
			int n, el[4] = { 0,0,0,0 };

			dataFile >> n >> name >> el[0] >> name >> el[1] >> name >> el[2] >> name >> el[3];

			element.ID[0] = el[0];
			element.ID[1] = el[1];
			element.ID[2] = el[2];
			element.ID[3] = el[3];
			
			//printf("%d %d %d %d %d\n", i, el[0], el[1], el[2], el[3]);
			grid->elements.push_back(element);
		}
		getline(dataFile, line);
		getline(dataFile, line);
		//cout << line;
	}
	dataFile.close();
}

//functions to read
double function1D(double x) {
	return 2 * x * x + 3 * x - 8;
}
double function2D(double ksi, double eta) {
	return -5 * ksi * ksi * eta + 2 * ksi * eta * eta + 10;
}
double func(double x) {
	return 3 * x * x - 6 * x + 1;
}

//integration without bounds, but with 2 dimentions
double integration(Integration* scheme, int numOfPoints, int numOfDimention) {
	double result = 0;

	if (numOfDimention == 1) {
		if (numOfPoints == 2) {
			for (int i = 0; i < numOfPoints; i++) result += function1D((*scheme).nodes2[i]) * (*scheme).weights2[i];
		}
		if (numOfPoints == 3) {
			for (int i = 0; i < numOfPoints; i++) result += function1D((*scheme).nodes3[i]) * (*scheme).weights3[i]; 
		}
	}
	if (numOfDimention == 2) {
		if (numOfPoints == 2) {
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < numOfPoints; j++) {
					result += function2D((*scheme).nodes2[i], (*scheme).nodes2[j]) * (*scheme).weights2[i] * (*scheme).weights2[j];
				}
			}
		}
		if (numOfPoints == 3) {
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < numOfPoints; j++) {
					result += function2D((*scheme).nodes3[i], (*scheme).nodes3[j]) * (*scheme).weights3[i] * (*scheme).weights3[j];
				}
			}
		}
	}
	return result;
}

//integration with bounds and 1 dimention
double integrationBounds(Integration* scheme, int numOfPoints, double x1, double x2) {
	double result = 0;
	if (numOfPoints == 2) {
		for (int i = 0; i < numOfPoints; i++) result += func(((1 - (*scheme).nodes2[i]) / 2) * x1 +
			(((*scheme).nodes2[i] + 1) / 2) * x2) * ((x2 - x1) / 2) * (*scheme).weights2[i];
	}
	if (numOfPoints == 3) {
		for (int i = 0; i < numOfPoints; i++) result += func(((1 - (*scheme).nodes3[i]) / 2) * x1 +
			(((*scheme).nodes3[i] + 1) / 2) * x2) * ((x2 - x1) / 2) * (*scheme).weights3[i];
	}
	return result;
}

//counting derivatives
void derivativesElem(Elem4* elem4, Elem9* elem9,Elem16* elem16, void(*derivativesScheme)(int,Elem4*,Elem9*,Elem16*, int), int numOfPoints) {

	if (numOfPoints == 2) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				derivativesScheme(i, elem4, elem9, elem16, 2);
				elem4->tabKsi[i][j] = nKsi[j];
			}
		}

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				derivativesScheme(i, elem4, elem9, elem16, 2);
				elem4->tabEta[i][j] = nEta[j];
			}
		}

		//elem4->printElem();
	}
	else if (numOfPoints == 3) {
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++) {
				derivativesScheme(i, elem4, elem9, elem16, 3);
				elem9->tabKsi[i][j] = nKsi[j];
			}
		}

		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++) {
				derivativesScheme(i, elem4, elem9, elem16, 3);
				elem9->tabEta[i][j] = nEta[j];
			}
		}
		//elem9->printElem();
	}
	else if (numOfPoints == 4) {
		for (int i = 0; i < 16; i++) {
			for (int j = 0; j < 4; j++) {
				derivativesScheme(i, elem4, elem9, elem16, 4);
				elem16->tabKsi[i][j] = nKsi[j];
			}
		}

		for (int i = 0; i < 16; i++) {
			for (int j = 0; j < 4; j++) {
				derivativesScheme(i, elem4, elem9, elem16, 4);
				elem16->tabEta[i][j] = nEta[j];
			}
		}
		//elem16->printElem();
	}
	else cout << "Wrong number of Points" << endl;
}

//counting Matrix H
void matrixH(int numOfPoints, Elem4* elem4, Elem9 *elem9, Elem16* elem16, Grid* grid, GlobalData* globaldata) {
	derivativesElem(elem4, elem9, elem16, &derivativesScheme, numOfPoints);
	numOfPoints *= numOfPoints;
	for (int i = 0; i < grid->nE; i++) {
		double x[4], y[4];
		int idN[4];
		for (int j = 0; j < 4; j++) {
			idN[j] = grid->elements[i].ID[j];
			x[j] = grid->nodes[idN[j] - 1].x;
			y[j] = grid->nodes[idN[j] - 1].y;
		}

		double** jacobian = new double* [numOfPoints];
		double** tabX = new double* [numOfPoints];
		double** tabY = new double* [numOfPoints];
		double** matrixH = new double* [numOfPoints];
		for (int i = 0; i < numOfPoints; i++) {
			jacobian[i] = new double[4];
			tabX[i] = new double[4];
			tabY[i] = new double[4];
			matrixH[i] = new double[4];
		}

		double* detJacobian = new double[numOfPoints];
		double* detJacobianMinus = new double[numOfPoints];

		for (int i = 0; i < numOfPoints; i++) {
			for (int j = 0; j < 4; j++) {
				jacobian[i][j] = 0;
				matrixH[i][j] = 0;
				detJacobian[i] = 0;
				detJacobianMinus[i] = 0;
			}
		}

		//jacobian
		if (numOfPoints == 4) {
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					jacobian[i][0] += elem4->tabKsi[i][j] * x[j];
					jacobian[i][1] += elem4->tabKsi[i][j] * y[j];
					jacobian[i][2] += elem4->tabEta[i][j] * x[j];
					jacobian[i][3] += elem4->tabEta[i][j] * y[j];
				}
			}
		}
		else if (numOfPoints == 9) {
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					jacobian[i][0] += elem9->tabKsi[i][j] * x[j];
					jacobian[i][1] += elem9->tabKsi[i][j] * y[j];
					jacobian[i][2] += elem9->tabEta[i][j] * x[j];
					jacobian[i][3] += elem9->tabEta[i][j] * y[j];
				}
			}
		}
		else if (numOfPoints == 16) {
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					jacobian[i][0] += elem16->tabKsi[i][j] * x[j];
					jacobian[i][1] += elem16->tabKsi[i][j] * y[j];
					jacobian[i][2] += elem16->tabEta[i][j] * x[j];
					jacobian[i][3] += elem16->tabEta[i][j] * y[j];
				}
			}
		}

		// 1/detJ
		for (int i = 0; i < numOfPoints; i++) {
			detJacobianMinus[i] = (jacobian[i][0] * jacobian[i][3] - jacobian[i][1] * jacobian[i][2]);
			detJacobian[i] = 1 / detJacobianMinus[i];
		}
		for (int i = 0; i < numOfPoints; i++) {
			for (int j = 0; j < 4; j++) {
				jacobian[i][j] *= detJacobian[i];
			}
		}

		//array with dX and dY
		if (numOfPoints == 4) {
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					tabX[i][j] = jacobian[i][0] * elem4->tabKsi[i][j] + jacobian[i][1] * elem4->tabEta[i][j];
				}
			}
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					tabY[i][j] = jacobian[i][2] * elem4->tabKsi[i][j] + jacobian[i][3] * elem4->tabEta[i][j];
				}
			}
		}
		else if (numOfPoints == 9) {
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					tabX[i][j] = jacobian[i][0] * elem9->tabKsi[i][j] + jacobian[i][1] * elem9->tabEta[i][j];
				}
			}
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					tabY[i][j] = jacobian[i][2] * elem9->tabKsi[i][j] + jacobian[i][3] * elem9->tabEta[i][j];
				}
			}
		}
		else if (numOfPoints == 16) {
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					tabX[i][j] = jacobian[i][0] * elem16->tabKsi[i][j] + jacobian[i][1] * elem16->tabEta[i][j];
				}
			}
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					tabY[i][j] = jacobian[i][2] * elem16->tabKsi[i][j] + jacobian[i][3] * elem16->tabEta[i][j];
				}
			}
		}
		
		//final Matrix H
		for (int i = 0; i < numOfPoints; i++) {
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					if (numOfPoints == 4) {
						matrixH[j][k] += (globaldata->lambda * (tabX[i][k] * tabX[i][j] + tabY[i][k] * tabY[i][j]) * detJacobianMinus[i]);
					}
					else if (numOfPoints == 9) {
						matrixH[j][k] += (globaldata->lambda * (tabX[i][k] * tabX[i][j] + tabY[i][k] * tabY[i][j]) * detJacobianMinus[i])
							* w1for3[i] * w2for3[i];
					}
					else if (numOfPoints == 16) {
						matrixH[j][k] += (globaldata->lambda * (tabX[i][k] * tabX[i][j] + tabY[i][k] * tabY[i][j]) * detJacobianMinus[i])
							* w1for4[i] * w2for4[i];
					}
				}
			}
		}

		//print
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << matrixH[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}

}

void aggregation() {
	double** globalH = new double* [16];
	for (int i = 0; i < 16; i++) globalH[i] = new double[16];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			
		}
	}

}

int main() {
	//structures and lists go brrrrrr
	Grid grid;
	GlobalData globaldata;
	Integration scheme;
	list<Element> listOfElements;
	Elem4 elem4; Elem9 elem9; Elem16 elem16;

	//reading the file
	readFile(&globaldata, &grid);

	//results of integration
	//solution integration(&scheme, numOfPoints, numOfDimention)

	//results of derivatives
	//solution derivativesElem(&elem4, &elem9, &elem16, &derivativesScheme, numOfPoints)

	//results of MatrixH
	//solution matrixH(numofPoints, &listOfNodes, &elem4, &elem9, &elem16);
	
	//matrixH with the data given
	matrixH(4, &elem4, &elem9, &elem16, &grid, &globaldata);
	
	return 0;
}