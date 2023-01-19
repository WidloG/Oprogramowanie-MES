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
const double eps = 1e-12;
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

	void setBC(bool BC) {
		bc = BC;
	}
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
	vector <Node> BC;
	
	void printGrid() {
		for (int i = 0; i < nodes.size(); i++) {
			printf("%.10f %.10f  %d\n", nodes[i].x, nodes[i].y, nodes[i].bc);
		}
	}
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

	double weights4[4] = { (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };
};

//structure for derivatives with 2points
struct Elem4 {
	double** tabKsi = new double*[4];
	double** tabEta = new double*[4];
	double** nKsiEta = new double* [4];

	double pcKsi[4] = { -1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3) };
	double pcEta[4] = { -1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3) };

	double pcKsiBC[8] = { -1 / sqrt(3), 1 / sqrt(3), 1, 1, 1 / sqrt(3), -1 / sqrt(3), -1, -1 };
	double pcEtaBC[8] = { -1, -1, -1 / sqrt(3), 1 / sqrt(3), 1,1,1 / sqrt(3),- 1 / sqrt(3) };

	double N[8][4];

	Elem4() {
		for (int i = 0; i < 4; i++) {
			tabKsi[i] = new double[4];
			tabEta[i] = new double[4];
			nKsiEta[i] = new double[4];
		}
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 4; j++) N[i][j] = 0.0;
		}

		for (int i = 0; i < 4; i++) {
			nKsiEta[i][0] = 0.25 * (1 - pcKsi[i]) * (1 - pcEta[i]);
			nKsiEta[i][1] = 0.25 * (1 + pcKsi[i]) * (1 - pcEta[i]);
			nKsiEta[i][2] = 0.25 * (1 + pcKsi[i]) * (1 + pcEta[i]);
			nKsiEta[i][3] = 0.25 * (1 - pcKsi[i]) * (1 + pcEta[i]);
			
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
	double** nKsiEta = new double* [9];

	double pcKsi[9] = { -sqrt(0.6), 0, sqrt(0.6), -sqrt(0.6), 0, sqrt(0.6), -sqrt(0.6), 0, sqrt(0.6) };
	double pcEta[9] = { -sqrt(0.6), -sqrt(0.6), -sqrt(0.6), 0, 0, 0, sqrt(0.6), sqrt(0.6), sqrt(0.6) };

	double pcKsiBC[12] = { -sqrt(0.6), 0, sqrt(0.6), 1, 1, 1, sqrt(0.6), 0 , -sqrt(0.6), -1 ,-1 ,-1 };
	double pcEtaBC[12] = { -1, -1, -1, -sqrt(0.6), 0, sqrt(0.6), 1, 1, 1, sqrt(0.6), 0 , -sqrt(0.6) };

	double N[12][4];

	Elem9() {
		for (int i = 0; i < 9; i++) {
			tabKsi[i] = new double[4];
			tabEta[i] = new double[4];
			nKsiEta[i] = new double[4];
		}
		for (int i = 0; i < 12; i++) {
			for (int j = 0; j < 4; j++) N[i][j] = 0.0;
		}

		for (int i = 0; i < 9; i++) {
			nKsiEta[i][0] = 0.25 * (1 - pcKsi[i]) * (1 - pcEta[i]);
			nKsiEta[i][1] = 0.25 * (1 + pcKsi[i]) * (1 - pcEta[i]);
			nKsiEta[i][2] = 0.25 * (1 + pcKsi[i]) * (1 + pcEta[i]);
			nKsiEta[i][3] = 0.25 * (1 - pcKsi[i]) * (1 + pcEta[i]);
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
	double** nKsiEta = new double* [16];

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

	double pcKsiBC[16] = { -sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), -sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)), sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)), sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)),
							1,1,1,1,  sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)),  sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),  -sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),  -sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), -1,-1,-1,-1 };
	double pcEtaBC[16] = { -1, -1, -1, -1,  -sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)),  -sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),  sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),  sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)),
							1,1,1,1, sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)), sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)), -sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)) , -sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)) };

	double N[16][4];

	Elem16() {
		for (int i = 0; i < 16; i++) {
			tabKsi[i] = new double[4];
			tabEta[i] = new double[4];
			nKsiEta[i] = new double[4];
		}
		for (int i = 0; i < 16; i++) {
			for (int j = 0; j < 4; j++) N[i][j] = 0.0;
		}

		for (int i = 0; i < 16; i++) {
			nKsiEta[i][0] = 0.25 * (1 - pcKsi[i]) * (1 - pcEta[i]);
			nKsiEta[i][1] = 0.25 * (1 + pcKsi[i]) * (1 - pcEta[i]);
			nKsiEta[i][2] = 0.25 * (1 + pcKsi[i]) * (1 + pcEta[i]);
			nKsiEta[i][3] = 0.25 * (1 - pcKsi[i]) * (1 + pcEta[i]);
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

//counting bc on borders + vectorP
struct Borders {
	int numOfWall = 0; int numOfPoints = 0;
	int alpha = 0;
	int tot = 0;

	Elem4 elem4; Elem9 elem9; Elem16 elem16;
	Integration scheme;
	
	double detJ[4] = { 0,0,0,0 };
	double** side1 = new double* [4];
	double** side2 = new double* [4];
	double** side3 = new double* [4];
	double** side4 = new double* [4];
	double* vectorP1{ new double [4] {} }; double* vectorP2{ new double [4] {} }; double* vectorP3{ new double [4] {} }; double* vectorP4{ new double [4] {} };
	
	double x[4], y[4];

	Borders(int numOfWall, int numOfPoints, double x[4], double y[4], int alpha, int tot) {
		this->numOfWall = numOfWall;
		this->numOfPoints = numOfPoints;
		for (int i = 0; i < 4; i++) {
			this->x[i] = x[i];
			this->y[i] = y[i];
		}
		this->alpha = alpha;
		this->tot = tot;
		for (int i = 0; i < 4; i++) {
			side1[i] = new double[4];
			side2[i] = new double[4];
			side3[i] = new double[4];
			side4[i] = new double[4];
			for (int j = 0; j < 4; j++) {
				side1[i][j] = 0.0; side2[i][j] = 0.0; side3[i][j] = 0.0; side4[i][j] = 0.0;
				vectorP1[i] = 0; vectorP2[i] = 0; vectorP3[i] = 0; vectorP4[i] = 0;
			}
		}
	}

	double **borders(){
		if (numOfPoints == 4) {
			for (int i = 0; i < 8; i++) {
				elem4.N[i][0] = 0.25 * (1 - elem4.pcKsiBC[i]) * (1 - elem4.pcEtaBC[i]);
				elem4.N[i][1] = 0.25 * (1 + elem4.pcKsiBC[i]) * (1 - elem4.pcEtaBC[i]);
				elem4.N[i][2] = 0.25 * (1 + elem4.pcKsiBC[i]) * (1 + elem4.pcEtaBC[i]);
				elem4.N[i][3] = 0.25 * (1 - elem4.pcKsiBC[i]) * (1 + elem4.pcEtaBC[i]);
			}

			if (numOfWall == 1) {
				detJ[0] = (sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side1[i][j] += alpha * (1 * (elem4.N[0][j] * elem4.N[0][i]) +
							1 * (elem4.N[1][j] * elem4.N[1][i])) * detJ[0];
					}
				}
				return side1;
			}
			else if (numOfWall == 2) {
				detJ[1] = (sqrt((x[2] - x[1]) * (x[2] - x[1]) + (y[2] - y[1]) * (y[2] - y[1]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side2[i][j] += alpha * (1 * (elem4.N[2][j] * elem4.N[2][i]) +
							1 * (elem4.N[3][j] * elem4.N[3][i])) * detJ[1];
					}
				}
				return side2;
			}
			else if (numOfWall == 3) {
				detJ[2] = (sqrt((x[3] - x[2]) * (x[3] - x[2]) + (y[3] - y[2]) * (y[3] - y[2]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side3[i][j] += alpha * (1 * (elem4.N[4][j] * elem4.N[4][i]) +
							1 * (elem4.N[5][j] * elem4.N[5][i])) * detJ[2];
					}
				}
				return side3;
			}
			else if (numOfWall == 4) {
				detJ[3] = (sqrt((x[0] - x[3]) * (x[0] - x[3]) + (y[0] - y[3]) * (y[0] - y[3]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side4[i][j] += alpha * (1 * (elem4.N[6][j] * elem4.N[6][i]) +
							1 * (elem4.N[7][j] * elem4.N[7][i])) * detJ[3];
					}
				}
				return side4;
			}

		}
		if (numOfPoints == 9) {
			for (int i = 0; i < 12; i++) {
				elem9.N[i][0] = 0.25 * (1 - elem9.pcKsiBC[i]) * (1 - elem9.pcEtaBC[i]);
				elem9.N[i][1] = 0.25 * (1 + elem9.pcKsiBC[i]) * (1 - elem9.pcEtaBC[i]);
				elem9.N[i][2] = 0.25 * (1 + elem9.pcKsiBC[i]) * (1 + elem9.pcEtaBC[i]);
				elem9.N[i][3] = 0.25 * (1 - elem9.pcKsiBC[i]) * (1 + elem9.pcEtaBC[i]);
			}

			if (numOfWall == 1) {
				detJ[0] = (sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side1[i][j] += alpha * ((scheme.weights3[0] * (elem9.N[0][j] * elem9.N[0][i])) +
							(scheme.weights3[1] * (elem9.N[1][j] * elem9.N[1][i])) +
							(scheme.weights3[2] * (elem9.N[2][j] * elem9.N[2][i]))) * detJ[0];
					}
				}
				return side1;
			}
			if (numOfWall == 2) {
				detJ[1] = (sqrt((x[2] - x[1]) * (x[2] - x[1]) + (y[2] - y[1]) * (y[2] - y[1]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side2[i][j] += alpha * ((scheme.weights3[0] * (elem9.N[3][j] * elem9.N[3][i])) +
							(scheme.weights3[1] * (elem9.N[4][j] * elem9.N[4][i])) +
							(scheme.weights3[2] * (elem9.N[5][j] * elem9.N[5][i]))) * detJ[1];
					}
				}
				return side2;
			}
			if (numOfWall == 3) {
				detJ[2] = (sqrt((x[3] - x[2]) * (x[3] - x[2]) + (y[3] - y[2]) * (y[3] - y[2]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side3[i][j] += alpha * ((scheme.weights3[0] * (elem9.N[6][j] * elem9.N[6][i])) +
							(scheme.weights3[1] * (elem9.N[7][j] * elem9.N[7][i])) +
							(scheme.weights3[2] * (elem9.N[8][j] * elem9.N[8][i]))) * detJ[2];
					}
				}
				return side3;
			}
			if (numOfWall == 4) {
				detJ[3] = (sqrt((x[0] - x[3]) * (x[0] - x[3]) + (y[0] - y[3]) * (y[0] - y[3]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side4[i][j] += alpha * ((scheme.weights3[0] * (elem9.N[9][j] * elem9.N[9][i])) +
							(scheme.weights3[1] * (elem9.N[10][j] * elem9.N[10][i])) +
							(scheme.weights3[2] * (elem9.N[11][j] * elem9.N[11][i]))) * detJ[3];
					}
				}
				return side4;
			}
		}
		if (numOfPoints == 16) {

			for (int i = 0; i < 16; i++) {
				elem16.N[i][0] = 0.25 * (1 - elem16.pcKsiBC[i]) * (1 - elem16.pcEtaBC[i]);
				elem16.N[i][1] = 0.25 * (1 + elem16.pcKsiBC[i]) * (1 - elem16.pcEtaBC[i]);
				elem16.N[i][2] = 0.25 * (1 + elem16.pcKsiBC[i]) * (1 + elem16.pcEtaBC[i]);
				elem16.N[i][3] = 0.25 * (1 - elem16.pcKsiBC[i]) * (1 + elem16.pcEtaBC[i]);
			}

			if (numOfWall == 1) {
				detJ[0] = (sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side1[i][j] += alpha * ((scheme.weights4[0] * (elem16.N[0][j] * elem16.N[0][i])) +
							(scheme.weights4[1] * (elem16.N[1][j] * elem16.N[1][i])) +
							(scheme.weights4[2] * (elem16.N[2][j] * elem16.N[2][i])) +
							(scheme.weights4[3] * (elem16.N[3][j] * elem16.N[3][i]))) * detJ[0];
					}
				}
				return side1;
			}
			if (numOfWall == 2) {
				detJ[1] = (sqrt((x[2] - x[1]) * (x[2] - x[1]) + (y[2] - y[1]) * (y[2] - y[1]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side2[i][j] += alpha * ((scheme.weights4[0] * (elem16.N[4][j] * elem16.N[4][i])) +
							(scheme.weights4[1] * (elem16.N[5][j] * elem16.N[5][i])) +
							(scheme.weights4[2] * (elem16.N[6][j] * elem16.N[6][i])) +
							(scheme.weights4[3] * (elem16.N[7][j] * elem16.N[7][i]))) * detJ[1];
					}
				}
				return side2;
			}
			if (numOfWall == 3) {
				detJ[2] = (sqrt((x[3] - x[2]) * (x[3] - x[2]) + (y[3] - y[2]) * (y[3] - y[2]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side3[i][j] += alpha * ((scheme.weights4[0] * (elem16.N[8][j] * elem16.N[8][i])) +
							(scheme.weights4[1] * (elem16.N[9][j] * elem16.N[9][i])) +
							(scheme.weights4[2] * (elem16.N[10][j] * elem16.N[10][i])) +
							(scheme.weights4[3] * (elem16.N[11][j] * elem16.N[11][i]))) * detJ[2];
					}
				}
				return side3;
			}
			if (numOfWall == 4) {
				detJ[3] = (sqrt((x[0] - x[3]) * (x[0] - x[3]) + (y[0] - y[3]) * (y[0] - y[3]))) / 2;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						side4[i][j] += alpha * ((scheme.weights4[0] * (elem16.N[12][j] * elem16.N[12][i])) +
							(scheme.weights4[1] * (elem16.N[13][j] * elem16.N[13][i])) +
							(scheme.weights4[2] * (elem16.N[14][j] * elem16.N[14][i])) +
							(scheme.weights4[3] * (elem16.N[15][j] * elem16.N[15][i]))) * detJ[3];
					}
				}

				return side4;
			}
		}
		return 0;
	}
	double* vectorP() {
		if (numOfPoints == 4) {
			for (int i = 0; i < 8; i++) {
				elem4.N[i][0] = 0.25 * (1 - elem4.pcKsiBC[i]) * (1 - elem4.pcEtaBC[i]);
				elem4.N[i][1] = 0.25 * (1 + elem4.pcKsiBC[i]) * (1 - elem4.pcEtaBC[i]);
				elem4.N[i][2] = 0.25 * (1 + elem4.pcKsiBC[i]) * (1 + elem4.pcEtaBC[i]);
				elem4.N[i][3] = 0.25 * (1 - elem4.pcKsiBC[i]) * (1 + elem4.pcEtaBC[i]);
			}

			if (numOfWall == 1) {
				detJ[0] = (sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP1[i] = alpha * (1 * (elem4.N[0][i] * tot) + 1 * (elem4.N[1][i] * tot)) * detJ[0];
				}
				return vectorP1;
			}
			if (numOfWall == 2) {
				detJ[1] = (sqrt((x[2] - x[1]) * (x[2] - x[1]) + (y[2] - y[1]) * (y[2] - y[1]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP2[i] = alpha * (1 * (elem4.N[2][i] * tot) + 1 * (elem4.N[3][i] * tot)) * detJ[1];
				}
				return vectorP2;
			}
			if (numOfWall == 3) {
				detJ[2] = (sqrt((x[3] - x[2]) * (x[3] - x[2]) + (y[3] - y[2]) * (y[3] - y[2]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP3[i] = alpha * (1 * (elem4.N[4][i] * tot) + 1 * (elem4.N[5][i] * tot)) * detJ[2];
				}
				return vectorP3;
			}
			if (numOfWall == 4) {
				detJ[3] = (sqrt((x[0] - x[3]) * (x[0] - x[3]) + (y[0] - y[3]) * (y[0] - y[3]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP4[i] = alpha * (1 * (elem4.N[6][i] * tot) + 1 * (elem4.N[7][i] * tot)) * detJ[3];
				}
				return vectorP4;
			}
		}
		if (numOfPoints == 9) {
			for (int i = 0; i < 12; i++) {
				elem9.N[i][0] = 0.25 * (1 - elem9.pcKsiBC[i]) * (1 - elem9.pcEtaBC[i]);
				elem9.N[i][1] = 0.25 * (1 + elem9.pcKsiBC[i]) * (1 - elem9.pcEtaBC[i]);
				elem9.N[i][2] = 0.25 * (1 + elem9.pcKsiBC[i]) * (1 + elem9.pcEtaBC[i]);
				elem9.N[i][3] = 0.25 * (1 - elem9.pcKsiBC[i]) * (1 + elem9.pcEtaBC[i]);
			}

			if (numOfWall == 1) {
				detJ[0] = (sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP1[i] = alpha * (scheme.weights3[0] * (elem9.N[0][i] * tot) + scheme.weights3[1] * (elem9.N[1][i] * tot) + scheme.weights3[2] * (elem9.N[2][i] * tot)) * detJ[0];
				}
				return vectorP1;
			}
			if (numOfWall == 2) {
				detJ[1] = (sqrt((x[2] - x[1]) * (x[2] - x[1]) + (y[2] - y[1]) * (y[2] - y[1]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP2[i] = alpha * (scheme.weights3[0] * (elem9.N[3][i] * tot) + scheme.weights3[1] * (elem9.N[4][i] * tot) + scheme.weights3[2] * (elem9.N[5][i] * tot)) * detJ[1];
				}
				return vectorP2;
			}
			if (numOfWall == 3) {
				detJ[2] = (sqrt((x[3] - x[2]) * (x[3] - x[2]) + (y[3] - y[2]) * (y[3] - y[2]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP3[i] = alpha * (scheme.weights3[0] * (elem9.N[6][i] * tot) + scheme.weights3[1] * (elem9.N[7][i] * tot) + scheme.weights3[2] * (elem9.N[8][i] * tot)) * detJ[2];
				}
				return vectorP3;
			}
			if (numOfWall == 4) {
				detJ[3] = (sqrt((x[0] - x[3]) * (x[0] - x[3]) + (y[0] - y[3]) * (y[0] - y[3]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP4[i] = alpha * (scheme.weights3[0] * (elem9.N[9][i] * tot) + scheme.weights3[1] * (elem9.N[10][i] * tot) + scheme.weights3[2] * (elem9.N[11][i] * tot)) * detJ[3];
				}
				return vectorP4;
			}
		}
		if (numOfPoints == 16) {
			for (int i = 0; i < 16; i++) {
				elem16.N[i][0] = 0.25 * (1 - elem16.pcKsiBC[i]) * (1 - elem16.pcEtaBC[i]);
				elem16.N[i][1] = 0.25 * (1 + elem16.pcKsiBC[i]) * (1 - elem16.pcEtaBC[i]);
				elem16.N[i][2] = 0.25 * (1 + elem16.pcKsiBC[i]) * (1 + elem16.pcEtaBC[i]);
				elem16.N[i][3] = 0.25 * (1 - elem16.pcKsiBC[i]) * (1 + elem16.pcEtaBC[i]);
			}

			if (numOfWall == 1) {
				detJ[0] = (sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP1[i] = alpha * (scheme.weights4[0] * (elem16.N[0][i] * tot) + scheme.weights4[1] * (elem16.N[1][i] * tot) +
						scheme.weights4[2] * (elem16.N[2][i] * tot) + scheme.weights4[3] * (elem16.N[3][i] * tot)) * detJ[0];
				}
				return vectorP1;
			}
			if (numOfWall == 2) {
				detJ[1] = (sqrt((x[2] - x[1]) * (x[2] - x[1]) + (y[2] - y[1]) * (y[2] - y[1]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP2[i] = alpha * (scheme.weights4[0] * (elem16.N[4][i] * tot) + scheme.weights4[1] * (elem16.N[5][i] * tot) +
						scheme.weights4[2] * (elem16.N[6][i] * tot) + scheme.weights4[3] * (elem16.N[7][i] * tot)) * detJ[1];
				}
				return vectorP2;
			}
			if (numOfWall == 3) {
				detJ[2] = (sqrt((x[3] - x[2]) * (x[3] - x[2]) + (y[3] - y[2]) * (y[3] - y[2]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP3[i] = alpha * (scheme.weights4[0] * (elem16.N[8][i] * tot) + scheme.weights4[1] * (elem16.N[9][i] * tot) +
						scheme.weights4[2] * (elem16.N[10][i] * tot) + scheme.weights4[3] * (elem16.N[11][i] * tot)) * detJ[2];
				}
				return vectorP3;
			}
			if (numOfWall == 4) {
				detJ[3] = (sqrt((x[0] - x[3]) * (x[0] - x[3]) + (y[0] - y[3]) * (y[0] - y[3]))) / 2;
				for (int i = 0; i < 4; i++) {
					vectorP4[i] = alpha * (scheme.weights4[0] * (elem16.N[12][i] * tot) + scheme.weights4[1] * (elem16.N[13][i] * tot) +
						scheme.weights4[2] * (elem16.N[14][i] * tot) + scheme.weights4[3] * (elem16.N[15][i] * tot)) * detJ[3];
				}
				return vectorP4;
			}
		}
		return 0;
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
	fstream dataFile;
	dataFile.open("Test3_31_31.txt", ios::in);

	if (dataFile.is_open()) {
		string name;

		//reading single values
		dataFile >> name >> globaldata->t >> name >> globaldata->st >> name >> globaldata->lambda >> name >> globaldata->alpha >> name >>
			globaldata->tot >> name >> globaldata->it >> name >> globaldata->d >> name >> globaldata->cp >> name >> name >> grid->nN >>
			name >> name >> grid->nE >> name;
		
		//reading nodes into vector
		Node node{};
		float x, y, n;
		for (int i = 0; i < grid->nN; i++) {
			
			dataFile >> n >> name >> x >> name >> y;
			node.x = x;
			node.y = y;
			
			//printf("%d %.10f  %.10f\n",i, x, y);
			grid->nodes.push_back(node);
		}
		getline(dataFile, line);
		getline(dataFile, line);
		
		//reading elements into vector
		Element element{};
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
		while (!dataFile.eof()) {
			int cond;
			dataFile >> cond >> name;
			grid->nodes[static_cast<std::vector<Node, std::allocator<Node>>::size_type>(cond) - 1].setBC(true);
			grid->BC.push_back(node);
		}
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
	//formula: f(ksi,eta)*wi*wj
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
	//formula: f(1 - N/2)*x1 + (((N+1)/2)*x2) * detJ * weight
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
	}
	else cout << "Wrong number of Points" << endl;
}

//Gauss method
void gauss(int n, double** AB, double* X){
	for (int j = 0; j < n - 1; j++){
		for (int i = j + 1; i < n; i++){
			double temp = AB[i][j] / AB[j][j];

			for (int k = 0; k < n + 1; k++)
				AB[i][k] -= AB[j][k] * temp;
		}
	}

	for (int i = n - 1; i >= 0; i--){
		double sum = 0;
		for (int j = i + 1; j < n; j++)
			sum += AB[i][j] * X[j];
		X[i] = (AB[i][n] - sum) / AB[i][i];
	}
}

//counting Matrix H and aggregation
void matrixH(int numOfPoints, Elem4* elem4, Elem9* elem9, Elem16* elem16, Grid* grid, GlobalData* globaldata) {
	derivativesElem(elem4, elem9, elem16, &derivativesScheme, numOfPoints);
	numOfPoints *= numOfPoints;

	double** globalH = new double* [grid->nN]; double** globalHW = new double* [grid->nN];
	double* globalP = new double[grid->nN]; double* globalT = new double[grid->nN]; double* globalPP = new double[grid->nN];
	double** globalC = new double* [grid->nN];

	for (int i = 0; i < grid->nN; i++) {
		globalH[i] = new double[grid->nN];
		globalC[i] = new double[grid->nN];
		globalHW[i] = new double[grid->nN+1];
		globalP[i] = 0; globalPP[i] = 0;
		globalT[i] = globaldata->it;
		for (int j = 0; j < grid->nN; j++) {
			globalH[i][j] = 0;
			globalC[i][j] = 0;
			globalHW[i][j+1] = 0;
		}
	}
	for (int s = 0; s < grid->nE; s++) {
		double x[4]{}, y[4]{};
		int idN[4]{}; bool bcN[4]{};
		//reading file into variables
		for (int j = 0; j < 4; j++) {
			idN[j] = grid->elements[s].ID[j];
			x[j] = grid->nodes[static_cast<std::vector<Node, std::allocator<Node>>::size_type>(idN[j]) - 1].x;
			y[j] = grid->nodes[static_cast<std::vector<Node, std::allocator<Node>>::size_type>(idN[j]) - 1].y;
			bcN[j] = grid->nodes[static_cast<std::vector<Node, std::allocator<Node>>::size_type>(idN[j]) - 1].bc;
		}

		double** jacobian = new double* [numOfPoints] {};
		double** tabX = new double* [numOfPoints]; double** tabY = new double* [numOfPoints];
		double** matrixH = new double* [numOfPoints]; double** matrixHBC = new double* [4]; double** matrixC = new double* [numOfPoints];
		double* detJacobian{ new double [numOfPoints] {} }; double* detJacobianMinus{ new double [numOfPoints] {} };
		double** side1 = new double* [4]; double** side2 = new double* [4]; double** side3 = new double* [4]; double** side4 = new double* [4];
		double* vectorP1{ new double [4] {} }; double* vectorP2{ new double [4] {} }; double* vectorP3{ new double [4] {} }; double* vectorP4{ new double [4] {} }; double* vectorP{ new double [4] {} };

		for (int i = 0; i < 4; i++) {
			side1[i] = new double[4];
			side2[i] = new double[4];
			side3[i] = new double[4];
			side4[i] = new double[4];
			matrixHBC[i] = new double[4];
			for (int j = 0; j < 4; j++) {
				matrixHBC[i][j] = 0.0;
				side1[i][j] = 0.0;
				side2[i][j] = 0.0;
				side3[i][j] = 0.0;
				side4[i][j] = 0.0;
			}
		}
		for (int i = 0; i < numOfPoints; i++) {
			jacobian[i] = new double[4] {};
			tabX[i] = new double[4];
			tabY[i] = new double[4];
			matrixH[i] = new double[4];
			matrixC[i] = new double[4];
			for (int j = 0; j < 4; j++) {
				jacobian[i][j] = 0;
				matrixH[i][j] = 0;
				matrixC[i][j] = 0;
				detJacobian[i] = 0;
				detJacobianMinus[i] = 0;
				vectorP1[j] = 0;
				vectorP2[j] = 0;
				vectorP3[j] = 0;
				vectorP4[j] = 0;
				vectorP[j] = 0;
			}
		}

		//jacobian
		if (numOfPoints == 4) {
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					jacobian[i][0] += elem4->tabEta[i][j] * y[j];
					jacobian[i][1] += elem4->tabKsi[i][j] * y[j] * (-1);
					jacobian[i][2] += elem4->tabEta[i][j] * x[j] * (-1);
					jacobian[i][3] += elem4->tabKsi[i][j] * x[j];
				}
			}
		}
		else if (numOfPoints == 9) {
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					jacobian[i][0] += elem9->tabEta[i][j] * y[j];
					jacobian[i][1] += elem9->tabKsi[i][j] * y[j] * (-1);
					jacobian[i][2] += elem9->tabEta[i][j] * x[j] * (-1);
					jacobian[i][3] += elem9->tabKsi[i][j] * x[j];
				}
			}
		}
		else if (numOfPoints == 16) {
			for (int i = 0; i < numOfPoints; i++) {
				for (int j = 0; j < 4; j++) {
					jacobian[i][0] += elem16->tabEta[i][j] * y[j];
					jacobian[i][1] += elem16->tabKsi[i][j] * y[j] * (-1);
					jacobian[i][2] += elem16->tabEta[i][j] * x[j] * (-1);
					jacobian[i][3] += elem16->tabKsi[i][j] * x[j];
				}
			}
		}
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
						matrixC[j][k] += globaldata->cp * globaldata->d * elem4->nKsiEta[i][k] * elem4->nKsiEta[i][j] * detJacobianMinus[i];
					}
					else if (numOfPoints == 9) {
						matrixH[j][k] += (globaldata->lambda * (tabX[i][k] * tabX[i][j] + tabY[i][k] * tabY[i][j]) * detJacobianMinus[i])
							* w1for3[i] * w2for3[i];
						matrixC[j][k] += (globaldata->cp * globaldata->d * elem9->nKsiEta[i][k] * elem9->nKsiEta[i][j] * detJacobianMinus[i]) * w1for3[i] * w2for3[i];
					}
					else if (numOfPoints == 16) {
						matrixH[j][k] += (globaldata->lambda * (tabX[i][k] * tabX[i][j] + tabY[i][k] * tabY[i][j]) * detJacobianMinus[i])
							* w1for4[i] * w2for4[i];
						matrixC[j][k] += (globaldata->cp * globaldata->d * elem16->nKsiEta[i][k] * elem16->nKsiEta[i][j] * detJacobianMinus[i]) * w1for4[i] * w2for4[i];
					}
				}
			}
		}

		//HBC
		if (grid->nodes[(static_cast<std::vector<Node, std::allocator<Node>>::size_type>(grid->elements[s].ID[0])) - 1].bc == 1 && grid->nodes[(static_cast<std::vector<Node, std::allocator<Node>>::size_type>(grid->elements[s].ID[1])) - 1].bc == 1) {
			side1 = Borders(1, numOfPoints, x, y, globaldata->alpha, globaldata->tot).borders();
			vectorP1 = Borders(1, numOfPoints, x, y, globaldata->alpha, globaldata->tot).vectorP();
		}
		if (grid->nodes[static_cast<std::vector<Node, std::allocator<Node>>::size_type>(grid->elements[s].ID[1]) - 1].bc == 1 && grid->nodes[static_cast<std::vector<Node, std::allocator<Node>>::size_type>(grid->elements[s].ID[2]) - 1].bc == 1) {
			side2 = Borders(2, numOfPoints, x, y, globaldata->alpha, globaldata->tot).borders();
			vectorP2 = Borders(2, numOfPoints, x, y, globaldata->alpha, globaldata->tot).vectorP();
		}
		if (grid->nodes[static_cast<std::vector<Node, std::allocator<Node>>::size_type>(grid->elements[s].ID[2]) - 1].bc == 1 && grid->nodes[static_cast<std::vector<Node, std::allocator<Node>>::size_type>(grid->elements[s].ID[3]) - 1].bc == 1) {
			side3 = Borders(3, numOfPoints, x, y, globaldata->alpha, globaldata->tot).borders();
			vectorP3 = Borders(3, numOfPoints, x, y, globaldata->alpha, globaldata->tot).vectorP();
		}
		if (grid->nodes[static_cast<std::vector<Node, std::allocator<Node>>::size_type>(grid->elements[s].ID[3]) - 1].bc == 1 && grid->nodes[static_cast<std::vector<Node, std::allocator<Node>>::size_type>(grid->elements[s].ID[0]) - 1].bc == 1) {
			side4 = Borders(4, numOfPoints, x, y, globaldata->alpha, globaldata->tot).borders();
			vectorP4 = Borders(4, numOfPoints, x, y, globaldata->alpha, globaldata->tot).vectorP();
		}

		for (int i = 0; i < 4; i++) {
			vectorP[i] = (vectorP1[i] + vectorP2[i] + vectorP3[i] + vectorP4[i]);
			for (int j = 0; j < 4; j++) {
				matrixHBC[i][j] = (side1[i][j] + side2[i][j] + side3[i][j] + side4[i][j]);
				matrixH[i][j] += matrixHBC[i][j];
			}
		}

		//aggregation
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				globalH[grid->elements[s].ID[i] - 1][grid->elements[s].ID[j] - 1] += matrixH[i][j];
				globalC[grid->elements[s].ID[i] - 1][grid->elements[s].ID[j] - 1] += matrixC[i][j];
			}
			globalP[grid->elements[s].ID[i] - 1] += vectorP[i];
		}
	}

	//simultaneous equations
	int size = grid->nN + 1;
	for (int g = 0; g < grid->nN; g++) {
		for (int p = 0; p < grid->nN; p++) {
			globalH[g][p] += (globalC[g][p] / globaldata->st);
		}
	}

	for (int t = globaldata->st; t <= globaldata->t; t += globaldata->st) {
		for (int g = 0; g < grid->nN; g++) {
			globalPP[g] = globalP[g];
			for (int p = 0; p < grid->nN; p++) globalPP[g] += globalC[g][p] / globaldata->st * globalT[p];
		}

		int size = grid->nN + 1;

		for (int i = 0; i < grid->nN; i++){
			for (int j = 0; j < grid->nN; j++) {
				globalHW[i][j] = globalH[i][j];
			}
			globalHW[i][size - 1] = globalPP[i];
		}

		gauss(grid->nN, globalHW, globalT);

		double min = globalT[0];
		double max = globalT[0];
		for (int i = 0; i < grid->nN; i++){
			if (globalT[i] > max) max = globalT[i];
			if (globalT[i] < min) min = globalT[i];
		}

		cout << "Time: " << t << "\t";
		cout << "Tmin: " << min << "\t";
		cout << "Tmax: " << max << std::endl;
	}
}

int main() {
	Grid grid;
	GlobalData globaldata{};
	Integration scheme;
	Elem4 elem4; Elem9 elem9; Elem16 elem16;
	int numOfPoints = 0;

	//reading the file
	readFile(&globaldata, &grid);

	cout << "Enter the number of points: ";
	cin >> numOfPoints;

	//matrixH with the data given
	if (numOfPoints == 2 || numOfPoints == 3 || numOfPoints == 4) matrixH(numOfPoints, &elem4, &elem9, &elem16, &grid, &globaldata);
	else cout << "Wrong number of points : (";

	return 0;
}