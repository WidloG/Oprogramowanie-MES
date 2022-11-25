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

//coordinates of single node
struct Node {
	float x, y;
	int n;	
	bool bc; 
};

//the id number and it's value
struct Element {
	int ID[5];
};

//grid with numbers of nodes and elements
struct Grid {
	int nN = 1; //number of nodes
	int nE = 1; //number of elements
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
	double weights3[3] = { 0.55555, 0.88888, 0.55555 };
};

//dynamic arrays for derivatives with 2points
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

//dynamic arrays for derivatives with 3points
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

//scheme for caculating derivatives with 2 points
void derivativesScheme4(int i) {
	double pcKsi[4] = { -1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3) };
	double pcEta[4] = { -1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3) };

	nKsi[0] = -0.25 * (1 - pcEta[i]); 
	nKsi[1] = 0.25 * (1 - pcEta[i]);
	nKsi[2] = 0.25 * (1 + pcEta[i]);
	nKsi[3] = -0.25 * (1 + pcEta[i]);

	nEta[0] = -0.25 * (1 - pcKsi[i]);
	nEta[1] = -0.25 * (1 + pcKsi[i]);
	nEta[2] = 0.25 * (1 + pcKsi[i]);
	nEta[3] = 0.25 * (1 - pcKsi[i]);
}

//scheme for caculating derivatives with 3 points
void derivativesScheme9(int i) {
	double pcKsi[9] = { -sqrt(0.6), 0, sqrt(0.6), -sqrt(0.6), 0, sqrt(0.6), -sqrt(0.6), 0 ,sqrt(0.6) };
	double pcEta[9] = { -sqrt(0.6), -sqrt(0.6), -sqrt(0.6), 0, 0, 0, sqrt(0.6), sqrt(0.6) , sqrt(0.6) };

	nKsi[0] = -0.25 * (1 - pcEta[i]);
	nKsi[1] = 0.25 * (1 - pcEta[i]);
	nKsi[2] = 0.25 * (1 + pcEta[i]);
	nKsi[3] = -0.25 * (1 + pcEta[i]);

	nEta[0] = -0.25 * (1 - pcKsi[i]);
	nEta[1] = -0.25 * (1 + pcKsi[i]);
	nEta[2] = 0.25 * (1 + pcKsi[i]);
	nEta[3] = 0.25 * (1 - pcKsi[i]);
}

//reading the file
void readFile(GlobalData* globaldata, Grid* grid, list <Node> *listOfNodes, list <Element>* listOfElements) {

	fstream dataFile;
	dataFile.open("Test1_4_4.txt", ios::in);
	string line;

	if (dataFile.is_open()) {
		string name;
		dataFile >> name >> globaldata->t >> name >> globaldata->st >> name >> globaldata->lambda >> name >> globaldata->alpha >> name >>
			globaldata->tot >> name >> globaldata->it >> name >> globaldata->d >> name >> globaldata->cp >> name >> name >> grid->nN >>
			name >> name >> grid->nE >> name;
		
		//wypisywanie
		/*
		cout << "Simulation Time: " << globaldata->t << endl;
		cout << "Simulation Step Time: " << globaldata->st << endl;
		cout << "Conductivity: " << globaldata->lambda << endl;
		cout << "Alpha: " << globaldata->alpha << endl;
		cout << "Tot: " << globaldata->tot << endl;
		cout << "Initial Temp: " << globaldata->it << endl;
		cout << "Density: " << globaldata->d << endl;
		cout << "Specific Heat: " << globaldata->cp << endl;
		cout << "Nodes Numer: " << grid->nN << endl;
		cout << "Elements Number: " << grid->nE << endl;
		cout << "\nNodes: " << endl;
		*/
		getline(dataFile, line);
		Node node;
		for (int i = 0; i < grid->nN; i++) {
			getline(dataFile, line);
			vector<string> useless;
			stringstream ss(line);

			while (ss.good()) {
				string useLess;
				getline(ss, useLess, ',');
				useless.push_back(useLess);
			}

			//konwersja ze stringa do liczby
			node.n = atoi(useless[0].c_str());
			node.x = atof(useless[1].c_str());
			node.y = atof(useless[2].c_str());

			//dodawanie gotowych nodeów do listy
			listOfNodes->insert(listOfNodes->end(), node);
		}

		/*
		list<Node>::iterator it;
		for (it = listOfNodes->begin(); it != listOfNodes->end(); ++it) printf("%d  %.10f  %.10f\n", it->n, it->x, it->y);
		*/
		//cout << "\nElements: " << endl;
		getline(dataFile, line);
		for (int i = 0; i < grid->nE; i++) {
			getline(dataFile, line);
			vector<string> useless;
			stringstream ss(line);

			while (ss.good()) {
				string useLess;
				getline(ss, useLess, ',');
				useless.push_back(useLess);
			}

			Element element;
			element.ID[0] = atoi(useless[0].c_str());
			element.ID[1] = atoi(useless[1].c_str());
			element.ID[2] = atoi(useless[2].c_str());
			element.ID[3] = atoi(useless[3].c_str());
			element.ID[4] = atoi(useless[4].c_str());

			listOfElements->insert(listOfElements->end(), element);

		}

		/*
		list<Element>::iterator it2;
		for (it2 = listOfElements->begin(); it2 != listOfElements->end(); ++it2) cout << it2->ID[0] << " " << it2->ID[1] << " " << it2->ID[2] << " " << it2->ID[3] << " " << it2->ID[4] << endl;
		*/
		getline(dataFile, line);
		vector<string> useless;
		stringstream ss(line);

		while (ss.good()) {
			string useLess;
			getline(ss, useLess, ',');
			useless.push_back(useLess);
		}
		// teraz w useless jest zaczytane liczby BC
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

//counting derivatives in two arrays for Elem4
void derivativesElem(Elem4* elem4, Elem9* elem9, void(*derivativesScheme4)(int), void(*derivativesScheme9)(int), int numOfPoints){

	if (numOfPoints == 2) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				derivativesScheme4(i);
				elem4->tabKsi[i][j] = nKsi[j];
			}
		}

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				derivativesScheme4(i);
				elem4->tabEta[i][j] = nEta[j];
			}
		}

		//elem4->printElem();
	}

	if (numOfPoints == 3) {
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++) {
				derivativesScheme9(i);
				elem9->tabKsi[i][j] = nKsi[j];
			}
		}

		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++) {
				derivativesScheme9(i);
				elem9->tabEta[i][j] = nEta[j];
			}
		}
		//elem9->printElem();
	}
}

void matrixH(int numOfPoints, list <Node>* listOfNodes, Elem4* elem4, Elem9 *elem9) {
	derivativesElem(elem4, elem9, &derivativesScheme4, &derivativesScheme9, numOfPoints);
	numOfPoints *= numOfPoints;

	double w1For3Pc[9] = { 0.555555556, 0.888888889, 0.555555556, 0.555555556, 0.888888889, 0.555555556, 0.555555556, 0.888888889, 0.555555556 };
	double w2For3Pc[9] = { 0.555555556, 0.55555555656, 0.555555556, 0.888888889, 0.888888889, 0.888888889, 0.555555556, 0.555555556, 0.555555556 };
	double w1For4Pc[16] = { 0.347855, 0.652145, 0.652145, 0.347855, 0.347855, 0.652145, 0.652145, 0.347855, 0.347855, 0.652145, 0.652145, 0.347855, 0.347855, 0.652145, 0.652145, 0.347855 };
	double w2For4Pc[16] = { 0.347855, 0.347855, 0.347855, 0.347855, 0.652145, 0.652145, 0.652145, 0.652145, 0.652145, 0.652145, 0.652145, 0.652145, 0.347855, 0.347855, 0.347855, 0.347855 };

	double x[4] = {0, 0.025, 0.025,0};
	double y[4] = {0, 0, 0.025, 0.025};
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

	for (int i = 0; i < numOfPoints; i++){
		for (int j = 0; j < 4; j++) {
			jacobian[i][j] = 0;
			matrixH[i][j] = 0;
			detJacobian[i] = 0;
			detJacobianMinus[i] = 0;
		}
	}

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
	else {
		for (int i = 0; i < numOfPoints; i++) {
			for (int j = 0; j < 4; j++) {
				jacobian[i][0] += elem9->tabKsi[i][j] * x[j];
				jacobian[i][1] += elem9->tabKsi[i][j] * y[j];
				jacobian[i][2] += elem9->tabEta[i][j] * x[j];
				jacobian[i][3] += elem9->tabEta[i][j] * y[j];
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
	else {
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
 
	for (int k = 0; k < numOfPoints; k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (numOfPoints == 4) {
					matrixH[i][j] += (30 * (tabX[k][j] * tabX[k][i] + tabY[k][j] * tabY[k][i]) * detJacobianMinus[k]);
				}
				else if (numOfPoints == 9) {
					matrixH[i][j] += (30 * (tabX[k][j] * tabX[k][i] + tabY[k][j] * tabY[k][i]) * detJacobianMinus[k])
						* w1For3Pc[k] * w2For3Pc[k];
				}

			}
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << matrixH[i][j] << " ";
		}
		cout << endl;
	}
}

int main() {

	//structures and lists go brrrrrr
	Grid grid;
	GlobalData globaldata;
	Integration scheme;
	list<Node> listOfNodes;
	list<Element> listOfElements;
	Elem4 elem4;
	Elem9 elem9;

	//reading the file
	readFile(&globaldata, &grid, &listOfNodes, &listOfElements);

	//results of integration
	//solution integration(&scheme, numOfPoints, numOfDimention)
	/*
	cout << "Integration without bounds for f(x) = 2x^2 + 3x - 8: \n" << endl;
	cout << "Integration for 1D and 2 points: " << integration(&scheme, 2, 1) << endl;
	cout << "Integration for 1D and 3 points: " << integration(&scheme, 3, 1) << endl << endl;;
	cout << "Integration without bounds for f(x) = -5ksi^2 * eta + 2 * ksi * eta^2 + 10: \n" << endl;
	cout << "Integration for 2D and 2 points: " << integration(&scheme, 2, 2) << endl;
	cout << "Integration for 2D and 3 points: " << integration(&scheme, 3, 2) << endl << endl;
	cout << "Integration with bounds [-3,6.5] for f(x) = 3x^2 - 6x + 1: \n" << endl;
	cout << "Integration for 1D and 2 points:" << integrationBounds(&scheme, 2, -3, 6.5) << endl;
	cout << "Integration for 1D and 3 points:" << integrationBounds(&scheme, 3, -3, 6.5) << endl << endl;
	*/

	//results of derivatives
	//solution derivativesElem(&elem4, &elem9, &derivativesScheme4, &derivativesScheme9, numOfPoints)
	//derivativesElem(&elem4, &elem9, &derivativesScheme4, &derivativesScheme9, 2);
	//derivativesElem(&elem4, &elem9, &derivativesScheme4, &derivativesScheme9, 3);

	//results of MatrixH
	//solution matrixH(numofPoints, &listOfNodes, &elem4, &elem9);
	//matrixH(2, &listOfNodes, &elem4, &elem9);
	//matrixH(3, &listOfNodes, &elem4, &elem9);

	return 0;

}