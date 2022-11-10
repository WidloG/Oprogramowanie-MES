#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>

using namespace std;

struct Node {
	float x; //wspolrzedna x
	float y; //wspolrzedna y
	int n;	//numer wezla
	bool bc; //warunek brzegowy
};

struct Element {
	int ID[5]; //numer + wartoœci
};

struct Grid {
	int nN = 1; //number of nodes
	int nE = 1; //number of elements
};

struct GlobalData {
	int t;		//SimulationTime
	int st;		//SimulationStepTime 
	int lambda;	//Conductivity
	int alpha;	//Alfa - wspolczynnik przewodzenia
	int tot;	//Tot - temperatura otoczenia
	int it;		//InitialTime
	int d;		//Density
	int cp;		//SpecificHeat
};

struct Integration {
	//punkty ca³kowania i ich wagi
	double nodes2[2] = { -1 / sqrt(3), 1 / sqrt(3) };
	double weights2[2] = { 1,1 };

	double nodes3[3] = { -sqrt(0.6), 0 , sqrt(0.6) };
	double weights3[3] = { 0.55555, 0.88888, 0.55555 };
};

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

double function1D(double x) {
	return 2 * x * x + 3 * x - 8;
}

double function2D(double ksi, double eta) {
	return -5 * ksi * ksi * eta + 2 * ksi * eta * eta + 10;
}

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


int main() {

	Grid grid;
	GlobalData globaldata;
	Integration scheme;
	list<Node> listOfNodes;
	list<Element> listOfElements;
	readFile(&globaldata, &grid, &listOfNodes, &listOfElements);

	cout << "Integration for 1D and 2 points: " << integration(&scheme, 2, 1) << endl;
	cout << "Integration for 1D and 3 points: " << integration(&scheme, 3, 1) << endl;
	cout << "Integration for 2D and 2 points: " << integration(&scheme, 2, 2) << endl;
	cout << "Integration for 2D and 3 points: " << integration(&scheme, 3, 2) << endl;


	system("pause");
}