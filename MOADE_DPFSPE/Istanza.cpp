#pragma once
#include "Istanza.h"

using namespace std;

Istanza::Istanza(string percorso) {

	ifstream file(percorso);

	if (file.is_open()) {
		file >> lavori >> fabbriche >> macchine;

		tp = new unsigned short* [lavori];
		s = new unsigned short* [lavori];

		pi = new double[macchine];
		ps = new double[macchine];

		velocita = new double[5];
		consumoPerVelocita = new double[5];

		for (unsigned short i = 0; i < 5; i++) {
			file >> velocita[i];
		}

		for (unsigned short i = 0; i < 5; i++) {
			file >> consumoPerVelocita[i];
			consumoPerVelocita[i] /= 60;
		}

		for (unsigned short i = 0; i < lavori; i++) {
			tp[i] = new unsigned short[macchine];
			for (unsigned short m = 0; m < macchine; m++) {
				file >> tp[i][m];
			}
		}

		for (unsigned short i = 0; i < lavori; i++) {
			s[i] = new unsigned short[macchine];
			for (unsigned short m = 0; m < macchine; m++) {
				file >> s[i][m];
			}
		}

		for (unsigned short m = 0; m < macchine; m++) {
			file >> pi[m];
			pi[m] /= 60;
		}


		for (unsigned short m = 0; m < macchine; m++) {
			file >> ps[m];
			ps[m] /= 60;
		}
	}
	else {
		cerr << "File non trovato!" << endl << endl;
		exit(-2);
	}
}

Istanza::~Istanza(){

	for (int i = 0; i < lavori; i++) {
		delete[] tp[i];
		delete[] s[i];
	}

	delete[] tp;
	delete[] s;

	delete[] consumoPerVelocita;
	delete[] velocita;
}