#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include "Globali.h"

using namespace std;

class Istanza {
	public:

		unsigned short macchine;
		unsigned short fabbriche;
		unsigned short lavori;
		unsigned short nVelocita = 5;

		unsigned short** tp; //tempi lavorazione
		unsigned short** s; //tempi di setup

		double* pi; //potenza idle
		double* ps; //potenza setup

		double* consumoPerVelocita;
		double* velocita;

		Istanza(string path);
		~Istanza();
};