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

		double** tp; //tempi lavorazione
		double** s; //tempi di setup

		unsigned short* pi; //potenza idle
		unsigned short* ps; //potenza setup

		double* consumoPerVelocita;
		double* velocita;

		Istanza(string path);
		~Istanza();
};