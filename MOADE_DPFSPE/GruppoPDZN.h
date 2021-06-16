#pragma once
#include "Permutazione.h"
#include "Modulo.h"
#include <iostream>

using namespace std;

class GruppoPDZN {

	public:
		GruppoPDZN(Permutazione*, Modulo*, Modulo*);
		~GruppoPDZN();

		GruppoPDZN& operator=(GruppoPDZN&);

		void somma(GruppoPDZN*);
		void differenza(GruppoPDZN*);
		virtual void prodotto(double);
		void random();
		void inversa();
		void identita();
		void stampa();
		void scambia(GruppoPDZN*);

		Permutazione* permutazione;
		Modulo* modulo1;
		Modulo* modulo2;
};