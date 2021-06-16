#pragma once
#include <iostream>
#include <cstring>
#include "Globali.h"
#include <algorithm>

using namespace std;

class Modulo {

	public:
		Modulo(unsigned short, unsigned short);
		Modulo(unsigned short*, unsigned short, unsigned short);
		Modulo(const Modulo& p);
		virtual ~Modulo();
		virtual Modulo& operator=(Modulo&);

		void somma(Modulo*);
		void differenza(Modulo*);
		void prodotto(double);
		void random();
		void randomOrdinato();
		void ordina();
		void inversa();
		void identita();
		void massimo();

		void scambia(Modulo*);

		void stampa();
		void randomZn(Modulo*, unsigned int, unsigned short* = nullptr);
		unsigned int distanzaID();
		unsigned int distanzaRevID();


		unsigned short* individuo;
		unsigned short dimensione;
		unsigned short modulo;
};