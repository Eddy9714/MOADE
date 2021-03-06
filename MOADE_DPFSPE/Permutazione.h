#pragma once
#include <cstring>
#include "Globali.h"
#include <iostream>

class Permutazione {

	public:
		Permutazione(unsigned short);
		Permutazione(unsigned short*, unsigned short);
		Permutazione(const Permutazione& p);
		virtual ~Permutazione();
		virtual Permutazione& operator=(Permutazione&);

		void somma(Permutazione*);
		void differenza(Permutazione*);
		virtual void prodotto(double) = 0;
		void random();
		void inversa();
		void identita();
		void stampa();
		void scambia(Permutazione*);
		virtual unsigned int distanzaID() = 0;
		virtual unsigned int distanzaRevID() = 0;

		unsigned short* individuo;
		const unsigned short dimensione;
};
