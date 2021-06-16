#pragma once
#include "IndiciRandom.h"
#include "Globali.h"
#include <chrono>
#include <iostream>
#include <algorithm>

using namespace std;

template <class T> class MOADE {

	protected:
		unsigned int seed;

		virtual void creaPopolazione(vector<T>&, unsigned short) = 0;
		virtual void inizializzaPopolazione(vector<T>&, unsigned int&, unsigned short, unsigned short, double, double, Coppia<double>&) = 0;
		virtual void eliminaPopolazione(vector<T>&) = 0;
		virtual void combina(vector<T>&, unsigned short, T&) = 0;
		virtual void aggiorna(vector<T>&, unsigned short, T&, Coppia<double>&, unsigned int&) = 0;
		virtual void ricercaLocale(T&, unsigned int&) = 0;
		virtual void ottimizza(vector<T>&, unsigned int&) = 0;
		virtual void valutaIndividuo(T&, unsigned int&, bool = true) = 0;
		virtual void stampa(vector<T>&) = 0;

		vector<T> esegui(unsigned short H, unsigned short Tsize, unsigned int numeroValutazioni, double Fmax, double alphaMin, double alphaMax, unsigned int s) {

			unsigned short nIndividui = H + 1;

			using orologio = std::chrono::system_clock;
			using sec = std::chrono::duration<double>;

			auto tempoIniziale = orologio::now();
			auto tempoAttuale = tempoIniziale;

			unsigned int valutazioniEffettuate = 0;
			unsigned int numeroGenerazioni = 0;

			seed = s;
			if(seed > 0) genRand.impostaSeed(seed);

			vector<T> popolazione;
			popolazione.reserve(nIndividui);

			vector<T> figlio;
			popolazione.reserve(1);

			Coppia<double> migliori = { DBL_MAX, DBL_MAX };

			creaPopolazione(popolazione, nIndividui);
			creaPopolazione(figlio, 1);

			inizializzaPopolazione(popolazione, valutazioniEffettuate, H, Tsize, alphaMin, alphaMax, migliori);

			while (valutazioniEffettuate < numeroValutazioni) {
				numeroGenerazioni++;

				sec intervallo = orologio::now() - tempoAttuale;
				if (intervallo.count() > 1) {
					cout << "G: " << numeroGenerazioni << "\n";
					tempoAttuale = orologio::now();
				}

				for (unsigned short i = 0; i < nIndividui; i++) {
					combina(popolazione, i, figlio[0]);	
					aggiorna(popolazione, i, figlio[0], migliori, valutazioniEffettuate);
				}

				if (valutazioniEffettuate > numeroValutazioni / 1.5) {
					for (unsigned short i = 0; i < popolazione.size(); i++) {
						ricercaLocale(popolazione[i], valutazioniEffettuate);
					}
				}
				
			}

			ottimizza(popolazione, valutazioniEffettuate);
			return popolazione;
		};
};
