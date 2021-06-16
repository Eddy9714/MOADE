
#include "MOADE.h"
#include "PermutazioneT.h"
#include "GruppoPDZN.h"
#include "Istanza.h"
#include "Globali.h"
#include <string>
#include <list>
#include <cfloat>

using namespace std;

struct Individuo {
	GruppoPDZN* rappresentazione;
	Coppia<double> punteggio;
	double alpha;
	vector<unsigned short> simili;
};

class MOADE_DPFSPE : MOADE<Individuo>  {

	public:

		const Istanza istanza;

		MOADE_DPFSPE(string percorso);

		vector<Individuo> esegui(unsigned short, unsigned short, unsigned int, double, double, double, unsigned int s);
		void eliminaPopolazione(vector<Individuo>&);
		void stampa(vector<Individuo>&);

	private:

		vector<double> c;
		vector<vector<Coppia<double>>> o;

		struct InfoInserzione {
			double makeSpan;
			unsigned short posizione;
		};

		Coppia<double> valutaIndividuoParziale(GruppoPDZN*, unsigned short, unsigned short, int = -1);
		void valutaIndividuo(Individuo&, unsigned int&, bool = true);

		void creaPopolazione(vector<Individuo>&, unsigned short);
		void inizializzaPopolazione(vector<Individuo>&, unsigned int&, unsigned short, unsigned short, double, double, Coppia<double>&);
		void combina(vector<Individuo>&, unsigned short, Individuo&);
		void aggiorna(vector<Individuo>&, unsigned short, Individuo&, Coppia<double>&, unsigned int&);

		void ENEH(Individuo&, unsigned int&);
		InfoInserzione miglioreInserzione(GruppoPDZN*, unsigned short, unsigned short, unsigned short);

		void ricercaLocale(Individuo&, unsigned int&);
		void IFLSI(Individuo&, vector<Coppia<unsigned short>>&, unsigned int&);
		void IFLSS(Individuo&, vector<Coppia<unsigned short>>&, unsigned int&);
		void EWFLSI(Individuo&, vector<Coppia<unsigned short>>&, unsigned int&);
		void EWFLSS(Individuo&, vector<Coppia<unsigned short>>&, unsigned int&);

		void ottimizza(vector<Individuo>&, unsigned int&);
		void ottimizzaEnergia(Individuo&, unsigned int&);
		void ottimizzaEnergiaParziale(Individuo&, Coppia<unsigned short>&);

		vector<Coppia<unsigned short>> calcolaInfoFabbriche(GruppoPDZN*);
		CoppiaM<unsigned short, double> trovaFabbricaPeggiore(GruppoPDZN*, vector<Coppia<unsigned short>>&);
};