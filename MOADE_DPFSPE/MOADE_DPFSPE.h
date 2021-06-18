
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

		struct InfoFabbrica {
			Coppia<unsigned short> coordinate;
			Coppia<double> valutazione;
		};

		vector<double> c;
		vector<vector<Coppia<double>>> o;
		
		Coppia<double> valutaIndividuoParziale(GruppoPDZN*, unsigned short, unsigned short, int = -1);
		void valutaIndividuo(Individuo&, unsigned int&, bool = true);

		void creaPopolazione(vector<Individuo>&, unsigned short);
		void inizializzaPopolazione(vector<Individuo>&, unsigned int&, unsigned short, unsigned short, double, double, Coppia<double>&);
		void combina(vector<Individuo>&, unsigned short, Individuo&);
		void aggiorna(vector<Individuo>&, unsigned short, Individuo&, Coppia<double>&, unsigned int&);

		void ENEH(Individuo&, unsigned int&);
		CoppiaM<double, unsigned short> miglioreInserzione(GruppoPDZN*, unsigned short, unsigned short, unsigned short);

		void ricercaLocale(vector<Individuo>&, unsigned int&);
		void IFLSI(Individuo&, CoppiaM<vector<InfoFabbrica>, unsigned short>&, unsigned int&, bool = false);
		void IFLSS(Individuo&, CoppiaM<vector<InfoFabbrica>, unsigned short>&, unsigned int&, bool = false);
		void EWFLSI(Individuo&, CoppiaM<vector<InfoFabbrica>, unsigned short>&, unsigned int&, bool = false);
		void EWFLSS(Individuo&, CoppiaM<vector<InfoFabbrica>, unsigned short>&, unsigned int&, bool = false);

		void ottimizza(vector<Individuo>&, unsigned int&);
		void ottimizzaEnergia(Individuo&, unsigned int&);
		void ottimizzaEnergiaParziale(Individuo&, Coppia<unsigned short>&);
		vector<unsigned short> paretoFront(vector<Individuo>&);

		vector<Coppia<unsigned short>> calcolaPosFabbriche(GruppoPDZN*);
		CoppiaM<vector<InfoFabbrica>, unsigned short> calcolaInfoFabbriche(GruppoPDZN*);
};