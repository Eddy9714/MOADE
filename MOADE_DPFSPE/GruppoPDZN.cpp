#include "GruppoPDZN.h"

GruppoPDZN::GruppoPDZN(Permutazione* p, Modulo* m1, Modulo* m2) {
	permutazione = p;
	modulo1 = m1;
	modulo2 = m2;
}

GruppoPDZN::~GruppoPDZN() {
	delete permutazione;
	delete modulo1;
	delete modulo2;
}

GruppoPDZN& GruppoPDZN::operator=(GruppoPDZN& g) {
	*permutazione = *(g.permutazione);
	*modulo1 = *(g.modulo1);
	*modulo2 = *(g.modulo2);
	return *this;
}

void GruppoPDZN::somma(GruppoPDZN* g) {
	permutazione->somma(g->permutazione);
	modulo1->somma(g->modulo1);
	modulo2->somma(g->modulo2);
}

void GruppoPDZN::differenza(GruppoPDZN* g) {
	permutazione->differenza(g->permutazione);
	modulo1->differenza(g->modulo1);
	modulo2->differenza(g->modulo2);
}

void GruppoPDZN::prodotto(double F) {

	if (F <= 1 && F > 0) {
		unsigned int dM1 = modulo1->distanzaID();
		unsigned int dM2 = modulo2->distanzaID();
		unsigned int dP = permutazione->distanzaID();
		unsigned int dM1Copia = dM1, dM2Copia = dM2, dPCopia = dP;

		unsigned int tot = dM1 + dM2 + dP;

		unsigned int numerodiOperazioni = ceil(tot * F);

		if (dM1Copia == 0 && dM2Copia == 0)
			permutazione->prodotto(F);
		else if (dM2Copia == 0 && dPCopia == 0)
			modulo1->prodotto(F);
		else if(dM1Copia == 0 && dPCopia == 0)
			modulo2->prodotto(F);
		else {

			unsigned int dM1F = 0, dM2F = 0, dPF = 0;

			double pm1, pm2, pp;
			double ran;

			for (unsigned int i = 0; i < numerodiOperazioni; i++, tot--) {
				pm1 = dM1 / (double)tot;
				pm2 = pm1 + (dM2 / (double) tot);
				pp = 1 - pm2;

				ran = genRand.randDouble(0, 1);
				if (ran < pm1)
					dM1F++, dM1--;
				else if (ran < pm2)
					dM2F++, dM2--;
				else dPF++, dP--;
			}

			modulo1->prodotto(dM1F / (double)dM1Copia);
			modulo2->prodotto(dM2F / (double)dM2Copia);
			permutazione->prodotto(dPF / (double)dPCopia);
		}
		
	}
}

void GruppoPDZN::identita() {
	modulo1->identita();
	modulo2->identita();
	permutazione->identita();
}

void GruppoPDZN::random() {
	modulo1->random();
	modulo2->random();
	permutazione->random();
}

void GruppoPDZN::inversa() {
	modulo1->inversa();
	modulo2->inversa();
	permutazione->inversa();
}

void GruppoPDZN::scambia(GruppoPDZN* g) {
	permutazione->scambia(g->permutazione);
	modulo1->scambia(g->modulo1);
	modulo2->scambia(g->modulo2);
}

void GruppoPDZN::stampa() {
	for (unsigned short i = 0; i < permutazione->dimensione; i++) {
		cout << permutazione->individuo[i] << " ";
	}
	cout << ": ";

	for (unsigned short i = 0; i < modulo1->dimensione; i++) {
		cout << modulo1->individuo[i] << " ";
	}

	cout << ": ";
	for (unsigned short i = 0; i < modulo2->dimensione; i++) {
		cout << modulo2->individuo[i] << " ";
	}

	cout << endl;
}