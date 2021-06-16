#include "Modulo.h"

Modulo::Modulo(const Modulo& m) : dimensione(m.dimensione), modulo(m.modulo) {
	individuo = new unsigned short[m.dimensione];
	memcpy(individuo, m.individuo, dimensione * sizeof(unsigned short));
}

Modulo::Modulo(unsigned short d, unsigned short r) : dimensione(d), modulo(r)  {
	individuo = new unsigned short[d];
}

Modulo::Modulo(unsigned short* m, unsigned short d, unsigned short r) : dimensione(d), modulo(r) {
	individuo = new unsigned short[dimensione];
	memcpy(individuo, m, dimensione * sizeof(unsigned short));
}

Modulo& Modulo::operator=(Modulo& m) {
	memcpy(this->individuo, m.individuo, sizeof(unsigned short) * m.dimensione);
	return *this;
}

Modulo::~Modulo() {
	delete[] individuo;
}


void Modulo::somma(Modulo* m) {
	for (unsigned short i = 0; i < dimensione; i++) {
		individuo[i] = ((unsigned int)individuo[i] + m->individuo[i]) % modulo;
	}
}

void Modulo::differenza(Modulo* m) {

	int accomulatore;

	for (unsigned short i = 0; i < dimensione; i++) {
		accomulatore = (int)individuo[i] - (int)m->individuo[i];
		if (accomulatore < 0)
			accomulatore = modulo + accomulatore;

		individuo[i] = accomulatore;
	}

}

void Modulo::inversa() {
	for (unsigned short i = 0; i < dimensione; i++) {
		individuo[i] = (modulo - individuo[i]) % modulo;
	}
}

void Modulo::identita() {
	memset(individuo, 0, sizeof(unsigned short) * dimensione);
}

void Modulo::massimo() {
	fill(individuo, individuo + dimensione, modulo - 1);
}

void Modulo::prodotto(double F) {
	
	if (F > 0) {

		unsigned short* ordinamento = new unsigned short[dimensione];

		for (unsigned short i = 0; i < dimensione; i++) {
			ordinamento[i] = i;
		}

		unsigned int numeroSomme = distanzaID();

		if (F <= 0.5) {
			Modulo copia(*this);
			copia.inversa();

			unsigned short operazioni = ceil(numeroSomme * F);

			this->identita();
			randomZn(&copia, operazioni, this->individuo);
		}
		else {
			unsigned short operazioni = ceil(numeroSomme * (1-F));
			randomZn(this, operazioni);
		}

		delete[] ordinamento;
	}
}

void Modulo::randomZn(Modulo* m, unsigned int numeroOperazioni, unsigned short* risultato) {

	unsigned short* permutazione = new unsigned short[dimensione];
	unsigned short limite = dimensione, posizione;
	int acc;
	bool ramo;

	unsigned short tmp;
	for (int i = dimensione - 1; i >= 0; i--) {
		if (m->individuo[i] == 0) {
			limite--;
			tmp = permutazione[i];
			permutazione[i] = permutazione[limite];
			permutazione[limite] = tmp;
		}
		else permutazione[i] = i;
	}

	for (unsigned short i = 0; i < numeroOperazioni; i++) {
		posizione = genRand.randIntU(0, limite - 1);

		ramo = (modulo % 2 == 0 && m->individuo[permutazione[posizione]] == modulo / 2 && genRand.randDouble(0, 1) < 0.5)
			|| (modulo % 2 == 1 && m->individuo[permutazione[posizione]] == modulo / 2);

		if (m->individuo[permutazione[posizione]] < modulo / 2 || ramo) {

			m->individuo[permutazione[posizione]]--;

			if (risultato) {
				acc = risultato[permutazione[posizione]] - 1;
				if (acc < 0) acc = modulo - 1;


				risultato[permutazione[posizione]] = acc;
			}
		}
		else {
			m->individuo[permutazione[posizione]] = (m->individuo[permutazione[posizione]] + 1) % modulo;

			if (risultato) {
				risultato[permutazione[posizione]]++;
			}
		}

		if (m->individuo[permutazione[posizione]] == 0) {
			limite--;
			tmp = permutazione[posizione];
			permutazione[posizione] = permutazione[limite];
			permutazione[limite] = tmp;
		}
	}

	delete[] permutazione;
}

unsigned int Modulo::distanzaID() {
	unsigned int contatore = 0;

	for (unsigned short i = 0; i < dimensione; i++) {
		contatore += min(individuo[i], (unsigned short)(modulo - individuo[i]));
	}

	return contatore;
}

unsigned int Modulo::distanzaRevID() {
	Modulo m(this->dimensione, this->modulo);
	fill(m.individuo, m.individuo + m.dimensione, ceil(this->modulo / 2.));
	m.differenza(this);
	return m.distanzaID();
}

void Modulo::random() {
	for (unsigned short i = 0; i < dimensione; i++) {
		individuo[i] = genRand.randIntU(0, modulo - 1);
	}
}

void Modulo::randomOrdinato() {

	unsigned short limite = modulo - 1;
	for (int i = dimensione - 1; i >= 0; i--) {
		individuo[i] = genRand.randIntU(0, limite);
		limite = individuo[i];
	}
}

void Modulo::ordina() {
	sort(individuo, individuo + dimensione);
}

void Modulo::scambia(Modulo* m) {
	unsigned short* tmpIndividuo = m->individuo;
	unsigned short tmpDimensione = m->dimensione;
	unsigned short tmpModulo = m->modulo;

	m->individuo = this->individuo;
	m->dimensione = this->dimensione;
	m->modulo = this->modulo;

	this->dimensione = tmpDimensione;
	this->individuo = tmpIndividuo;
	this->modulo = tmpModulo;
}

void Modulo::stampa() {
	for (unsigned i = 0; i < dimensione; i++) {
		cout << individuo[i] << " ";
	}
	cout << endl;
}
