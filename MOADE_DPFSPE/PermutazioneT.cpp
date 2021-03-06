#include "PermutazioneT.h"

using namespace std;

PermutazioneT::PermutazioneT(const PermutazioneT& p) : Permutazione(p) {}

PermutazioneT::PermutazioneT(unsigned short d) : Permutazione(d) {}

PermutazioneT::PermutazioneT(unsigned short* p, unsigned short d) : Permutazione(p, d) {}

void PermutazioneT::prodotto(double F) {

	if (F >= 0) {

		Coppia<unsigned int>* coppie = new Coppia<unsigned int>[dimensione];
		unsigned short spostamenti;

		if (F < 1) {
			PermutazioneT copia(*this);
			copia.inversa();

			randomSS(&copia, spostamenti, F, coppie);

			this->identita();
		}
		else {
			randomMergeSS(this, spostamenti, F, coppie);
		}

		unsigned short tmp;
		for (unsigned short i = 0; i < spostamenti; i++) {
			tmp = this->individuo[coppie[i].x];
			this->individuo[coppie[i].x] = this->individuo[coppie[i].y];
			this->individuo[coppie[i].y] = tmp;
		}

		delete[] coppie;
	}
}

void PermutazioneT::randomSS(PermutazioneT* p, unsigned short& spostamenti, double F, Coppia<unsigned int>* coppie) {

	unsigned short** c = new unsigned short*[p->dimensione];
	unsigned short* lc = new unsigned short [p->dimensione];
	bool* v = new bool[p->dimensione];

	memset(v, 0, sizeof(bool) * p->dimensione);
	memset(lc, 0, sizeof(unsigned short) * p->dimensione);

	unsigned short nc1 = 0, nc = 0;

	int i, j, k, nit = 0, lim, r, nexc = 0, psum, t, lc1, lc2;

	for (i = 0; i < p->dimensione; i++) {
		c[i] = new unsigned short[p->dimensione];
	}

	for (i = 0; i < p->dimensione; i++) {
		if (v[i]) continue;
		j = i;
		while (!v[j]) {
			c[nc][lc[nc]++] = j;
			v[j] = true;
			j = p->individuo[j];
		}
		nc1++; 

		if (lc[nc] > 1) {
			nexc += lc[nc] * (lc[nc] - 1) / 2;
			nc++;
		}
		else
			lc[nc] = 0;
	}

	lim = ceil(F * (p->dimensione - nc1)); 

	spostamenti = 0;
	while (nc > 0 && nit < lim) {

		r = genRand.randIntU(0, nexc-1);
		psum = 0;
		for (k = 0; k < nc; k++) {
			psum += lc[k] * (lc[k] - 1) / 2;
			if (r < psum) break;
		}

		genRand.dueIndiciRandom(lc[k], i, j);

		coppie[spostamenti].x = c[k][i];
		coppie[spostamenti].y = c[k][j];
		spostamenti++;

		nexc -= lc[k] * (lc[k] - 1) / 2; 
		if (i > j) {
			t = i;
			i = j;
			j = t;
		}

		lc1 = i - j + lc[k];	
		lc2 = j - i;			

		if (lc1 > 1 && lc2 > 1) { 

			memcpy(c[nc], c[k] + (i + 1), sizeof(unsigned short) * lc2);
			lc[nc] = lc2;
			nc++;
			
			t = lc[k] - j - 1;	
			if (t > 0) memmove(c[k] + (i + 1), c[k] + (j + 1), sizeof(unsigned short) * t);
			lc[k] = lc1;

			nexc += (lc1 * (lc1 - 1) + lc2 * (lc2 - 1)) / 2;	
		}
		else if (lc1 > 1 && lc2 == 1) { 

			t = lc[k] - j - 1;	
			if (t > 0) memmove(c[k] + (i + 1), c[k] + (j + 1), sizeof(unsigned short) * t);
			lc[k] = lc1;
	
			nexc += lc1 * (lc1 - 1) / 2;
		}
		else if (lc2 > 1) { 

			memmove(c[k], c[k] + (i + 1), sizeof(unsigned short) * lc2);
			lc[k] = lc2;
		
			nexc += lc2 * (lc2 - 1) / 2;	
		}
		else { 

			nc--;
			if (k != nc) {
				memcpy(c[k], c[nc], sizeof(unsigned short) * lc[nc]);
				lc[k] = lc[nc];
			}
			
		}

		nit++;
	}
	delete[] v;
	delete[] lc;

	for (i = 0; i < p->dimensione; i++) {
		delete[] c[i];
	}
	delete[] c;
}

void PermutazioneT::randomMergeSS(PermutazioneT* p, unsigned short& spostamenti, double F, Coppia<unsigned int>* coppie) {

	unsigned short** c = new unsigned short* [p->dimensione];
	unsigned short* lc = new unsigned short[p->dimensione];
	bool* v = new bool[p->dimensione];

	memset(v, 0, sizeof(bool) * p->dimensione);
	memset(lc, 0, sizeof(unsigned short) * p->dimensione);

	int i, j, k, nc = 0, nit = 0, lim, r, p1sum = 0, psum, t;

	for (i = 0; i < p->dimensione; i++) {
		c[i] = new unsigned short[p->dimensione];
	}

	for (i = 0; i < p->dimensione; i++) {
		if (v[i]) continue;
		j = i;
		while (!v[j]) {
			c[nc][lc[nc]++] = j;
			v[j] = true;
			j = p->individuo[j];
		}
		p1sum += lc[nc] * (p->dimensione - lc[nc]);
		nc++;
	}

	k = ceil(F * (p->dimensione - nc)); 
	if (k > p->dimensione - 1) k = p->dimensione - 1; 
	lim = k - p->dimensione + nc; 

	spostamenti = 0;
	while (nc > 1 && nit < lim) {
		
		r = genRand.randIntU(0, p1sum - 1);
		psum = 0;
		for (i = 0; i < nc; i++) {
			psum += lc[i] * (p->dimensione - lc[i]);
			if (r < psum) break;
		}
		
		r = genRand.randIntU(0, p->dimensione - lc[i] - 1);
		psum = 0;
		for (j = 0; j < nc; j++) {
			if (j != i) psum += lc[j];
			if (r < psum) break;
		}
		
		coppie[spostamenti].x = c[i][genRand.randIntU(0, lc[i] - 1)];
		coppie[spostamenti].y = c[j][genRand.randIntU(0, lc[j] - 1)];
		spostamenti++;

		p1sum -= lc[i] * (p->dimensione - lc[i]) + lc[j] * (p->dimensione - lc[j]);
		
		if (lc[i] < lc[j]) {
			t = i;
			i = j;
			j = t;
		}

		memcpy(c[i] + lc[i], c[j], sizeof(unsigned short) * lc[j]);
		lc[i] += lc[j];
		
		nc--;
		if (j != nc) {
			memcpy(c[j], c[nc], sizeof(unsigned short) * lc[nc]);
			lc[j] = lc[nc];
		}

		p1sum += lc[i] * (p->dimensione - lc[i]);
		nit++;
	}

	delete[] v;
	delete[] lc;

	for (i = 0; i < p->dimensione; i++) {
		delete[] c[i];
	}
	delete[] c;
}

unsigned int PermutazioneT::nCycles(PermutazioneT* p) {
	bool* v = new bool[p->dimensione];

	int count = 0, i, j;
	memset(v, 0, sizeof(bool) * p->dimensione);

	for (i = 0; i < p->dimensione; i++) {
		if (v[i]) continue;
		j = i;
		while (!v[j]) {
			v[j] = true;
			j = p->individuo[j];
		}
		count++;
	}

	delete[] v;
	return count;
}


unsigned int PermutazioneT::excWeight(PermutazioneT* p) {
	return p->dimensione - nCycles(p);
}

unsigned int PermutazioneT::distanzaID() {
	return excWeight(this);
}
unsigned int PermutazioneT::distanzaRevID() {
	return this->dimensione - 1 - excWeight(this);
}