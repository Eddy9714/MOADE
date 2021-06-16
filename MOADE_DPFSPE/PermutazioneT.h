#include "Permutazione.h"
#include <cstring>
#include "Random.h"
#include "Globali.h"

class PermutazioneT : public Permutazione {

public:
	PermutazioneT(unsigned short);
	PermutazioneT(unsigned short*, unsigned short);
	PermutazioneT(const PermutazioneT& p);

	void prodotto(double);

	unsigned int distanzaID();
	unsigned int distanzaRevID();

private:
	unsigned int excWeight(PermutazioneT*);
	unsigned int nCycles(PermutazioneT*);
	void randomSS(PermutazioneT*, unsigned short&, double, Coppia<unsigned int>*);
	void randomMergeSS(PermutazioneT*, unsigned short&, double, Coppia<unsigned int>*);
};