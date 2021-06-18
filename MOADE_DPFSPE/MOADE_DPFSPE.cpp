#include "MOADE_DPFSPE.h"

MOADE_DPFSPE::MOADE_DPFSPE(string percorso) : istanza(percorso) {
	c = vector<double>(istanza.macchine);

	o = vector<vector<Coppia<double>>>(istanza.lavori);

	for (unsigned short i = 0; i < o.size(); i++) {
		o[i] = vector<Coppia<double>>(istanza.macchine);
	}
}

vector<Individuo> MOADE_DPFSPE::esegui(unsigned short H, unsigned short T, unsigned int numeroValutazioni, double Fmax, double alphaMin, double alphaMax, unsigned int s) {
	return MOADE::esegui(H, T, numeroValutazioni, Fmax, alphaMin, alphaMax, s);
}

void MOADE_DPFSPE::stampa(vector<Individuo>& popolazione) {
	for (unsigned short i = 0; i < popolazione.size(); i++) {
		popolazione[i].rappresentazione->stampa();
		cout << popolazione[i].punteggio.x << ":" << popolazione[i].punteggio.y << endl << endl;
	}
}

void MOADE_DPFSPE::creaPopolazione(vector<Individuo>& popolazione, unsigned short nIndividui) {
	GruppoPDZN* rappresentazione;

	for (unsigned short i = 0; i < nIndividui; i++) {
		rappresentazione = new GruppoPDZN(new PermutazioneT(istanza.lavori), new Modulo(istanza.fabbriche - 1, istanza.lavori - 1),
			new Modulo(istanza.lavori * istanza.macchine, istanza.nVelocita));

		popolazione.push_back({rappresentazione});
	}
}

void MOADE_DPFSPE::inizializzaPopolazione(vector<Individuo>& popolazione, unsigned int& valutazioniEffettuate, unsigned short H, unsigned short T, 
	double alphaMin, double alphaMax, Coppia<double>& migliori) {

	GruppoPDZN* rappresentazione;

	popolazione[0].rappresentazione->modulo2->identita();
	popolazione[0].alpha = alphaMin;
	popolazione[0].simili.reserve(T);
	ENEH(popolazione[0], valutazioniEffettuate);

	popolazione[popolazione.size() - 1].rappresentazione->modulo2->massimo();
	popolazione[popolazione.size() - 1].alpha = alphaMax;
	popolazione[popolazione.size() - 1].simili.reserve(T);
	ENEH(popolazione[popolazione.size() - 1], valutazioniEffettuate);

	for (unsigned short i = 1; i < popolazione.size() - 1; i++) {
		popolazione[i].rappresentazione->modulo2->random();
		popolazione[i].alpha = alphaMin + ((alphaMax - alphaMin) * (double)i / H);
		popolazione[i].simili.reserve(T);
		ENEH(popolazione[i], valutazioniEffettuate);
	}

	int prima, dopo;

	for (unsigned short i = 0; i < popolazione.size(); i++) {
	
		prima = (int)i - 1;
		dopo = (int)i + 1;

		while (true) {

			if (prima >= 0) {
				popolazione[i].simili.push_back(prima);
				prima--;

				if (popolazione[i].simili.size() == T) break;
			}

			if (dopo < popolazione.size()) {
				popolazione[i].simili.push_back(dopo);
				dopo++;

				if (popolazione[i].simili.size() == T) break;
			}
		}
	}

	migliori.x = popolazione[popolazione.size() - 1].punteggio.x;
	migliori.y = popolazione[0].punteggio.y;
}

void MOADE_DPFSPE::combina(vector<Individuo>& popolazione, unsigned short indice, Individuo& risultato) {
	int i1, i2;

	Individuo& individuo = popolazione[indice];

	genRand.dueIndiciRandom(individuo.simili.size(), i1, i2);

	i1 = individuo.simili[i1];
	i2 = individuo.simili[i2];

	*(risultato.rappresentazione) = *(popolazione[i1].rappresentazione);
	risultato.rappresentazione->differenza(popolazione[i2].rappresentazione);
	risultato.rappresentazione->prodotto(genRand.randDouble(0,1));
	risultato.rappresentazione->somma(individuo.rappresentazione);
	risultato.rappresentazione->modulo1->ordina();

	risultato.alpha = individuo.alpha;
}

void MOADE_DPFSPE::aggiorna(vector<Individuo>& popolazione, unsigned short indice, Individuo& risultato, 
	Coppia<double>& migliori, unsigned int& valutazioniEffettuate) {

	valutaIndividuo(risultato, valutazioniEffettuate);

	if (risultato.punteggio.x < migliori.x)
		migliori.x = risultato.punteggio.x;

	if (risultato.punteggio.y < migliori.y)
		migliori.y = risultato.punteggio.y;
	
	vector<unsigned short>& simili = popolazione[indice].simili;

	Individuo* target;

	double p1, p2;

	for (unsigned short i = 0; i < simili.size(); i++) {
		target = &popolazione[simili[i]];
		p1 = max(target->alpha * target->punteggio.x - migliori.x, (1 - target->alpha) * target->punteggio.y - migliori.y);
		p2 = max(target->alpha * risultato.punteggio.x - migliori.x, (1 - target->alpha) * risultato.punteggio.y - migliori.y);

		if (p2 < p1) {				
			target->punteggio = risultato.punteggio;
			*(target->rappresentazione) = *(risultato.rappresentazione);
		}
	}
}


void MOADE_DPFSPE::eliminaPopolazione(vector<Individuo>& popolazione) {
	for (unsigned short i = 0; i < popolazione.size(); i++) {
		delete popolazione[i].rappresentazione;
	}
}

void MOADE_DPFSPE::valutaIndividuo(Individuo& individuo, unsigned int& valutazioniEffettuate, bool conta) {

	GruppoPDZN* g = individuo.rappresentazione;

	if (conta) valutazioniEffettuate++;

	vector<Coppia<unsigned short>> infoFabbriche = calcolaPosFabbriche(g);
	Coppia<double> coppia, risultato = { 0., 0. };

	for (unsigned short i = 0; i < infoFabbriche.size(); i++) {
		coppia = valutaIndividuoParziale(g, infoFabbriche[i].x, infoFabbriche[i].y);
		risultato.x = max(risultato.x, coppia.x);
		risultato.y += coppia.y;
	}

	individuo.punteggio = risultato;
}

Coppia<double> MOADE_DPFSPE::valutaIndividuoParziale(GruppoPDZN* g, unsigned short inizioFabbrica, unsigned short lunghezza, int posDaEscludere) {

	unsigned short* fabbrica = g->permutazione->individuo + inizioFabbrica;

	unsigned int pos;
	double tempoReale, energiaConsumata = 0.;

	for (unsigned short i = 0; i < lunghezza; i++) {

		if (i != posDaEscludere) {
			for (unsigned short j = 0; j < istanza.macchine; j++) {

				pos = fabbrica[i] * istanza.macchine + j;
				tempoReale = istanza.tp[fabbrica[i]][j] / istanza.velocita[g->modulo2->individuo[pos]];
				energiaConsumata += tempoReale * istanza.consumoPerVelocita[g->modulo2->individuo[pos]];

				tempoReale += istanza.s[fabbrica[i]][j];
				energiaConsumata += istanza.ps[j] * istanza.s[fabbrica[i]][j];

				if (i == 0 && j == 0) {
					c[j] = tempoReale;
				}
				else if (i == 0 && j != 0) {
					c[j] = tempoReale + c[j - 1];
				}
				else if (j != 0) {

					if (c[j - 1] > c[j])
						energiaConsumata += istanza.pi[j] * (c[j - 1] - c[j]);

					c[j] = tempoReale + max(c[j - 1], c[j]);
				}
				else {
					c[0] = tempoReale + c[0];
				}
			}
		}
	}

	double makeSpan = c[istanza.macchine - 1];
	
	return { makeSpan , energiaConsumata};
}


vector<unsigned short> MOADE_DPFSPE::paretoFront(vector<Individuo>& popolazione) {

	vector<unsigned short> indici;

	sort(popolazione.begin(), popolazione.end(), [](Individuo& a, Individuo& b){
		return a.punteggio.x < b.punteggio.x;
	});

	unsigned short i = 0, j=0;

	while (i < popolazione.size()) {
		indici.push_back(i);

		j = i + 1;
		while (j < popolazione.size()) {
			if (popolazione[j].punteggio.y < popolazione[i].punteggio.y) break;
			j++;
		}

		i = j;
	}

	return indici;
}


void MOADE_DPFSPE::ricercaLocale(vector<Individuo>& popolazione, unsigned int& valutazioniEffettuate, unsigned int numeroValutazioni) {

	if (valutazioniEffettuate >= numeroValutazioni / 1.5) {

		vector<unsigned short> indici = paretoFront(popolazione);

		unsigned short indice;
		for (unsigned short i = 0; i < popolazione.size(); i++) {
			indice = indici[genRand.randIntU(0, indici.size() - 1)];
			GruppoPDZN* g = popolazione[indice].rappresentazione;

			auto infoFabbriche = calcolaInfoFabbriche(g);

			unsigned short ricerca = genRand.randIntU(0, 3);

			switch (ricerca) {
			case 0:
				//Inserzione intra-fabbrica peggiore rispetto al makespan
				IFLSI(popolazione[indice], infoFabbriche, valutazioniEffettuate, true);
				break;
			case 1:
				//Swap intra-fabbrica peggiore rispetto al makespan
				IFLSS(popolazione[indice], infoFabbriche, valutazioniEffettuate, true);
				break;
			case 2:
				//Inserzione extra-fabbrica peggiore rispetto al makespan
				EWFLSI(popolazione[indice], infoFabbriche, valutazioniEffettuate, true);
				break;
			case 3:
				//Swap extra-fabbrica peggiore rispetto al makespan
				EWFLSS(popolazione[indice], infoFabbriche, valutazioniEffettuate, true);
				break;
			}
		}
	}
}

void MOADE_DPFSPE::IFLSI(Individuo& individuo, CoppiaM<vector<InfoFabbrica>, unsigned short>& infoFabbriche, 
	unsigned int& valutazioniEffettuate, bool dominazioneRichiesta) {

	GruppoPDZN* g = individuo.rappresentazione;

	InfoFabbrica& infoFabbricaPeggiore = infoFabbriche.x[infoFabbriche.y];
	auto coordinateFabbricaPeggiore = infoFabbricaPeggiore.coordinate;

	unsigned short ran = genRand.randIntU(0, infoFabbricaPeggiore.coordinate.y - 1);
	unsigned short posLavoro = infoFabbricaPeggiore.coordinate.x + ran;
	unsigned short nuovoPosLavoro;
	unsigned short lavoro = g->permutazione->individuo[posLavoro];
	//rimuovi elemento in posizione i e trova la migliore posizione possibile

	if (ran >= coordinateFabbricaPeggiore.y / 2) {
		unsigned short posizioniDaCopiare = coordinateFabbricaPeggiore.y - ran - 1;

		if(posizioniDaCopiare > 0)
			memcpy(g->permutazione->individuo + posLavoro, g->permutazione->individuo + posLavoro + 1,
				sizeof(unsigned short) * posizioniDaCopiare);

		auto info = miglioreInserzione(g, coordinateFabbricaPeggiore.x, coordinateFabbricaPeggiore.y - 1, lavoro);

		nuovoPosLavoro = coordinateFabbricaPeggiore.x + info.y;

		memcpy(g->permutazione->individuo + nuovoPosLavoro + 1, g->permutazione->individuo + nuovoPosLavoro, sizeof(unsigned short) *
			(coordinateFabbricaPeggiore.y - 1 - info.y));
	}
	else {

		if (ran > 0)
			memcpy(g->permutazione->individuo + coordinateFabbricaPeggiore.x + 1, g->permutazione->individuo + coordinateFabbricaPeggiore.x,
				sizeof(unsigned short) * ran);

		auto info = miglioreInserzione(g, coordinateFabbricaPeggiore.x + 1, coordinateFabbricaPeggiore.y - 1, lavoro);

		nuovoPosLavoro = coordinateFabbricaPeggiore.x + info.y;

		memcpy(g->permutazione->individuo + coordinateFabbricaPeggiore.x, g->permutazione->individuo + coordinateFabbricaPeggiore.x + 1, sizeof(unsigned short) *
			(info.y));
	}

	g->permutazione->individuo[nuovoPosLavoro] = lavoro;

	if (nuovoPosLavoro != posLavoro) {
		auto valutazione =
			valutaIndividuoParziale(individuo.rappresentazione, coordinateFabbricaPeggiore.x, coordinateFabbricaPeggiore.y);

		double variazioneEnergia = (valutazione.y - infoFabbricaPeggiore.valutazione.y);

		if (dominazioneRichiesta && variazioneEnergia > 0) {
			//rollback
			//rimuovi da nuovoPosLavoro e inserisci in posLavoro
			if (posLavoro < nuovoPosLavoro) {
				memcpy(g->permutazione->individuo + posLavoro + 1, g->permutazione->individuo + posLavoro,
					sizeof(unsigned short) * (nuovoPosLavoro - posLavoro));
			}
			else {
				memcpy(g->permutazione->individuo + nuovoPosLavoro, g->permutazione->individuo + nuovoPosLavoro + 1,
					sizeof(unsigned short) * (posLavoro - nuovoPosLavoro));
			}

			g->permutazione->individuo[posLavoro] = lavoro;
		}
		else {
			individuo.punteggio.y += variazioneEnergia;
			infoFabbricaPeggiore.valutazione = valutazione;

			double makespanPeggiore = 0, makespan;
			for (unsigned short i = 0; i < infoFabbriche.x.size(); i++) {
				makespan = infoFabbriche.x[i].valutazione.x;
				if (makespan > makespanPeggiore) {
					makespanPeggiore = makespan;
					individuo.punteggio.x = makespanPeggiore;
				}
			}
		}
	}

	valutazioniEffettuate++;
}

void MOADE_DPFSPE::IFLSS(Individuo& individuo, CoppiaM<vector<InfoFabbrica>, unsigned short>& infoFabbriche, 
	unsigned int& valutazioniEffettuate, bool dominanzaRichiesta) {

	GruppoPDZN* g = individuo.rappresentazione;

	InfoFabbrica& infoFabbricaPeggiore = infoFabbriche.x[infoFabbriche.y];
	auto coordinateFabbricaPeggiore = infoFabbricaPeggiore.coordinate;

	if (coordinateFabbricaPeggiore.y < 2) return;

	int i1, i2;
	genRand.dueIndiciRandom(coordinateFabbricaPeggiore.y, i1, i2);

	i1 += coordinateFabbricaPeggiore.x;
	i2 += coordinateFabbricaPeggiore.x;

	unsigned short tmp = g->permutazione->individuo[i1];
	g->permutazione->individuo[i1] = g->permutazione->individuo[i2];
	g->permutazione->individuo[i2] = tmp;

	Coppia<double> valutazione = valutaIndividuoParziale(g, coordinateFabbricaPeggiore.x, coordinateFabbricaPeggiore.y);
	double variazioneEnergia = (valutazione.y - infoFabbricaPeggiore.valutazione.y);

	if (valutazione.x > infoFabbricaPeggiore.valutazione.x || (dominanzaRichiesta && variazioneEnergia > 0)) {
		tmp = g->permutazione->individuo[i1];
		g->permutazione->individuo[i1] = g->permutazione->individuo[i2];
		g->permutazione->individuo[i2] = tmp;		
	}
	else {

		individuo.punteggio.y += variazioneEnergia;
		infoFabbricaPeggiore.valutazione = valutazione;

		double makespanPeggiore = 0, makespan;
		for (unsigned short i = 0; i < infoFabbriche.x.size(); i++) {
			makespan = infoFabbriche.x[i].valutazione.x;
			if (makespan > makespanPeggiore) {
				makespanPeggiore = makespan;
				individuo.punteggio.x = makespanPeggiore;
			}
		}
	}

	valutazioniEffettuate++;
}

void MOADE_DPFSPE::EWFLSI(Individuo& individuo, CoppiaM<vector<InfoFabbrica>, unsigned short>& infoFabbriche, 
	unsigned int& valutazioniEffettuate, bool dominanzaRichiesta) {

	GruppoPDZN* g = individuo.rappresentazione;

	InfoFabbrica& infoFabbricaPeggiore = infoFabbriche.x[infoFabbriche.y];
	auto& coordinateFabbricaPeggiore = infoFabbricaPeggiore.coordinate;

	unsigned short indiceFabbricaScelta = genRand.randIntU(0, istanza.fabbriche - 2);

	if (indiceFabbricaScelta >= infoFabbriche.y)
		indiceFabbricaScelta++;

	InfoFabbrica& infoFabbricaScelta = infoFabbriche.x[indiceFabbricaScelta];
	auto& coordinateFabbricaScelta = infoFabbricaScelta.coordinate;

	unsigned short ran = genRand.randIntU(0, coordinateFabbricaPeggiore.y - 1);
	unsigned short posLavoro = coordinateFabbricaPeggiore.x + ran;
	unsigned short lavoro = g->permutazione->individuo[posLavoro];

	auto info = miglioreInserzione(g, coordinateFabbricaScelta.x, coordinateFabbricaScelta.y, lavoro);

	if (info.x <= infoFabbricaPeggiore.valutazione.x) {

		//Verifichiamo che il makespan della fabbrica peggiore non sia peggiorato
		Coppia<double> val = valutaIndividuoParziale(g, coordinateFabbricaPeggiore.x, coordinateFabbricaPeggiore.y, ran);

		if (val.x <= infoFabbricaPeggiore.valutazione.x) {
			
			//applica modifiche

			unsigned short nuovoPosLavoro = coordinateFabbricaScelta.x + info.y;

			if (posLavoro < nuovoPosLavoro) {
				nuovoPosLavoro--;

				memcpy(g->permutazione->individuo + posLavoro, g->permutazione->individuo + posLavoro + 1,
					sizeof(unsigned short) * (nuovoPosLavoro - posLavoro));

				g->permutazione->individuo[nuovoPosLavoro] = lavoro;

				for (unsigned short k = infoFabbriche.y; k < indiceFabbricaScelta; k++) {
					g->modulo1->individuo[k]--;
				}

				coordinateFabbricaScelta.x--;
			}
			else {
				memcpy(g->permutazione->individuo + nuovoPosLavoro + 1, g->permutazione->individuo + nuovoPosLavoro,
					sizeof(unsigned short) * (posLavoro - nuovoPosLavoro));

				g->permutazione->individuo[nuovoPosLavoro] = lavoro;

				for (unsigned short k = indiceFabbricaScelta; k < infoFabbriche.y; k++) {
					g->modulo1->individuo[k]++;
				}

			}

			coordinateFabbricaPeggiore.y--;
			coordinateFabbricaScelta.y++;

			auto val2 = valutaIndividuoParziale(g, coordinateFabbricaScelta.x, coordinateFabbricaScelta.y);
			
			double variazioneEnergia = (val.y - infoFabbricaPeggiore.valutazione.y) + (val2.y - infoFabbricaScelta.valutazione.y);


			if (dominanzaRichiesta && variazioneEnergia > 0) {
				//rollback
				if (posLavoro < nuovoPosLavoro) {
					memcpy(g->permutazione->individuo + posLavoro + 1, g->permutazione->individuo + posLavoro,
						sizeof(unsigned short) * (nuovoPosLavoro - posLavoro));

					g->permutazione->individuo[posLavoro] = lavoro;

					for (unsigned short k = infoFabbriche.y; k < indiceFabbricaScelta; k++) {
						g->modulo1->individuo[k]++;
					}
				}
				else {

					memcpy(g->permutazione->individuo + nuovoPosLavoro, g->permutazione->individuo + nuovoPosLavoro + 1,
						sizeof(unsigned short) * (posLavoro - nuovoPosLavoro));

					g->permutazione->individuo[posLavoro] = lavoro;

					for (unsigned short k = indiceFabbricaScelta; k < infoFabbriche.y; k++) {
						g->modulo1->individuo[k]--;
					}

					coordinateFabbricaScelta.x++;
				}

				coordinateFabbricaPeggiore.y++;
				coordinateFabbricaScelta.y--;
			}
			else {
				individuo.punteggio.y += variazioneEnergia;

				infoFabbricaPeggiore.valutazione = val;
				infoFabbricaScelta.valutazione = val2;

				double makespanPeggiore = 0, makespan;
				for (unsigned short i = 0; i < infoFabbriche.x.size(); i++) {
					makespan = infoFabbriche.x[i].valutazione.x;
					if (makespan > makespanPeggiore) {
						makespanPeggiore = makespan;
						individuo.punteggio.x = makespanPeggiore;
					}
				}
			}
		}
	}
	valutazioniEffettuate++;
}

void MOADE_DPFSPE::EWFLSS(Individuo& individuo, CoppiaM<vector<InfoFabbrica>, unsigned short>& infoFabbriche, 
	unsigned int& valutazioniEffettuate, bool dominazioneRichiesta) {

	GruppoPDZN* g = individuo.rappresentazione;

	InfoFabbrica& infoFabbricaPeggiore = infoFabbriche.x[infoFabbriche.y];
	auto coordinateFabbricaPeggiore = infoFabbricaPeggiore.coordinate;

	unsigned short indiceFabbricaScelta = genRand.randIntU(0, istanza.fabbriche - 2);

	if (indiceFabbricaScelta >= infoFabbriche.y)
		indiceFabbricaScelta++;

	InfoFabbrica& infoFabbricaScelta = infoFabbriche.x[indiceFabbricaScelta];
	auto coordinateFabbricaScelta = infoFabbricaScelta.coordinate;

	if (coordinateFabbricaScelta.y == 0) return;

	unsigned short posLavoroFabbricaPeggiore = 
		coordinateFabbricaPeggiore.x + genRand.randIntU(0, coordinateFabbricaPeggiore.y - 1);

	unsigned short posLavoroFabbricaScelta = coordinateFabbricaScelta.x +
		genRand.randIntU(0, coordinateFabbricaScelta.y - 1);

	//Scambia lavori
	unsigned short tmp = g->permutazione->individuo[posLavoroFabbricaPeggiore];
	g->permutazione->individuo[posLavoroFabbricaPeggiore] = g->permutazione->individuo[posLavoroFabbricaScelta];
	g->permutazione->individuo[posLavoroFabbricaScelta] = tmp;

	Coppia<double> val1 = valutaIndividuoParziale(g, coordinateFabbricaPeggiore.x, coordinateFabbricaPeggiore.y);
	Coppia<double> val2 = valutaIndividuoParziale(g, coordinateFabbricaScelta.x, coordinateFabbricaScelta.y);

	double variazioneEnergia = (val1.y - infoFabbricaPeggiore.valutazione.y) + (val2.y - infoFabbricaScelta.valutazione.y);

	if (val1.x > infoFabbricaPeggiore.valutazione.x || val2.x > infoFabbricaPeggiore.valutazione.x 
		|| (dominazioneRichiesta && variazioneEnergia > 0)) {
		tmp = g->permutazione->individuo[posLavoroFabbricaPeggiore];
		g->permutazione->individuo[posLavoroFabbricaPeggiore] = g->permutazione->individuo[posLavoroFabbricaScelta];
		g->permutazione->individuo[posLavoroFabbricaScelta] = tmp;
	}
	else {

		individuo.punteggio.y += variazioneEnergia;
		infoFabbricaPeggiore.valutazione = val1;
		infoFabbricaScelta.valutazione = val2;

		double makespanPeggiore = 0, makespan;
		for (unsigned short i = 0; i < infoFabbriche.x.size(); i++) {
			makespan = infoFabbriche.x[i].valutazione.x;
			if (makespan > makespanPeggiore) {
				makespanPeggiore = makespan;
				individuo.punteggio.x = makespanPeggiore;
			}
		}
	}

	valutazioniEffettuate++;
}


void MOADE_DPFSPE::ottimizza(vector<Individuo>& popolazione, unsigned int& numerValutazioni) {
	for (unsigned short i = 0; i < popolazione.size(); i++) {
		ottimizzaEnergia(popolazione[i], numerValutazioni);
	}
}

void MOADE_DPFSPE::ottimizzaEnergia(Individuo& individuo, unsigned int& numeroValutazioni) {

	vector<Coppia<unsigned short>> posFabbriche = calcolaPosFabbriche(individuo.rappresentazione);

	for (unsigned short i = 0; i < posFabbriche.size(); i++) {
		ottimizzaEnergiaParziale(individuo, posFabbriche[i]);
	}

	valutaIndividuo(individuo, numeroValutazioni);
}

void MOADE_DPFSPE::ottimizzaEnergiaParziale(Individuo& individuo, Coppia<unsigned short>& posFabbriche) {

	unsigned short pos, vel;
	double tempoReale, tempoEstensione, tempoDelta;

	GruppoPDZN* g = individuo.rappresentazione; 
	unsigned short* fabbrica = g->permutazione->individuo + posFabbriche.x;

	for (unsigned short j = 0; j < istanza.macchine; j++) {
		for (unsigned short i = 0; i < posFabbriche.y; i++) {
			
			pos = fabbrica[i] * istanza.macchine + j;
			tempoReale = istanza.tp[fabbrica[i]][j] / istanza.velocita[g->modulo2->individuo[pos]] + istanza.s[fabbrica[i]][j];

			if (i == 0 && j == 0) {
				o[0][0].x = 0;
				o[0][0].y = tempoReale;
			}
			else if (i == 0 && j != 0) {
				o[0][j].x = o[0][j - 1].y;
				o[0][j].y = o[0][j].x + tempoReale;
			}
			else if (j == 0) {
				o[i][0].x = o[i - 1][0].y;
				o[i][0].y = o[i][0].x + tempoReale;
			}
			else {
				o[i][j].x = max(o[i][j - 1].y, o[i - 1][j].y);
				o[i][j].y = o[i][j].x + tempoReale;
			}

			if (i != posFabbriche.y - 1 && j != 0) {

				pos = fabbrica[i] * istanza.macchine + j - 1;

				vel = g->modulo2->individuo[pos];

				if (vel > 0) {
					tempoEstensione = min(o[i + 1][j - 1].x - o[i][j - 1].y, o[i][j].x - o[i][j - 1].y);

					if (tempoEstensione != 0) {
						tempoDelta = istanza.tp[fabbrica[i]][j - 1] * (1 / istanza.velocita[vel - 1] - 1 / istanza.velocita[vel]);
						if (tempoEstensione >= tempoDelta) {
							g->modulo2->individuo[pos]--;
						}
					}
				}
			}
		}
	}
}

void MOADE_DPFSPE::ENEH(Individuo& individuo, unsigned int& valutazioniEffettuate) {

	GruppoPDZN* g = individuo.rappresentazione;

	double* p = new double[istanza.lavori];
	unsigned short* ordinamento = new unsigned short[istanza.lavori];

	unsigned short pos;
	//Sommiamo i pij, per i=1...lavori, j = 1...macchine
	for (unsigned short i = 0; i < istanza.lavori; i++) {
		p[i] = 0;
		ordinamento[i] = i;
		for (unsigned short j = 0; j < istanza.macchine; j++) {
			pos = i * istanza.macchine + j;
			p[i] += istanza.tp[i][j] / istanza.velocita[g->modulo2->individuo[pos]] + istanza.s[i][j];
		}
	}

	//Ordiniamo i job in ordine decrescente (rispetto ai pi)
	sort(ordinamento, ordinamento + istanza.lavori, [p](unsigned short i, unsigned short j) {
		return p[i] > p[j];
	});

	//Inizializziamo l'individuo che contiene i primi f lavori
	for (unsigned short f = 0; f < istanza.fabbriche; f++) {
		g->permutazione->individuo[f] = ordinamento[f];

		if(f != istanza.fabbriche - 1)
			g->modulo1->individuo[f] = f;
	}

	unsigned short lunghezzaIndividuo = istanza.fabbriche;

	//Inizializziamo gli indici delle fabbriche

	vector<Coppia<unsigned short>> indiciFabbriche;
	indiciFabbriche.reserve(istanza.fabbriche);

	for (unsigned short k = 0; k < istanza.fabbriche; k++) {
		indiciFabbriche.push_back({k, 1});
	}

	for (unsigned short i = istanza.fabbriche; i < istanza.lavori; i++) {
		int fMigliore = -1;
		CoppiaM<double, unsigned short> migliore = { DBL_MAX, 0 };

		for (unsigned short f = 0; f < istanza.fabbriche; f++) {
			//inserisci ordinamento[i] in ogni posizione della fabbrica ff
			auto info = miglioreInserzione(g, indiciFabbriche[f].x, indiciFabbriche[f].y, ordinamento[i]);

			//valuta se è stata trovata la migliore inserzione possibile fino ad ora ed eventualmente
			//registra fabbrica e posizione migliore
			if (info.x < migliore.x) {
				migliore = info;
				fMigliore = f;
			}
		}

		//Inserisci ordinamento[i] nella posizione migliore della fabbrica migliore mai trovata
		for (unsigned short i = lunghezzaIndividuo++; i > indiciFabbriche[fMigliore].x + migliore.y; i--) {
			g->permutazione->individuo[i] = g->permutazione->individuo[i - 1];
		}
		g->permutazione->individuo[indiciFabbriche[fMigliore].x + migliore.y] = ordinamento[i];

		//Aggiorna indici fabbriche
		for (unsigned short f = fMigliore + 1; f < istanza.fabbriche; f++) {
			indiciFabbriche[f].x++;
		}

		indiciFabbriche[fMigliore].y++;

		//Aggiorna tagli 
		
		for (unsigned short f = fMigliore; f < istanza.fabbriche - 1; f++) {
			g->modulo1->individuo[f]++;
		}
		
	}

	ottimizzaEnergia(individuo, valutazioniEffettuate);
	
	delete[] p;
	delete[] ordinamento;
}


CoppiaM<double, unsigned short> MOADE_DPFSPE::miglioreInserzione(GruppoPDZN* g, unsigned short inizioFabbrica, 
	unsigned short lunghezzaFabbrica, unsigned short lavoroDaInserire) {

	if (lunghezzaFabbrica == 0) {
		unsigned short fabbrica[1] = { lavoroDaInserire };
		return {valutaIndividuoParziale(g, inizioFabbrica, 1).x , 0};
	}

	unsigned short* fabbrica = g->permutazione->individuo + inizioFabbrica;

	//Accelerazione di Taillard, guardare l'articolo relativo
	double* e = new double[istanza.macchine];
	double* f = new double[istanza.macchine];
	double** q = new double* [istanza.lavori];

	double miglioreMakespan = DBL_MAX;
	unsigned short migliorePosizione = 0;

	double tpReale;
	unsigned short pos;

	for (int i = lunghezzaFabbrica - 1; i >= 0; i--) {
		q[i] = new double[istanza.macchine];

		for (int j = istanza.macchine - 1; j >= 0; j--) {

			pos = fabbrica[i] * istanza.macchine + j;
			tpReale = istanza.tp[fabbrica[i]][j] / istanza.velocita[g->modulo2->individuo[pos]];
			tpReale += istanza.s[fabbrica[i]][j];

			if (i != lunghezzaFabbrica - 1 && j != istanza.macchine - 1) {
				q[i][j] = max(q[i][j + 1], q[i + 1][j]) + tpReale;
			}
			else if (i == lunghezzaFabbrica - 1 && j != istanza.macchine - 1) {
				q[i][j] = q[i][j + 1] + tpReale;
			}
			else if (i != lunghezzaFabbrica - 1 && j == istanza.macchine - 1) {
				q[i][j] = q[i + 1][j] + tpReale;
			}
			else {
				q[i][j] = tpReale;
			}
		}
	}

	double tpReale2;
	unsigned short pos2;

	for (unsigned short i = 0; i <= lunghezzaFabbrica; i++) {

		for (unsigned short j = 0; j < istanza.macchine; j++) {
			
			if (i != lunghezzaFabbrica) {
				pos = fabbrica[i] * istanza.macchine + j;
				tpReale = istanza.tp[fabbrica[i]][j] / istanza.velocita[g->modulo2->individuo[pos]];
				tpReale += istanza.s[fabbrica[i]][j];
			}

			pos2 = lavoroDaInserire * istanza.macchine + j;
			tpReale2 = istanza.tp[lavoroDaInserire][j] / istanza.velocita[g->modulo2->individuo[pos2]];
			tpReale2 += istanza.s[lavoroDaInserire][j];

			if (i == 0 && j == 0) {
				f[0] = tpReale2;
				e[0] = tpReale;
			}
			else if (i == 0 && j != 0) {
				f[j] = f[j - 1] + tpReale2;
				e[j] = e[j - 1] + tpReale;
			}
			else if (j == 0) {
				f[0] = e[0] + tpReale2;

				if (i != lunghezzaFabbrica)
					e[0] = e[0] + tpReale;
			}
			else {
				f[j] = max(f[j - 1], e[j]) + tpReale2;
				if (i != lunghezzaFabbrica)
					e[j] = max(e[j - 1], e[j]) + tpReale;
			}
		}

		double makeSpanParziale = 0;

		for (unsigned short j = 0; j < istanza.macchine; j++) {

			if (i != lunghezzaFabbrica)
				makeSpanParziale = max(makeSpanParziale, f[j] + q[i][j]);
			else
				makeSpanParziale = max(makeSpanParziale, f[j]);
		}

		if (makeSpanParziale < miglioreMakespan) {
			miglioreMakespan = makeSpanParziale;
			migliorePosizione = i;
		}

	}

	delete[] e;
	delete[] f;

	for (unsigned short i = 0; i < lunghezzaFabbrica; i++) {
		delete[] q[i];
	}

	delete[] q;

	return { miglioreMakespan, migliorePosizione };
}

vector<Coppia<unsigned short>> MOADE_DPFSPE::calcolaPosFabbriche(GruppoPDZN* g) {

	vector<Coppia<unsigned short>> posFabbriche;
	posFabbriche.reserve(istanza.fabbriche);

	unsigned short pos = 0, lunghezza;

	Modulo* m = g->modulo1;

	for (unsigned short i = 0; i < m->dimensione; i++) {

		lunghezza = m->individuo[i] - pos + 1;

		posFabbriche.push_back({ pos, lunghezza });

		pos = m->individuo[i] + 1;
	}

	posFabbriche.push_back({ pos, (unsigned short)(istanza.lavori - pos) });

	return posFabbriche;
}

CoppiaM<vector<MOADE_DPFSPE::InfoFabbrica>, unsigned short> MOADE_DPFSPE::calcolaInfoFabbriche(GruppoPDZN* g) {

	auto posFabbriche = calcolaPosFabbriche(g);

	auto info = vector<InfoFabbrica>(istanza.fabbriche);

	CoppiaM<unsigned short, double> fabbricaPeggiore = { 0, 0. };

	Coppia<double> valutazioneFabbrica;

	for (unsigned short i = 0; i < posFabbriche.size(); i++) {
		valutazioneFabbrica = valutaIndividuoParziale(g, posFabbriche[i].x, posFabbriche[i].y);

		info[i].coordinate = posFabbriche[i];
		info[i].valutazione = valutazioneFabbrica;

		if (valutazioneFabbrica.x > fabbricaPeggiore.y) {
			fabbricaPeggiore.x = i;
			fabbricaPeggiore.y = valutazioneFabbrica.x;
		}
	}

	return {info, fabbricaPeggiore.x};
}