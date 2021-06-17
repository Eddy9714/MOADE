#include "Main.h"

using namespace std;
using chrono::microseconds;
using chrono::system_clock;
using chrono::duration_cast;

//namespace fs = filesystem;
namespace fs = experimental::filesystem;

int main(int argc, char* argv[])
{

	unsigned int seed = 0;
	unsigned short H = 99;
	unsigned short T = 10;
	unsigned int numeroValutazioni = 100000 - (H + 1);
	double Fmax = 1.;
	double alphaMin = 0.4;
	double alphaMax = 0.8;
	string percorso;

	switch (argc) {
		case 8:
			seed = stoi(argv[7]);
		case 7:
			alphaMax = stoi(argv[6]);
		case 6:
			alphaMin = stoi(argv[5]);
		case 5:
			T = stoi(argv[4]);
		case 4:
			H = stoi(argv[3]);
		case 3:
			numeroValutazioni = stoi(argv[2]) - (H + 1);
		case 2:
			percorso = argv[1];
			break;
		case 1:
			cerr << "run with params [instance_path] [evaluations*] [H*] [T*] [AlphaMin*] [AlphaMax*] [seed*]" << endl;
			exit(-1);
			break;
	}

	MOADE_DPFSPE esecutore(percorso);

	auto inizio = chrono::high_resolution_clock::now();
	vector<Individuo> risultato = esecutore.esegui(H, T, numeroValutazioni, Fmax, alphaMin, alphaMax, seed);
	auto fine = chrono::high_resolution_clock::now();
	chrono::duration<double> durata = fine - inizio;

	int posizione = percorso.find_last_of(".");
	string percorsoNoExt = percorso.substr(0, posizione);

	auto tempo = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	posizione = percorsoNoExt.find_last_of("/");
	string nomeFile = percorsoNoExt.substr(posizione + 1);

	percorsoNoExt = percorsoNoExt + "-" + to_string((numeroValutazioni + H + 1)) + "_" + to_string(H) + "_" + to_string(T) 
		+ "_" + to_string(alphaMin) + "_" + to_string(alphaMax) + "_" + to_string(seed);

	string nomeReport = percorsoNoExt + "/" + to_string(tempo) + "-report.csv";

	if (!fs::exists(percorsoNoExt))
		fs::create_directory(percorsoNoExt);

	ofstream file(nomeReport, ios_base::app);

	if (file.is_open()) {
		
		file << "Evaluations;H;T;Fmax; Alpha Min;Alpha Max;seed;Time" << endl;
		file << numeroValutazioni + (H + 1) << ";" << H << ";" << T << ";" << Fmax << ";" << alphaMin << ";" << alphaMax << ";" << seed << ";" << durata.count() << endl << endl;

		for (unsigned short i = 0; i < risultato.size(); i++) {
			file << risultato[i].punteggio.x << ";";
			file << risultato[i].punteggio.y << endl;
		}

		file.close();
	}

	esecutore.eliminaPopolazione(risultato);
	return 0;
}
