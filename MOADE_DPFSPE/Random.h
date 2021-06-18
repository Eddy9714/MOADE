#pragma once
#include <random>

using namespace std;

class Random {

	public:
		mt19937 gen;

		Random() {
			random_device rd;
			gen.seed(rd());
		}

		Random(unsigned int seed) {
			gen.seed(seed);
		}

		int randIntU(int a, int b) {
			uniform_int_distribution<int> distrib(a, b);
			return distrib(gen);
		}

		void dueIndiciRandom(unsigned int dimensione, int& i1, int& i2) {
			i1 = randIntU(0, dimensione - 1);
			i2 = (i1 + 1 + randIntU(0, dimensione - 2)) % dimensione;
		}

		void treIndiciRandom(unsigned int dimensione, int& i1, int& i2, int& i3) {

			i1 = (dimensione + 1 + randIntU(0, dimensione - 1)) % dimensione;
			int min, med, max;

			min = i1;
			max = dimensione;

			//2
			i2 = (min + 1 + randIntU(0, dimensione - 2)) % dimensione;
			if (i2 >= max || i2 < min)
				i2 = (i2 + 1) % dimensione;
			if (i2 < min)
			{
				med = min;
				min = i2;
			}
			else if (i2 > max)
			{
				med = max;
				max = i2;
			}
			else
				med = i2;
			//3
			i3 = (min + 1 + randIntU(0, dimensione - 3)) % dimensione;
			if (i3 < min || i3 >= max - 1)
				i3 = (i3 + 2) % dimensione;
			else if (i3 >= med) i3++; 
		}

		double cauchy(double a, double b) {
			std::cauchy_distribution<double> distrib(a, b);
			return distrib(gen);
		}

		double normale(double m, double s) {
			std::normal_distribution<double> distrib(m, s);
			return distrib(gen);
		}

		double randDouble(double a, double b) {
			uniform_real_distribution<double> distrib(a, b);
			return distrib(gen);
		}

		void impostaSeed(unsigned int seed) {
			gen.seed(seed);
		}
};