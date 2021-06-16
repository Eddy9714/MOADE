#pragma once
#include "Random.h"

extern Random genRand;

template <typename T>
struct Coppia {
	T x;
	T y;
};

template<typename T1, typename T2>
struct CoppiaM {
	T1 x;
	T2 y;
};

template <typename T>
struct Tripla {
	T x;
	T y;
	T z;
};