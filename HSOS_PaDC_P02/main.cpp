//============================================================================
// Name        : Main.cpp
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Description : Sieb des Eratosthenes. Umsetzung mit tbb.
//============================================================================

#include <iostream>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/tick_count.h>

#define DEBUG 0
#define N 1 << 20

typedef long long ll;
typedef unsigned long long ull;

using namespace tbb;

/**
 * Gibt den übergebenen Vector in der Console aus.
 * Gibt darüber hinaus nur Primzahlen aus, solange der Parameter all false ist.
 */
void print_vector(concurrent_vector<bool> &primes, bool all = false) {
	for (ull i = 0; i < N; i++) {
		if (primes[i]) {
			std::cout << i + 2 << " ";
		}
		else if (all) {
			std::cout << i + 2 << " ";
		}
	}
}

/**
 * Setzt alle folgenden Mehrfache von i auf im Vector an der Stelle i auf false.
 */
void eliminate_primes(concurrent_vector<bool> &primes, ull i) {
	if (primes[i - 2]) {
		for (auto j = 2 * i; j <= N; j += i) {
			primes[j - 2] = false;
		}
	}
}

/**
 * Leert und initialisiert den übergebenen Vector in dem alle Einträge auf true gesetzt werden.
 */
void clear_and_initialize_primes(concurrent_vector<bool> &primes) {
	primes.clear();
	for (ull i = 2; i < N; i++) {
		primes.push_back(true);
	}
}

/**
 * Funktions-Objekt zur Parallelisierung des Sieb des Eratosthenes.
 */
class ParallelEratosthenes {
	concurrent_vector<bool> &primes;
public:
	ParallelEratosthenes(concurrent_vector<bool> &primes) : primes(primes) {}
    void operator() (const blocked_range<ull>& range) const {
    	for (ull i = range.begin(); i * i <= N; i++) {
    		eliminate_primes(primes, i);
    	}
    }
};

/**
 * Main-Methode.
 */
int main(int argc, char** argv) {
	tick_count start, end;
	concurrent_vector<bool> primes;


	// Primes (Sequential)
	clear_and_initialize_primes(primes);
	start = tick_count::now();
	for (ull i = 2; i * i <= N; i++) {
		eliminate_primes(primes, i);
	}
	end = tick_count::now();
#if DEBUG
	print_vector(primes, 0);
#endif
	std::cout << std::endl << "Sequential:\t" << (end - start).seconds() << " s" << std::endl << std::endl;


	// Primes (Parallel)
	clear_and_initialize_primes(primes);
	start = tick_count::now();
	parallel_for(blocked_range<ull>(2, N), ParallelEratosthenes(primes));
	end = tick_count::now();
#if DEBUG
	print_vector(primes, 0);
#endif
	std::cout << std::endl << "Parallel:\t" << (end - start).seconds() << " s" << std::endl << std::endl;

    return 0;
}
