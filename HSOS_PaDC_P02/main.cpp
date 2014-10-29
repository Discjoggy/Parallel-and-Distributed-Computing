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
#include <tbb/task_scheduler_init.h>
#include <tbb/tick_count.h>

#define DEBUG 0
#define USE_SEQ 0
#define N 1000000000  // 10^9  =  50.847.534
//#define N 10000000000 // 10^10 = 455.052.511

typedef long long ll;
typedef unsigned long long ull;

using namespace tbb;

/**
 * Gibt den übergebenen Vector in der Console aus.
 * Gibt darüber hinaus nur Primzahlen aus, solange der Parameter all false ist.
 */
void print_vector(concurrent_vector<bool> &primes, bool all = false) {
	ull primeCount = 0;
	for (ull i = 0; i < N; i++) {
		if (primes[i]) {
			primeCount++;
#if DEBUG
			std::cout << i + 2 << " ";
#endif
		}
		else if (all) {
			primeCount++;
#if DEBUG
			std::cout << i + 2 << " ";
#endif
		}
	}
	std::cout << "primeCount: " << primeCount << std::endl;
}

/**
 * Setzt alle folgenden Mehrfache von i auf im Vector an der Stelle i auf false.
 */
void eliminate_primes(concurrent_vector<bool> &primes, ull i) {
	if (primes[i - 2]) {
		for (ull j = 2 * i; j <= N; j += i) {
			primes[j - 2] = false;
		}
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
	//task_scheduler_init init(2);
	tick_count start, end;
	tick_count::interval_t dif1, dif2;
	concurrent_vector<bool> primes(N, true);

	// Primes (Sequential)
#if USE_SEQ
	for (ull i = 2; i * i <= N; i++) {
		eliminate_primes(primes, i);
	}
	end= tick_count::now();
	dif1 = end - start;
	print_vector(primes, 0);
	std::cout << std::endl << "Sequential:\t" << dif1.seconds() << " s, primeCount: \t" << primeCount << std::endl << std::endl;
#endif

	// Primes (Parallel)
	start = tick_count::now();
	parallel_for(blocked_range<ull>(2, N), ParallelEratosthenes(primes));
	end = tick_count::now();
	dif2 = end - start;
	print_vector(primes, 0);
	std::cout << std::endl << "Parallel:\t" << dif2.seconds() << " s" << std::endl << std::endl;

    return 0;
}
