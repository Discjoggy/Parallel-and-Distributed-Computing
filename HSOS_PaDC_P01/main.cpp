//============================================================================
// Name        : Main.cpp
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Description : Beispiel für das Parallelisieren einer for-Schleife mit tbb.
//============================================================================

#include <iostream>
#include <tbb/tbb.h>

//#pragma warning(disable: 588)

#define DEBUG 0
#define RAND 100
#define N 1 << 26

typedef long long ll;
typedef unsigned long long ull;

using namespace tbb;

/**
 * Reduziert positive zahlen zu 1, negative zu -1.
 */
void make_binary(ll &value) {
    if (value > 0) {
        value = 1;
    }
    else if (value < 0) {
        value = -1;
    }
#if DEBUG
    std::cout << value << "\t";
#endif
}

/**
 * Kopiert die Werte vom source-Array in das destination-Array.
 */
void copy_values(ll* source, ll* destination) {
	for (ll i = 0; i < N; ++i) {
		destination[i] = source[i];
	}
}

/**
 * Generiert Zufallszahlen und initialisiert diese in einem Array.
 */
void randomize(ll* values) {
    srand(time(NULL));
    for (ull i = 0; i < N; ++i) {
        values[i] = (rand() % RAND - (RAND >> 1));
#if DEBUG
        std::cout << values[i] << "\t";
#endif
    }
#if DEBUG
    std::cout << std::endl << std::endl;
#endif
}

/**
 * Funktions-Objekt zur Parallelisierung der Prüfung und Zuweisung.
 */
class ParallelBinaryMaker {
	ll* values;
public:
	ParallelBinaryMaker(ll* a) : values(a) {}
    void operator() (const blocked_range<ll>& range) const {
        for (ll i = range.begin(); i != range.end(); i++) {
        	make_binary(values[i]);
        }
    }
};

/**
 * Main-Methode.
 */
int main(int argc, char** argv) {
	tick_count start, end;
    ll *randoms = new ll[N];
    ll *backup_randoms = new ll[N];
    randomize(backup_randoms);


    // Sequential loop
    copy_values(backup_randoms, randoms);
    start = tick_count::now();
    for (ull i = 0; i < N; ++i) {
    	make_binary(randoms[i]);
    }
    end = tick_count::now();
    std::cout << std::endl << "Sequential:\t" << (end - start).seconds() << " s" << std::endl << std::endl;


    // Parallel loop
    copy_values(backup_randoms, randoms);
    start = tick_count::now();
    parallel_for(blocked_range<ll>(0, N), ParallelBinaryMaker(randoms));
    end = tick_count::now();
    std::cout << std::endl << "Parallel:\t" << (end - start).seconds() << " s" << std::endl << std::endl;


    // Parallel loop (with Lambda-Expression)
    copy_values(backup_randoms, randoms);
    start = tick_count::now();
    parallel_for(
		blocked_range<ll>(0, N), [=](const blocked_range<ll> range) {
            for (ll i = range.begin(); i != range.end(); ++i) {
            	make_binary(randoms[i]);
    		}
		}
    );
    end = tick_count::now();
    std::cout << std::endl << "Parallel (L):\t" << (end - start).seconds() << " s" << std::endl << std::endl;


    return 0;
}
