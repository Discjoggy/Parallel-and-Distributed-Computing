//============================================================================
// Name        : Main.cpp
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Description : Sieb des Eratosthenes. Umsetzung mit OpenMP.
//============================================================================

#include <iostream>
#include <math.h>
#include <omp.h>				// OpenMP
#include <tbb/tick_count.h>

typedef unsigned long Number;

//             10 =>           4
//            100 =>          25
//          1,000 =>         168
//         10,000 =>       1,229
//        100,000 =>       9,592
//      1,000,000 =>      78,498
//     10,000,000 =>     664,579
//    100,000,000 =>   5,761,455
//  1,000,000,000 =>  50,847,534
// 10,000,000,000 => 455,052,511
#define N 1000000000

Number parallel_eratosthenes(Number lastNumber) {
	// enable/disable OpenMP
	omp_set_num_threads(omp_get_num_procs());

	// instead of i * i <= lastNumber we write i <= lastNumberSquareRoot to help OpenMP
	const Number lastNumberSqrt = (Number) sqrt((double) lastNumber);
	Number memorySize = (lastNumber - 1) >> 1;
	bool* isPrime = new bool[memorySize + 1];

	#pragma omp parallel for
	for (Number i = 0; i <= memorySize; ++i) {
		isPrime[i] = 1;
	}

	#pragma omp parallel for schedule(dynamic)
	for (Number i = 3; i <= lastNumberSqrt; i += 2) {
		if (isPrime[i >> 1]) {
			for (Number j = i * i; j <= lastNumber; j += (i << 1)) {
				isPrime[j >> 1] = 0;
			}
		}
	}

	// Should be atomic, but this is already given in pragma-for-clause
	Number found = lastNumber >= 2 ? 1 : 0;
	#pragma omp parallel for reduction(+:found)
	for (Number i = 1; i <= memorySize; ++i) {
		found += isPrime[i];
	}

	delete[] isPrime;
	return found;
}

/**
 * Main-Methode.
 */
int main(int argc, char** argv) {
	tbb::tick_count t0, t1;

	t0 = tbb::tick_count::now();
	Number primes = parallel_eratosthenes(N);
	t1 = tbb::tick_count::now();
	std::cout << "\nParallel:\t" << (t1 - t0).seconds() << "s, primeCount: \t" << primes << "\n\n";

	return 0;
}
