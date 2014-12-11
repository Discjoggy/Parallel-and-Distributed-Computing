//============================================================================
// Name        : Main.cpp
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
//				 Jens Overmoeller <Jens.Overmoeller@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Created on  : 15.11.2014
// Description : Algorithmus: Strassen. Parallelisierung mit tbb.
//============================================================================

#include "Definitions.h"
#include "Helper.h"
//#include "Log.h"
#include "Matrix.h"
#include "Strassen.h"
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#include <tbb/task.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tick_count.h>
#include <vector>

using namespace tbb;

/**
*  @brief  Main-Methode zum Ausfuehren der verschiedenen Algorithmen.
*/
int main(/*int argc, char* argv[]*/) {
	std::cout << sysconf(_SC_LEVEL1_DCACHE_SIZE) / 1024 << "\n\n";
	//Log::ReportingLevel() = logDEBUG4;
	if (!isPowerOfTwo(SIZE)) {
		std::cout << "Value of SIZE is not power of two\n";
		return EXIT_FAILURE;
	}

#if USE_SPECIFIC_THREAD_COUNT
	tbb::task_scheduler_init init(NO_THREADS);
	std::cout << "Using " << NO_THREADS << " thread(s)\n";
	//LOG(logDEBUG) << "Using " << NO_THREADS << " threads";
#endif

	tick_count t0, t1;
	Matrix A(std::vector<std::vector<M_VAL_TYPE> >(SIZE, std::vector<M_VAL_TYPE>(SIZE)));
	Matrix B(std::vector<std::vector<M_VAL_TYPE> >(SIZE, std::vector<M_VAL_TYPE>(SIZE)));
	Matrix C1(std::vector<std::vector<M_VAL_TYPE> >(SIZE, std::vector<M_VAL_TYPE>(SIZE)));
	Matrix C2(std::vector<std::vector<M_VAL_TYPE> >(SIZE, std::vector<M_VAL_TYPE>(SIZE)));

	initRandomizer();
	initializeRandpriomMatrix(A, SIZE);
	initializeRandpriomMatrix(B, SIZE);

	printMatrix(A, "A");
	printMatrix(B, "B");

	// Naiv
//	resetValuesMatrix(C2, SIZE);
//	t0 = tick_count::now();
//	MatrixMultSeq(A, B, C2, 0, SIZE, 0, SIZE, 0, SIZE);
//	t1 = tick_count::now();
//	printMatrix(C2, "C2 = A * B");
//	std::cout << "SMULT:  Time was " << (t1 - t0).seconds() << "s - Naiv\n";
//	//LOG(logINFO) << "SMULT:  Time was " << (t1 - t0).seconds() << "s - Naiv";

	// Naiv-Parallel
//	resetValuesMatrix(C2, SIZE);
//	t0 = tick_count::now();
//	parallel_for(blocked_range2d<size_t>(0, SIZE, 0, SIZE), MatrixMultPBody(C2, A, B, 0, SIZE));
//	t1 = tick_count::now();
//	printMatrix(C2, "C2 = A * B");
//	std::cout << "PMULT:  Time was " << (t1 - t0).seconds() << "s - Naiv-Parallel\n";
//	//LOG(logINFO) << "PMULT:  Time was " << (t1 - t0).seconds() << "s - Naiv-Parallel";

	// Strassen-Algorithmus: Non-Tasks
	resetValuesMatrix(C1, SIZE);
	t0 = tick_count::now();
	strassenRecursive(C1, A, B, SIZE);
	t1 = tick_count::now();
	printMatrix(C1, "C1 = A * B");
	std::cout << "RSMULT: Time was " << (t1 - t0).seconds() << "s - Strassen-Alg.: Non-Tasks\n";
	//LOG(logINFO) << "RSMULT: Time was " << (t1 - t0).seconds() << "s - Strassen-Alg.: Non-Tasks";
	compareMatrices(C1, C2, SIZE) ? std::cout << "RSMULT: Matrizen ungleich\n" : std::cout << "RSMULT: OK\n";

	// Strassen-Algorithmus: Tasks
	resetValuesMatrix(C1, SIZE);
	t0 = tick_count::now();
	task::spawn_root_and_wait(*new (task::allocate_root()) Strassen(C1, A, B, SIZE));
	t1 = tick_count::now();
	printMatrix(C1, "C1 = A * B");
	std::cout << "PSMULT: Time was " << (t1 - t0).seconds() << "s - Strassen-Alg.: Tasks\n";
	//LOG(logINFO) << "PSMULT: Time was " << (t1 - t0).seconds() << "s - Strassen-Alg.: Tasks";
	compareMatrices(C1, C2, SIZE) ? std::cout << "PSMULT: Matrizen ungleich\n" : std::cout << "PSMULT: OK\n";

	std::cout << "\n\nEND\n" ;
	return EXIT_SUCCESS;
}
