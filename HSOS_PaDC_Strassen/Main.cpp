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
#include "Log.h"
#include "Matrix.h"
#include "Strassen.h"
#include <vector>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tick_count.h>

using namespace tbb;

/**
*  @brief  Main-Methode zum Ausfuehren der verschiedenen Algorithmen.
*/
int main(/*int argc, char* argv[]*/) {
	Log::ReportingLevel() = logDEBUG4;
	tick_count t0, t1;

#if USE_SPECIFIC_THREAD_COUNT
	LOG(logDEBUG) << "Using " << NO_THREADS << " threads";
	tbb::task_scheduler_init init(NO_THREADS);
#endif

	Matrix A(std::vector<std::vector<M_VAL_TYPE> >(SIZE, std::vector<M_VAL_TYPE>(SIZE, 0)));
	Matrix B(std::vector<std::vector<M_VAL_TYPE> >(SIZE, std::vector<M_VAL_TYPE>(SIZE, 0)));
	Matrix C1(std::vector<std::vector<M_VAL_TYPE> >(SIZE, std::vector<M_VAL_TYPE>(SIZE, 0)));
	Matrix C2(std::vector<std::vector<M_VAL_TYPE> >(SIZE, std::vector<M_VAL_TYPE>(SIZE, 0)));

	initRandomizer();
	initializeRandpriomMatrix(A, SIZE);
	initializeRandpriomMatrix(B, SIZE);

	printMatrix(A, "A");
	printMatrix(B, "B");

//	// Naiv
//	resetValuesMatrix(C2, SIZE);
//	t0 = tick_count::now();
//	MatrixMultSeq(A, B, C2, 0, SIZE, 0, SIZE, 0, SIZE);
//	printMatrix(C2, "C2 = A * B");
//	t1 = tick_count::now();
//	LOG(logINFO) << "SMULT: Time was " << (t1 - t0).seconds() << "s" << " - Naiv";

	// Naiv-Parallel
	resetValuesMatrix(C2, SIZE);
	t0 = tick_count::now();
	parallel_for(blocked_range2d<size_t>(0, SIZE, 0, SIZE), MatrixMultPBody(A, B, C2, 0, SIZE));
	printMatrix(C2, "C2 = A * B");
	t1 = tick_count::now();
	LOG(logINFO) << "PMULT: Time was " << (t1 - t0).seconds() << "s" << " - Naiv-Parallel";

	// Strassen-Algorithmus: Non-Tasks
	resetValuesMatrix(C1, SIZE);
	t0 = tick_count::now();
	strassenRecursive(A, B, C1, SIZE);
	printMatrix(C1, "C1 = A * B");
	t1 = tick_count::now();
	LOG(logINFO) << "RSMULT: Time was " << (t1 - t0).seconds() << "s"<< " - Strassen-Alg.: Non-Tasks";
	compareMatrices(C1, C2, SIZE) ? std::cout << "RSMULT: Matrizen ungleich\n" : std::cout << "RSMULT: OK\n";

	// Strassen-Algorithmus: Tasks
	resetValuesMatrix(C1, SIZE);
	t0 = tick_count::now();
	task::spawn_root_and_wait(*new (task::allocate_root()) Strassen(A, B, C1, SIZE));
	printMatrix(C1, "C1 = A * B");
	t1 = tick_count::now();
	LOG(logINFO) << "PSMULT: Time was " << (t1 - t0).seconds() << "s" << " - Strassen-Alg.: Tasks";
	compareMatrices(C1, C2, SIZE) ? std::cout << "PSMULT: Matrizen ungleich\n" : std::cout << "PSMULT: OK\n";

	std::cout << "\n\nEND\n" ;
	return EXIT_SUCCESS;
}
