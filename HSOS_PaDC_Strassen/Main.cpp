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
int main(int argc, char* argv[]) {
	NO_THREADS = tbb::task_scheduler_init::default_num_threads();
	const int result = init_arguments(argc, argv);
	if (result != 0) {
		return result;
	}
	tbb::task_scheduler_init init(NO_THREADS);
	std::cout << "Threads:\t" << NO_THREADS << "\n";
	std::cout << "Dimension:\t" << M_SIZE << " x " << M_SIZE << "\n";
	std::cout << "Cut-Off:\t" << CUT_OFF << "\n";

	tick_count t0, t1;
	Matrix A(M_SIZE, InnerArray(M_SIZE));
	Matrix B(M_SIZE, InnerArray(M_SIZE));
	Matrix C1(M_SIZE, InnerArray(M_SIZE));
	Matrix C2(M_SIZE, InnerArray(M_SIZE));

	initRandomizer();
	initializeRandpriomMatrix(A, M_SIZE);
	initializeRandpriomMatrix(B, M_SIZE);

	printMatrix(A, "A");
	printMatrix(B, "B");

	// Naiv
	if (RUN_NAIV_SEQ != 0) {
		resetValuesMatrix(C2, M_SIZE);
		t0 = tick_count::now();
		matrixMultSeq(C2, A, B, M_SIZE);
		t1 = tick_count::now();
		printMatrix(C2, "C2 = A * B");
		std::cout << "Naiv Seq:\tTime was " << (t1 - t0).seconds() << "s - Naiv\n";
	}

	// Naiv-Parallel
	if (RUN_NAIV_PAR != 0) {
		resetValuesMatrix(C2, M_SIZE);
		t0 = tick_count::now();
		parallel_for(blocked_range2d<M_SIZE_TYPE>(0, M_SIZE, 0, M_SIZE), MatrixMultPBody(C2, A, B, M_SIZE));
		t1 = tick_count::now();
		printMatrix(C2, "C2 = A * B");
		std::cout << "Naiv Par:\tTime was " << (t1 - t0).seconds() << "s - Naiv-Parallel\n";
	}

	// Strassen-Algorithmus: Non-Tasks
	if (RUN_STRASSEN_SEQ != 0) {
		resetValuesMatrix(C1, M_SIZE);
		t0 = tick_count::now();
		strassenRecursive(C1, A, B, M_SIZE);
		t1 = tick_count::now();
		printMatrix(C1, "C1 = A * B");
		std::cout << "Strassen Seq:\tTime was " << (t1 - t0).seconds() << "s - Strassen-Alg.: Non-Tasks\n";
		if (RUN_NAIV_PAR) {
			compareMatrices(C1, C2, M_SIZE) ? std::cout << "Num stability:\tSome differences!\n" : std::cout << "Num stability:\tOK\n";
		}
	}

	// Strassen-Algorithmus: Tasks
	if (RUN_STRASSEN_PAR != 0) {
		resetValuesMatrix(C1, M_SIZE);
		t0 = tick_count::now();
		task::spawn_root_and_wait(*new (task::allocate_root()) Strassen(C1, A, B, M_SIZE));
		t1 = tick_count::now();
		printMatrix(C1, "C1 = A * B");
		std::cout << "Strassen Par:\tTime was " << (t1 - t0).seconds() << "s - Strassen-Alg.: Tasks\n";
		if (RUN_NAIV_PAR) {
			compareMatrices(C1, C2, M_SIZE) ? std::cout << "Num stability:\tSome differences!\n" : std::cout << "Num stability:\tOK\n";
		}
	}

	std::cout << "\n\nEND\n" ;
	return 0;
}
