//============================================================================
// Name        : Main.cpp
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Created on  : 17.11.2014
// Description : Algorithmus: Langford. Parallelisierung mit tbb.
//============================================================================

#include "Definitions.h"
#include "Langford.h"
#include <iostream>
#include <tbb/atomic.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tick_count.h>

using namespace tbb;

tbb::atomic<size_t> countAtomic;

/**
*  @brief  Main-Methode zum Ausfuehren der verschiedenen Algorithmen.
*/
int main () {
	size_t ln = 12;
	size_t size = ln << 1;
	tick_count t0, t1;
	size_t count = 0;
	if ((ln - (ln >> 1)) % 2) {
		std::cout << "Fuer " << ln << " gibt es keine Loesung!";
		return 0;
	}
	std::cout << "TEST: L(2, " << ln << ")\n";

#if USE_SPECIFIC_THREAD_COUNT
	LOG(logDEBUG) << "Using " << NO_THREADS << " threads";
	tbb::task_scheduler_init init(NO_THREADS);
#endif

	count = 0;
	t0 = tick_count::now();
	LangfordRecursiveBit(0, ln, size, count);
	t1 = tick_count::now();
	std::cout << "Seq. Bit: Time was " << (t1 - t0).seconds() << "s" << " - Non-Tasks\n";
	std::cout << "Count: " << count << "\n\n";

	t0 = tick_count::now();
	task::spawn_root_and_wait(*new (task::allocate_root()) Langford(0, ln, size, countAtomic));
	t1 = tick_count::now();
	std::cout << "Par: Time was " << (t1 - t0).seconds() << "s" << " - Tasks\n";
	std::cout << "Count: " << countAtomic << "\n\n";

	std::cout << "\n\nEND\n" ;
	return 0;
}
