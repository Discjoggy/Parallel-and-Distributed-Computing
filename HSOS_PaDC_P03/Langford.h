//============================================================================
// Name        : Langford.h
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Created on  : 17.11.2014
// Description : Algorithmus: Langford. Parallelisierung mit tbb.
//============================================================================

#ifndef LANGFORD_H_
#define LANGFORD_H_

#include "Definitions.h"
#include <iostream>
#include <tbb/atomic.h>
#include <tbb/task.h>
#include <vector>

void printBits(size_t& x, size_t& size);

size_t reverse(size_t x, size_t size);

size_t nextClearBit(bool* arr, size_t n, size_t i);

void print(bool* arr, size_t n);

void LangfordRecursive(bool* __arr, size_t lvl, size_t& size, size_t& count);

void LangfordRecursiveBit(size_t tree, size_t lvl, size_t& size, size_t& count);

/**
*  @brief  Repraesentiert eine Klasse, welche das Langford-
*  Problem rekursiv löst.
*/
class Langford : public tbb::task {
	size_t tree;
	size_t lvl;
	size_t& size;
	tbb::atomic<size_t>& count;

public:
	Langford(size_t __tree, size_t __lvl, size_t& __size, tbb::atomic<size_t>& __count);
	tbb::task* execute();
};

#endif /* LANGFORD_H_ */
