//============================================================================
// Name        : Definitions.h
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
//				 Jens Overmoeller <Jens.Overmoeller@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Created on  : 15.11.2014
// Description : Algorithmus: Strassen. Parallelisierung mit tbb.
//============================================================================

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <stdint.h>
#include <vector>

typedef uint_fast32_t M_SIZE_TYPE;		// Groessentyp der Matrizen, Schleifenzaehler usw.
typedef double M_VAL_TYPE;				// Typ der Werte in den Matrizen (Gut: int_least32_t)

#define ALL_TYPE 'S'					// Allocator-Typ
#if ALL_TYPE == 'S'
#include <tbb/scalable_allocator.h>
typedef std::vector<M_VAL_TYPE, tbb::scalable_allocator<M_VAL_TYPE> > InnerArray;
typedef std::vector<InnerArray, tbb::scalable_allocator<InnerArray> > Matrix;
#elif ALL_TYPE == 'C'
#include <tbb/cache_aligned_allocator.h>
typedef std::vector<M_VAL_TYPE, tbb::cache_aligned_allocator<M_VAL_TYPE> > InnerArray;
typedef std::vector<InnerArray, tbb::cache_aligned_allocator<InnerArray> > Matrix;
#else
typedef std::vector<M_VAL_TYPE> InnerArray;
typedef std::vector<InnerArray> Matrix;
#endif

#define USE_PARTITIONS 1				// Aktiviert partitionierte Strassen-Algorithmen (bspw. Half-And-Half)
#define DEBUG 0							// Debuggen? (Z. B. Verwendung von Consolen-Ausgaben, Konstanten Werten usw.)
#define MAX_RAND_VAL RAND_MAX / 50		// Zufallszahlen bis (21474836472147483647 / X) z.B. 50 oder 750
#define STD_WIDTH 9						// Matrixausgabe: Indexbreite
#define STD_PRECISION 5					// Matrixausgabe: Genauigkeit bei Gleitkommawerten
#define THRESHOLD 0.001					// Max. Abweichung als Ungenauigkeit der Gleitkommawerte

// extern - globals
extern int RUN_NAIV_SEQ;				// Naiven Algorithmus sequentiell ausfuehren
extern int RUN_NAIV_PAR;				// Naiven Algorithmus parallel ausfuehren
extern int RUN_STRASSEN_SEQ;			// Strassen-Algorithmus sequentiell ausfuehren
extern int RUN_STRASSEN_PAR;			// Strassen-Algorithmus parallel ausfuehren

extern M_SIZE_TYPE M_SIZE;				// Dimension der Matrix (SIZE x SIZE)
extern M_SIZE_TYPE CUT_OFF;				// Ab welcher Dimension soll naiver Algorithmus eingesetzt werden? Min. CUT_OFF x CUT_OFF = 4! //(2 << (POW >> 1));
extern unsigned NO_THREADS;				// Anzahl zu nutzender Threads

#endif
