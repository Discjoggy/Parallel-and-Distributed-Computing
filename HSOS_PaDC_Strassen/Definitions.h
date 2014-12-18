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

#include <vector>
#include <valarray>
#include <stdint.h>

typedef uint_fast32_t M_SIZE_TYPE;					// Groessentyp der Matrizen, Schleifenzaehler usw.
typedef double M_VAL_TYPE;							// Typ der Werte in den Matrizen (Gut: int_least32_t)

#define ALL_TYPE 'S'								// Allocator-Typ
#if ALL_TYPE == 'S'
#include <tbb/scalable_allocator.h>
typedef std::vector<M_VAL_TYPE, tbb::scalable_allocator<M_VAL_TYPE>> InnerArray;
typedef std::vector<InnerArray, tbb::scalable_allocator<InnerArray>> Matrix;
#elif ALL_TYPE == 'C'
#include <tbb/cache_aligned_allocator.h>
typedef std::vector<M_VAL_TYPE, tbb::cache_aligned_allocator<M_VAL_TYPE>> InnerArray;
typedef std::vector<InnerArray, tbb::cache_aligned_allocator<InnerArray>> Matrix;
#else
typedef std::vector<M_VAL_TYPE> InnerArray;
typedef std::vector<InnerArray> Matrix;
#endif

#define DEBUG 1										// Debuggen? (Z. B. Verwendung von Consolen-Ausgaben, Konstanten Werten usw.)
#define MAX_RAND_VAL RAND_MAX / 50					// Zufallszahlen bis (21474836472147483647 / X) z.B. 50 oder 750
#define STD_WIDTH 9									// Matrixausgabe: Indexbreite
#define STD_PRECISION 5								// Matrixausgabe: Genauigkeit bei Gleitkommawerten
#define THRESHOLD 0.001								// Max. Abweichung als Ungenauigkeit der Gleitkommawerte
#define NO_BITWISE 1								// Keine arithmetischen operationen mit Shiftings loesen (Nur Ganzzahlen!)

static int RUN_NAIV_SEQ 		= 0;				// Naiven Algorithmus sequentiell ausfuehren
static int RUN_NAIV_PAR 		= 0;				// Naiven Algorithmus parallel ausfuehren
static int RUN_STRASSEN_SEQ 	= 1;				// Strassen-Algorithmus sequentiell ausfuehren
static int RUN_STRASSEN_PAR 	= 1;				// Strassen-Algorithmus parallel ausfuehren

#define POW 3										// 2^POW
static M_SIZE_TYPE M_SIZE	= (2 << POW);			// Dimension der Matrix (SIZE x SIZE)
static M_SIZE_TYPE CUT_OFF 	= 64;					// Ab welcher Dimension soll naiver Algorithmus eingesetzt werden? Min. CUT_OFF x CUT_OFF = 4! //(2 << (POW >> 1));
static unsigned NO_THREADS	= 0;					// Anzahl zu nutzender Threads

#endif
