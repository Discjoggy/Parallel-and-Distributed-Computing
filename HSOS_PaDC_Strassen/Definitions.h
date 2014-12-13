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
#include <stdint.h>
#include <tbb/scalable_allocator.h>

typedef uint_fast32_t M_SIZE_TYPE;					// Groessentyp der Matrizen, Schleifenzaehler usw.
typedef double M_VAL_TYPE;							// Typ der Werte in den Matrizen (Gut: int_least32_t)
typedef std::vector<M_VAL_TYPE, tbb::scalable_allocator<M_VAL_TYPE>> InnerArray;
typedef std::vector<InnerArray, tbb::scalable_allocator<InnerArray>> Matrix;

#define DEBUG 0										// Debuggen? (Z. B. Verwendung von Consolen-Ausgaben, Konstanten Werten usw.)
#define MAX_RAND_VAL RAND_MAX / 50					// Zufallszahlen bis (21474836472147483647 / X) z.B. 50 oder 750
#define STD_WIDTH 9									// Matrixausgabe: Indexbreite
#define STD_PRECISION 5								// Matrixausgabe: Genauigkeit bei Gleitkommawerten
#define THRESHOLD 0.001								// Max. Abweichung als Ungenauigkeit der Gleitkommawerte
#define NO_BITWISE 1								// Keine arithmetischen operationen mit Shiftings loesen (Nur Ganzzahlen!)

#define NO_THREADS 2								// Anzahl zu nutzender Threads (Falls aktiviert)
#define USE_SPECIFIC_THREAD_COUNT 0					// Soll spezifische Thread-Anzahl genutzt werden?

#define POW 10										// 2^POW
const M_SIZE_TYPE M_SIZE	= (2 << POW);			// Dimension der Matrix (SIZE x SIZE)
const M_SIZE_TYPE CUT_OFF 	= 64;//(2 << (POW >> 1));	// Ab welcher Dimension soll naiver Algorithmus eingesetzt werden? Min. CUT_OFF x CUT_OFF = 4!

#endif
