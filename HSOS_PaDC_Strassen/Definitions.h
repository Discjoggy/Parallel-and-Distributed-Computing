//============================================================================
// Name        : Main.cpp
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

typedef int M_VAL_TYPE;
typedef std::vector< std::vector<M_VAL_TYPE> > Matrix;

#define DEBUG 1							// Debuggen? (z.B. Verwendung von Consolen-Ausgaben)
#define MAX_RAND_VAL RAND_MAX / 50		// Zufallszahlen bis (21474836472147483647 / X) z.B. 50 oder 750
#define STD_WIDTH 9						// Matrixausgabe: Indexbreite
#define STD_PRECISION 5					// Matrixausgabe: Genauigkeit bei Gleitkommawerten
#define THRESHOLD 0.001					// Max. Abweichung als Ungenauigkeit der Gleitkommawerte

#define USE_SEQ_IN_STRASSEN 1			// Schleifen im Non-Task-Algorithmus sequentiell ausfuehren?
#define USE_SEQ_IN_STRASSEN_TASKS 1		// Schleifen im Task-Algorithmus sequentiell ausfuehren?

#define POW 9							// 2^POW
#define SIZE 2 << POW					// 4 - Dimension der Matrix (SIZE x SIZE)
#define SEQ_SIZE 2 << (POW / 2)			// Ab welcher Dimension soll naiver Algorithmus eingesetzt werden? Min. SEQ_SIZE x SEQ_SIZE = 4!
//#define SIZE 2 << 9					// 1024
//#define SIZE 1024						// 2^10
//#define SIZE 2048						// 2^11
//#define SIZE 4097						// 2^12
//#define SIZE 32769					// 2^15
//#define SIZE 262145					// 2^18

#define NO_THREADS 4					// Anzahl zu nutzender Threads (Falls aktiviert)
#define USE_SPECIFIC_THREAD_COUNT 0		// Soll spezifische Thread-Anzahl genutzt werden?

#endif /* DEFINITIONS_H_ */
