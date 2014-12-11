//============================================================================
// Name        : Helper.h
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
//				 Jens Overmoeller <Jens.Overmoeller@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Created on  : 15.11.2014
// Description : Algorithmus: Strassen. Parallelisierung mit tbb.
//============================================================================

#ifndef SRC_HELPER_H_
#define SRC_HELPER_H_

#include "Definitions.h"
#include <iomanip> 			// Formatierung für Matrix-Ausgabe
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>			// fabs

static inline int isPowerOfTwo(const size_t& x) {
  return ((x != 0) && ((x & (~x + 1)) == x));
}

/**
*  @brief  Initialisiert den Zufallsgenerator.
*/
inline void initRandomizer() {
	time_t t;
	time(&t);
	srand((unsigned int) t);
}

/**
*  @brief  Initialisiert eine Matrix mit Zufallswerten.
*  @param     M  Matrix M.
*  @param  size  Matrixdimension (NxN).
*/
inline void initializeRandpriomMatrix(Matrix& M, const size_t& size) {
	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j) {
#if DEBUG
			M[i][j] = 2;
#else
			M[i][j] = (M_VAL_TYPE) rand() / (MAX_RAND_VAL);
#endif
		}
	}
}

/**
*  @brief  Setzt alle Werte einer Matrix auf den uebergebenen Wert.
*  @param      M  Matrix M.
*  @param   size  Matrixdimension (NxN).
*  @param  value  Zu setzender Wert (optional)
*/
inline void resetValuesMatrix(Matrix& M, const size_t& size, const M_VAL_TYPE& value = 0) {
	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j) {
			M[i][j] = value;
		}
	}
}

/**
*  @brief  Vergleicht zwei Matrizen miteinander. Weicht ein Differenzwert
*  staerker als die THRESHOLD definierte Konstante ab, so werden die
*  Matrizen als ungleich angesehen.
*  @param           M  Matrix M.
*  @param  matrixName  Name der Matrix (optional).
*/
inline void printMatrix(const Matrix& M, const char* matrixName = NULL) {
	std::cout << "\n";
	size_t mSize = M.size();
	if (M.size() > 20) {
#if DEBUG
		std::cout << "Dimension of matrix ";
		if (matrixName != NULL) {
			std::cout << matrixName;
		}
		std::cout << " are bigger than 20x20!\n";
#endif
		return;
	}
	if (matrixName != NULL) {
		std::cout << "Matrix " << matrixName << ":\n";
	}
	for (size_t i = 0; i < mSize; ++i) {
		for (size_t j = 0; j < mSize; ++j) {
			std::cout << std::setw(STD_WIDTH) << std::setfill(' ') << std::setprecision(STD_PRECISION) << static_cast<M_VAL_TYPE>(M[i][j]);
		}
		std::cout << "\n";
	}
}

/**
*  @brief  Vergleich zwei Matrizen miteinander. Weicht ein Differenzwert
*  staerker als die THRESHOLD definierte Konstante ab, so werden die
*  Matrizen als ungleich angesehen.
*  @param  A  Matrix A.
*  @param  B  Matrix B.
*  @param  n  Matrixdimension (NxN).
*  @return 0 falls beide Matrizen gleich sind, andernfalls < 0
*/
inline int compareMatrices(const Matrix& A, const Matrix& B, const size_t& n) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			if (fabs(A[i][j] - B[i][j]) > THRESHOLD) {
				std::cout << A[i][j] << ":" << B[i][j] << "\n";
				return -1;
			}
		}
	}
	return 0;
}

#endif
