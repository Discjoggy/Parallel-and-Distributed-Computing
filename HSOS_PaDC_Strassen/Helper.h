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
#include <math.h>			// fabs
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

inline int isPowerOfTwo(const M_SIZE_TYPE& x) {
  return ((x != 0) && ((x & (~x + 1)) == x));
}

inline int show_usage(const char* name) {
    std::cerr << "Usage: " << name << " -<option> <value>\n"
			  << "Options:\n"
			  << "\t-h\tShow this help message\n"
			  << "\t-n\tDimension of the matrices (n X n)\n"
			  << "\t-c\tCut-Off\n"
			  << "\t-r\tRuns (e.g. 1111 for all)\n"
	  	  	  << "\t-t\tNumber of threads\n";
    return 1;
}

inline int init_arguments(const int argc, char* argv[]) {
	if (argc > 0) {
		if ((argc -1) % 2 != 0) {
			return show_usage(argv[0]);
		}
		for (int i = 1; i < argc; i++) {
			if (i % 2 != 0) {
				if (strlen(argv[i]) != 2) {
					return show_usage(argv[0]);
				}
				const char *options = "hnctr";
				if (strchr(options, argv[i][1]) == NULL) {
					return show_usage(argv[0]);
				}
				int tmp;
				switch (strchr(options, argv[i][1])[0]) {
				case 'h':
					return show_usage(argv[0]);
				case 'n':
					if (i + 1 >= argc) {
						return show_usage(argv[0]);
					}
					tmp = atoi(argv[i + 1]);
					if (tmp <= 0) {
						return show_usage(argv[0]);
					}
					if (!isPowerOfTwo(tmp)) {
						std::cerr << "Matrix dimension n has to be a value of power of two\n";
						return 1;
					}
					M_SIZE = tmp;
					break;
				case 'c':
					if (i + 1 >= argc) {
						return show_usage(argv[0]);
					}
					tmp = atoi(argv[i + 1]);
					if (tmp <= 0) {
						return show_usage(argv[0]);
					}
					if (!isPowerOfTwo(tmp)) {
						std::cerr << "Cut-Off c has to be a value of power of two\n";
						return 1;
					}
					CUT_OFF = tmp;
					break;
				case 'r':
					if (i + 1 >= argc || strlen(argv[i + 1]) != 4) {
						return show_usage(argv[0]);
					}
					RUN_NAIV_SEQ 		= argv[i + 1][0] - 48;
					RUN_NAIV_PAR 		= argv[i + 1][1] - 48;
					RUN_STRASSEN_SEQ 	= argv[i + 1][2] - 48;
					RUN_STRASSEN_PAR 	= argv[i + 1][3] - 48;
					break;
				case 't':
					if (i + 1 >= argc) {
						return show_usage(argv[0]);
					}
					if (atoi(argv[i + 1]) <= 0) {
						return show_usage(argv[0]);
					}
					NO_THREADS = atoi(argv[i + 1]);
					break;
				default:
					return show_usage(argv[0]);
				}
			}
		}

	}
	return 0;
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
inline void initializeRandpriomMatrix(Matrix& M, const M_SIZE_TYPE& size) {
	for (M_SIZE_TYPE i = 0; i < size; ++i) {
		for (M_SIZE_TYPE j = 0; j < size; ++j) {
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
*  @param  value  Zu setzender Wert (optional).
*/
inline void resetValuesMatrix(Matrix& M, const M_SIZE_TYPE& size, const M_VAL_TYPE& value = 0) {
	for (M_SIZE_TYPE i = 0; i < size; ++i) {
		for (M_SIZE_TYPE j = 0; j < size; ++j) {
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
	M_SIZE_TYPE mSize = M.size();
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
	for (M_SIZE_TYPE i = 0; i < mSize; ++i) {
		for (M_SIZE_TYPE j = 0; j < mSize; ++j) {
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
*  @return 0 falls beide Matrizen gleich sind, andernfalls != 0.
*/
inline int compareMatrices(const Matrix& A, const Matrix& B, const M_SIZE_TYPE& n) {
	for (M_SIZE_TYPE i = 0; i < n; ++i) {
		for (M_SIZE_TYPE j = 0; j < n; ++j) {
			if (fabs(A[i][j] - B[i][j]) > THRESHOLD) {
				std::cout << "A[" << i << "][" << j << "](" << A[i][j] << ") != B[" << i << "][" << j << "](" << B[i][j] << ") ";
				return 1;
			}
		}
	}
	return 0;
}

#endif
