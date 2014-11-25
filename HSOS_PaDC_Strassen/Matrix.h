//============================================================================
// Name        : Matrix.h
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
//				 Jens Overmoeller <Jens.Overmoeller@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Created on  : 15.11.2014
// Description : Algorithmus: Strassen. Parallelisierung mit tbb.
//============================================================================

#ifndef MATRIX_H_
#define MATRIX_H_

#include "Definitions.h"
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

/**
*  @brief  Subtrahiert Matrix B von Matrix A sequentiell.
*  @param  C  Matrix C (Ergebnismatrix).
*  @param  A  Matrix A.
*  @param  B  Matrix B.
*  @param  n  Matrixdimension (NxN).
*/
inline void MatrixSubSeq(Matrix& C, Matrix& A, Matrix& B, size_t n) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			C[i][j] = A[i][j] - B[i][j];
		}
	}
}

/**
*  @brief  Addiert Matrix B auf Matrix A sequentiell.
*  @param  C  Matrix C (Ergebnismatrix).
*  @param  A  Matrix A.
*  @param  B  Matrix B.
*  @param  n  Matrixdimension (NxN).
*/
inline void MatrixAddSeq(Matrix& C, Matrix& A, Matrix& B, size_t n) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}

/**
*  @brief  Multipliziert Matrix B mit Matrix A sequentiell.
*  @param  C      Matrix C (Ergebnismatrix).
*  @param  A      Matrix A.
*  @param  B      Matrix B.
*  @param  aFrom  Ab welcher Zeile.
*  @param  aTo    Bis welche Zeile.
*  @param  bFrom  Ab welcher Spalte.
*  @param  bTo    Bis welche Spalte.
*  @param  cFrom  Ab welcher Zeile (Partitioniertes multiplizieren).
*  @param  cTo    Bis welche Zeile (Partitioniertes multiplizieren).
*/
inline void MatrixMultSeq(Matrix& C, Matrix& A, Matrix& B, size_t aFrom, size_t aTo, size_t bFrom, size_t bTo, size_t cFrom, size_t cTo) {
	for (size_t i = aFrom; i < aTo; ++i) {
		for (size_t j = bFrom; j < bTo; ++j) {
			for (size_t k = cFrom; k < cTo; ++k) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

/**
*  @brief  Addiert vier Teilmatrizen auf eine Matrix.
*  @param  C    Matrix C (Ergebnismatrix).
*  @param  C11  Matrix C11 Teilmatrix (1. Quadrant).
*  @param  C12  Matrix C11 Teilmatrix (2. Quadrant).
*  @param  C21  Matrix C11 Teilmatrix (3. Quadrant).
*  @param  C22  Matrix C11 Teilmatrix (5. Quadrant).
*  @param  n    Matrixdimension (NxN).
*/
inline void MatrixIncreaseAndCopySeq(Matrix& C, Matrix& C11, Matrix& C12, Matrix& C21, Matrix& C22, size_t n) {
	for (size_t i = 0; i < n; ++i) {
		size_t iDif = i + n;
		for (size_t j = 0; j < n; ++j) {
			size_t jDif = j + n;
			C[i]	[j] 	+= 	C11[i][j];
			C[i]	[jDif]	+=	C12[i][j];
			C[iDif]	[j] 	+= 	C21[i][j];
			C[iDif]	[jDif]	+=	C22[i][j];
		}
	}
}

/**
*  @brief  Klasse für parallelisierbare Funktionsobjekte. Hinweis:
*  Wird von den nachfolgenden Klassen vererbt.
*/
struct MatrixPBody {
	Matrix& C;
	Matrix& A;
	Matrix& B;

	MatrixPBody(Matrix& __C, Matrix& __A, Matrix& __B) : C(__C), A(__A), B(__B) {}
};

/**
*  @brief  Funktionsobjekt zum parallelisierten Subtrahieren.
*/
struct MatrixSubPBody : public MatrixPBody {
	MatrixSubPBody(Matrix& __C, Matrix& __A, Matrix& __B) : MatrixPBody(__C, __A, __B) {}

	void operator() (const tbb::blocked_range2d<size_t>& range) const {
		for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
			for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
				C[i][j] = A[i][j] - B[i][j];
			}
		}
	}
};

/**
*  @brief  Funktionsobjekt zum parallelisierten Addieren.
*/
struct MatrixAddPBody : public MatrixPBody {
	MatrixAddPBody(Matrix& __C, Matrix& __A, Matrix& __B) : MatrixPBody(__C, __A, __B) {}

	void operator() (const tbb::blocked_range2d<size_t>& range) const {
		for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
			for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
				C[i][j] = A[i][j] + B[i][j];
			}
		}
	}
};

/**
*  @brief  Funktionsobjekt zum parallelisierten Multiplizieren.
*/
struct MatrixMultPBody : public MatrixPBody {
	size_t cFrom;
	size_t cTo;

	MatrixMultPBody(Matrix& __C, Matrix& __A, Matrix& __B, size_t __cFrom, size_t __cTo) : MatrixPBody(__C, __A, __B), cFrom(__cFrom), cTo(__cTo) {}

	void operator() (const tbb::blocked_range2d<size_t>& range) const {
		for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
			for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
				for (size_t k = cFrom; k < cTo; ++k) {	// Zeile X Spalte
					C[i][j] += A[i][k] * B[k][j];
				}
			}
		}
	}
};

inline void MatrixMult(Matrix& C, Matrix& A, Matrix& B, size_t aFrom, size_t aTo, size_t bFrom, size_t bTo, size_t cFrom, size_t cTo) {
#if USE_SEQ_IN_STRASSEN
	MatrixMultSeq(C, A, B, aFrom, aTo, bFrom, bTo, cFrom, cTo);
#else
	tbb::parallel_for(tbb::blocked_range2d<size_t>(0, cTo, 0, cTo), MatrixMultPBody(C, A, B, 0, cTo));
#endif
}

inline void MatrixSub(Matrix& C, Matrix& A, Matrix& B, size_t n) {
#if USE_SEQ_IN_STRASSEN
	MatrixSubSeq(C, A, B, n);
#else
	tbb::parallel_for(tbb::blocked_range2d<size_t>(0, n, 0, n), MatrixSubPBody(C, A, B));
#endif
}

inline void MatrixAdd(Matrix& C, Matrix& A, Matrix& B, size_t n) {
#if USE_SEQ_IN_STRASSEN
	MatrixAddSeq(C, A, B, n);
#else
	tbb::parallel_for(tbb::blocked_range2d<size_t>(0, n, 0, n), MatrixAddPBody(C, A, B));
#endif
}

#endif /* MATRIX_H_ */
