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

/**
*  @brief  Subtrahiert Matrix B von Matrix A sequentiell.
*  @param  A  Matrix A.
*  @param  B  Matrix B.
*  @param  B  Matrix C (Ergebnismatrix).
*  @param  n  Matrixdimension (NxN).
*/
inline void MatrixSubSeq(Matrix& A, Matrix& B, Matrix& C, size_t n) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			C[i][j] = A[i][j] - B[i][j];
		}
	}
}

/**
*  @brief  Addiert Matrix B auf Matrix A sequentiell.
*  @param  A  Matrix A.
*  @param  B  Matrix B.
*  @param  B  Matrix C (Ergebnismatrix).
*  @param  n  Matrixdimension (NxN).
*/
inline void MatrixAddSeq(Matrix& A, Matrix& B, Matrix& C, size_t n) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}

/**
*  @brief  Multipliziert Matrix B mit Matrix A sequentiell.
*  @param  A      Matrix A.
*  @param  B      Matrix B.
*  @param  B      Matrix C (Ergebnismatrix).
*  @param  aFrom  Ab welcher Zeile.
*  @param  aTo    Bis welche Zeile.
*  @param  bFrom  Ab welcher Spalte.
*  @param  bTo    Bis welche Spalte.
*  @param  cFrom  Ab welcher Zeile (Partitioniertes multiplizieren).
*  @param  cTo    Bis welche Zeile (Partitioniertes multiplizieren).
*/
inline void MatrixMultSeq(Matrix& A, Matrix& B, Matrix& C, size_t aFrom, size_t aTo, size_t bFrom, size_t bTo, size_t cFrom, size_t cTo) {
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
*  @param  AC   Matrix C (Ergebnismatrix).
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
	Matrix& A;
	Matrix& B;
	Matrix& C;

	MatrixPBody(Matrix& __A, Matrix& __B, Matrix& __C) : A(__A), B(__B), C(__C) {}
};

/**
*  @brief  Funktionsobjekt zum parallelisierten Subtrahieren.
*/
struct MatrixSubPBody : public MatrixPBody {
	MatrixSubPBody(Matrix& __A, Matrix& __B, Matrix& __C) : MatrixPBody(__A, __B, __C) {}

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
	MatrixAddPBody(Matrix& __A, Matrix& __B, Matrix& __C) : MatrixPBody(__A, __B, __C) {}

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

	MatrixMultPBody(Matrix& __A, Matrix& __B, Matrix& __C, size_t __cFrom, size_t __cTo) : MatrixPBody(__A, __B, __C), cFrom(__cFrom), cTo(__cTo) {}

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

/**
*  @brief  Funktionsobjekt zum parallelisierten Addieren von vier Teilmatrizen.
*/
struct MatrixIncreaseAndCopyBody : public MatrixPBody {
	Matrix& C21;
	Matrix& C22;
	size_t n;

	MatrixIncreaseAndCopyBody(Matrix& __C, Matrix& __C11, Matrix& __C12, Matrix& __C21, Matrix& __C22, size_t n) : MatrixPBody(__C11, __C12, __C), C21(__C21), C22(__C22), n(n)  {}

	void operator() (const tbb::blocked_range2d<size_t>& range) const {
		for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
			size_t iDif = i + n;
			for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
				size_t jDif = j + n;
				C[i]	[j] 	+= 	A[i][j];
				C[i]	[jDif]	+=	B[i][j];
				C[iDif]	[j] 	+= 	C21[i][j];
				C[iDif]	[jDif]	+=	C22[i][j];
			}
		}
	}
};

#endif /* MATRIX_H_ */
