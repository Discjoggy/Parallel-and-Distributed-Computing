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
inline void matrixSubSeq(Matrix& C, const Matrix& A, const Matrix& B, const M_SIZE_TYPE& n) {
	for (M_SIZE_TYPE i = 0; i < n; ++i) {
		for (M_SIZE_TYPE j = 0; j < n; ++j) {
			C[i][j] = A[i][j] - B[i][j];
		}
	}
}

/**
 *  @brief  Addiert Matrix B zu Matrix A sequentiell.
 *  @param  C  Matrix C (Ergebnismatrix).
 *  @param  A  Matrix A.
 *  @param  B  Matrix B.
 *  @param  n  Matrixdimension (NxN).
 */
inline void matrixAddSeq(Matrix& C, const Matrix& A, const Matrix& B, const M_SIZE_TYPE& n) {
	for (M_SIZE_TYPE i = 0; i < n; ++i) {
		for (M_SIZE_TYPE j = 0; j < n; ++j) {
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}

/**
 *  @brief  Multipliziert Matrix B mit Matrix A sequentiell.
 *  @param  C  Matrix C (Ergebnismatrix).
 *  @param  A  Matrix A.
 *  @param  B  Matrix B.
 *  @param  n  Matrixdimension (NxN).
 */
inline void matrixMultSeq(Matrix& C, const Matrix& A, const Matrix& B, const M_SIZE_TYPE& n) {
#if USE_IKJ
	for (M_SIZE_TYPE i = 0; i < n; ++i) {
		for (M_SIZE_TYPE k = 0; k < n; ++k) {
			const M_VAL_TYPE a = A[i][k];
			for (M_SIZE_TYPE j = 0; j < n; ++j) {
				C[i][j] += a * B[k][j];
			}
		}
	}
#else
	for (M_SIZE_TYPE i = 0; i < n; ++i) {
		for (M_SIZE_TYPE j = 0; j < n; ++j) {
			for (M_SIZE_TYPE k = 0; k < n; ++k) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
#endif
}

/**
 *  @brief  Klasse für parallelisierbare Funktionsobjekte. Hinweis:
 *  Wird von den nachfolgenden Klassen vererbt.
 */
struct MatrixPBody {
	Matrix& C;
	const Matrix& A;
	const Matrix& B;

	MatrixPBody(Matrix& __C, const Matrix& __A, const Matrix& __B) :
			C(__C), A(__A), B(__B) {
	}
};

/**
 *  @brief  Funktionsobjekt zum parallelisierten Subtrahieren.
 */
struct MatrixSubPBody: public MatrixPBody {
	MatrixSubPBody(Matrix& __C, const Matrix& __A, const Matrix& __B) : MatrixPBody(__C, __A, __B) { }

	void operator()(const tbb::blocked_range2d<M_SIZE_TYPE>& range) const {
		for (M_SIZE_TYPE i = range.rows().begin(); i != range.rows().end(); ++i) {
			for (M_SIZE_TYPE j = range.cols().begin(); j != range.cols().end(); ++j) {
				C[i][j] = A[i][j] - B[i][j];
			}
		}
	}
};

/**
 *  @brief  Funktionsobjekt zum parallelisierten Addieren.
 */
struct MatrixAddPBody: public MatrixPBody {
	MatrixAddPBody(Matrix& __C, const Matrix& __A, const Matrix& __B) : MatrixPBody(__C, __A, __B) { }

	void operator()(const tbb::blocked_range2d<M_SIZE_TYPE>& range) const {
		for (M_SIZE_TYPE i = range.rows().begin(); i != range.rows().end(); ++i) {
			for (M_SIZE_TYPE j = range.cols().begin(); j != range.cols().end(); ++j) {
				C[i][j] = A[i][j] + B[i][j];
			}
		}
	}
};

/**
 *  @brief  Funktionsobjekt zum parallelisierten Multiplizieren.
 */
struct MatrixMultPBody: public MatrixPBody {
	const M_SIZE_TYPE& n;

	MatrixMultPBody(Matrix& __C, const Matrix& __A, const Matrix& __B, const M_SIZE_TYPE& __n) : MatrixPBody(__C, __A, __B), n(__n) { }

	void operator()(const tbb::blocked_range2d<M_SIZE_TYPE>& range) const {
		for (M_SIZE_TYPE i = range.rows().begin(); i != range.rows().end(); ++i) {
			for (M_SIZE_TYPE k = 0; k < n; ++k) {
				const M_VAL_TYPE a = A[i][k];
				for (M_SIZE_TYPE j = range.cols().begin(); j != range.cols().end(); ++j) {
					C[i][j] += a * B[k][j];
				}
			}
		}
	}
};

#endif
