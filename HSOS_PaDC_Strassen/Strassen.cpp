//============================================================================
// Name        : Strassen.cpp
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
//				 Jens Overmoeller <Jens.Overmoeller@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Created on  : 15.11.2014
// Description : Algorithmus: Strassen. Parallelisierung mit tbb.
//============================================================================

#include "Strassen.h"
#include "Helper.h"
#include "Matrix.h"
#include <vector>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

/**
*  @brief  Konstruktor der Strassen-Klasse.
*  @param  __A  Matrix A.
*  @param  __B  Matrix B.
*  @param  __C  Matrix C (Ergebnismatrix).
*  @param  __n  Matrixdimension (NxN).
*/
Strassen::Strassen(Matrix& __A, Matrix& __B, Matrix& __C, size_t __n) : A(__A), B(__B), C(__C), n(__n) {
}

/**
*  @brief  Vererbte und ueberschriebene Methode der tbb::task-Klasse,
*  welche die erbende Klasse in einem Task ausfuehren laesst. Der
*  Strassen-Algorithmus wird in dieser Methode rekursiv ausgefuerht.
*/
tbb::task* Strassen::execute() {
	size_t newN = n >> 1;

	if (n <= SEQ_SIZE) {
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixMultSeq(A, B, C, 0, n, 0, n, 0, n);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, n, 0, n), MatrixMultPBody(A, B, C, 0, n));
#endif
	}
	else {
		// Devide & Conquer
		Matrix A11(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix A12(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix A21(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix A22(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));

		Matrix B11(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix B12(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix B21(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix B22 (std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));

		Matrix C11(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix C12(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix C21(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix C22(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));

		for (size_t i = 0; i < newN; i++) {
			for (size_t j = 0; j < newN; j++) {
				A11[i][j] = A[i][j];
				A12[i][j] = A[i][j + newN];
				A21[i][j] = A[i + newN][j];
				A22[i][j] = A[i + newN][j + newN];

				B11[i][j] = B[i][j];
				B12[i][j] = B[i][j + newN];
				B21[i][j] = B[i + newN][j];
				B22[i][j] = B[i + newN][j + newN];

				C11[i][j] = C[i][j];
				C12[i][j] = C[i][j + newN];
				C21[i][j] = C[i + newN][j];
				C22[i][j] = C[i + newN][j + newN];
			}
		}


		// M2 = (A21 + A22) * B11
		Matrix M2(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp1M2(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixAddSeq(A21, A22, tmp1M2, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(A21, A22, tmp1M2));
#endif
		spawn(*new (allocate_child()) Strassen(tmp1M2, B11, M2, newN));

		// M3 = A11 * (B12 - B22)
		Matrix M3(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp1M3(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixSubSeq(B12, B22, tmp1M3, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixSubPBody(B12, B22, tmp1M3));
#endif
		spawn(*new (allocate_child()) Strassen(A11, tmp1M3, M3, newN));

		// M4 = A22 * (B21 - B11)
		Matrix M4(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp1M4(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixSubSeq(B21, B11, tmp1M4, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixSubPBody(B21, B11, tmp1M4));
#endif
		spawn(*new (allocate_child()) Strassen(A22, tmp1M4, M4, newN));

		// M5 = (A11 + A12) * B22128
		Matrix M5(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp1M5(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixAddSeq(A11, A12, tmp1M5, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(A11, A12, tmp1M5));
#endif
		set_ref_count(5);
		spawn_and_wait_for_all(*new (allocate_child()) Strassen(tmp1M5, B22, M5, newN));

		for (size_t i = 0; i < newN; ++i) {
			for (size_t j = 0; j < newN; ++j) {
				C11[i][j] = M4[i][j] - M5[i][j];
				C12[i][j] = M3[i][j] + M5[i][j];
				C21[i][j] = M2[i][j] + M4[i][j];
				C22[i][j] = M3[i][j] - M2[i][j];
			}
		}

		M2.clear(); M2.shrink_to_fit();
		M3.clear(); M3.shrink_to_fit();
		M4.clear(); M4.shrink_to_fit();
		M5.clear(); M5.shrink_to_fit();
		tmp1M2.clear(); tmp1M2.shrink_to_fit();
		tmp1M3.clear(); tmp1M3.shrink_to_fit();
		tmp1M4.clear(); tmp1M4.shrink_to_fit();
		tmp1M5.clear(); tmp1M5.shrink_to_fit();

		// M1 = (A11 + A22) * (B11 + B22)
		Matrix M1(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp1M1(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp2M1(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixAddSeq(A11, A22, tmp1M1, newN);
		MatrixAddSeq(B11, B22, tmp2M1, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(A11, A22, tmp1M1));
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(B11, B22, tmp2M1));
#endif
		spawn(*new (allocate_child()) Strassen(tmp1M1, tmp2M1, M1, newN));

		// M6 = (A21 - A11) * (B11 + B12)
		Matrix M6(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp1M6(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp2M6(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixSubSeq(A21, A11, tmp1M6, newN);
		MatrixAddSeq(B11, B12, tmp2M6, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixSubPBody(A21, A11, tmp1M6));
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(B11, B12, tmp2M6));
#endif
		spawn(*new (allocate_child()) Strassen(tmp1M6, tmp2M6, M6, newN));

		// M7 = (A12 - A22) * (B21 + B22)
		Matrix M7(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp1M7(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp2M7(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixSubSeq(A12, A22, tmp1M7, newN);
		MatrixAddSeq(B21, B22, tmp2M7, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixSubPBody(A12, A22, tmp1M7));
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(B21, B22, tmp2M7));
#endif
		set_ref_count(4);
		spawn_and_wait_for_all(*new (allocate_child()) Strassen(tmp1M7, tmp2M7, M7, newN));

		for (size_t i = 0; i < newN; i++) {
			for (size_t j = 0; j < newN; j++) {
				C11[i][j] += M1[i][j] + M7[i][j];
				C22[i][j] += M1[i][j] + M6[i][j];
			}
		}

		MatrixIncreaseAndCopy(C, C11, C12, C21, C22, newN);

		M1.clear(); M1.shrink_to_fit();
		M6.clear(); M6.shrink_to_fit();
		M7.clear(); M7.shrink_to_fit();

		A11.clear(); A11.shrink_to_fit();
		A12.clear(); A12.shrink_to_fit();
		A21.clear(); A21.shrink_to_fit();
		A22.clear(); A22.shrink_to_fit();

		B11.clear(); B11.shrink_to_fit();
		B12.clear(); B12.shrink_to_fit();
		B21.clear(); B21.shrink_to_fit();
		B22.clear(); B22.shrink_to_fit();

		C11.clear(); C11.shrink_to_fit();
		C12.clear(); C12.shrink_to_fit();
		C21.clear(); C21.shrink_to_fit();
		C22.clear(); C22.shrink_to_fit();

		tmp1M1.clear(); tmp1M1.shrink_to_fit();
		tmp1M6.clear(); tmp1M6.shrink_to_fit();
		tmp1M7.clear(); tmp1M7.shrink_to_fit();
		tmp2M1.clear(); tmp2M1.shrink_to_fit();
		tmp2M6.clear(); tmp2M6.shrink_to_fit();
		tmp2M7.clear(); tmp2M7.shrink_to_fit();
 	}

	return NULL;
}

/**
*  @brief  Rekursive Methode, welche zwei Matrizen nach dem
*  Strassen-Algorithmus errechnet.
*  @param  A  Matrix A.
*  @param  B  Matrix B.
*  @param  C  Matrix C (Ergebnismatrix).
*  @param  n  Matrixdimension (NxN).
*/
void strassenRecursive(Matrix& A, Matrix& B, Matrix& C, size_t n) {
	size_t newN = n >> 1;

	if (n <= SEQ_SIZE) {
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixMultSeq(A, B, C, 0, n, 0, n, 0, n);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, n, 0, n), MatrixMultPBody(A, B, C, 0, n));
#endif
	}
	else {
		// Devide & Conquer
		Matrix A11(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix A12(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix A21(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix A22(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));

		Matrix B11(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix B12(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix B21(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix B22 (std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));

		Matrix C11(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix C12(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix C21(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix C22(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));

		for (size_t i = 0; i < newN; i++) {
			for (size_t j = 0; j < newN; j++) {
				A11[i][j] = A[i][j];
				A12[i][j] = A[i][j + newN];
				A21[i][j] = A[i + newN][j];
				A22[i][j] = A[i + newN][j + newN];

				B11[i][j] = B[i][j];
				B12[i][j] = B[i][j + newN];
				B21[i][j] = B[i + newN][j];
				B22[i][j] = B[i + newN][j + newN];

				C11[i][j] = C[i][j];
				C12[i][j] = C[i][j + newN];
				C21[i][j] = C[i + newN][j];
				C22[i][j] = C[i + newN][j + newN];
			}
		}


		// M2 = (A21 + A22) * B11
		Matrix tmp1(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix M2(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixAddSeq(A21, A22, tmp1, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(A21, A22, tmp1));
#endif
		strassenRecursive(tmp1, B11, M2, newN);

		// M3 = A11 * (B12 - B22)
		Matrix M3(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixSubSeq(B12, B22, tmp1, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixSubPBody(B12, B22, tmp1));
#endif
		strassenRecursive(A11, tmp1, M3, newN);

		// M4 = A22 * (B21 - B11)
		Matrix M4(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixSubSeq(B21, B11, tmp1, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixSubPBody(B21, B11, tmp1));
#endif
		strassenRecursive(A22, tmp1, M4, newN);

		// M5 = (A11 + A12) * B22128
		Matrix M5(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixAddSeq(A11, A12, tmp1, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(A11, A12, tmp1));
#endif
		strassenRecursive(tmp1, B22, M5, newN);




		// M1 = (A11 + A22) * (B11 + B22)
		Matrix tmp2(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix M1(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixAddSeq(A11, A22, tmp1, newN);
		MatrixAddSeq(B11, B22, tmp2, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(A11, A22, tmp1));
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(B11, B22, tmp2));
#endif
		strassenRecursive(tmp1, tmp2, M1, newN);

		// M6 = (A21 - A11) * (B11 + B12)
		Matrix M6(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixSubSeq(A21, A11, tmp1, newN);
		MatrixAddSeq(B11, B12, tmp2, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixSubPBody(A21, A11, tmp1));
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(B11, B12, tmp2));
#endif
		strassenRecursive(tmp1, tmp2, M6, newN);

		// M7 = (A12 - A22) * (B21 + B22)
		Matrix M7(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
#if USE_SEQ_IN_STRASSEN_TASKS
		MatrixSubSeq(A12, A22, tmp1, newN);
		MatrixAddSeq(B21, B22, tmp2, newN);
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixSubPBody(A12, A22, tmp1));
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), MatrixAddPBody(B21, B22, tmp2));
#endif
		strassenRecursive(tmp1, tmp2, M7, newN);

		for (size_t i = 0; i < newN; ++i) {
			for (size_t j = 0; j < newN; ++j) {
				C11[i][j] = M1[i][j] + M4[i][j] - M5[i][j] + M7[i][j];
				C12[i][j] = M3[i][j] + M5[i][j];
				C21[i][j] = M2[i][j] + M4[i][j];
				C22[i][j] = M1[i][j] - M2[i][j] + M3[i][j] + M6[i][j];
			}
		}

		MatrixIncreaseAndCopy(C, C11, C12, C21, C22, newN);

		M1.clear(); M1.shrink_to_fit();
		M2.clear(); M2.shrink_to_fit();
		M3.clear(); M3.shrink_to_fit();
		M4.clear(); M4.shrink_to_fit();
		M5.clear(); M5.shrink_to_fit();
		M6.clear(); M6.shrink_to_fit();
		M7.clear(); M7.shrink_to_fit();

		A11.clear(); A11.shrink_to_fit();
		A12.clear(); A12.shrink_to_fit();
		A21.clear(); A21.shrink_to_fit();
		A22.clear(); A22.shrink_to_fit();

		B11.clear(); B11.shrink_to_fit();
		B12.clear(); B12.shrink_to_fit();
		B21.clear(); B21.shrink_to_fit();
		B22.clear(); B22.shrink_to_fit();

		C11.clear(); C11.shrink_to_fit();
		C12.clear(); C12.shrink_to_fit();
		C21.clear(); C21.shrink_to_fit();
		C22.clear(); C22.shrink_to_fit();

		tmp1.clear(); tmp1.shrink_to_fit();
		tmp2.clear(); tmp2.shrink_to_fit();
	}
}
