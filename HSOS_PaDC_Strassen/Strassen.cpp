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

/**
*  @brief  Vererbte und ueberschriebene Methode der tbb::task-Klasse,
*  welche die erbende Klasse in einem Task ausfuehren laesst. Der
*  Strassen-Algorithmus wird in dieser Methode rekursiv ausgefuerht.
*  @return tbb::task.
*/
tbb::task* Strassen::execute() {
	if (n <= CUT_OFF) {
		MatrixMult(C, A, B, 0, n, 0, n, 0, n);
	}
	else {
		const size_t newN = n >> 1;
		// Devide & Conquer
		Matrix A11(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix A12(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix A21(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix A22(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));

		Matrix B11(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix B12(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix B21(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix B22(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));

#if USE_SEQ_IN_STRASSEN
		for (size_t i = 0; i < newN; ++i) {
			size_t iPlusNewN = i + newN;
			for (size_t j = 0; j < newN; ++j) {
				size_t jPlusNewN = j + newN;
				A11[i][j] = A[i][j];
				A12[i][j] = A[i][jPlusNewN];
				A21[i][j] = A[iPlusNewN][j];
				A22[i][j] = A[iPlusNewN][jPlusNewN];
				B11[i][j] = B[i][j];
				B12[i][j] = B[i][jPlusNewN];
				B21[i][j] = B[iPlusNewN][j];
				B22[i][j] = B[iPlusNewN][jPlusNewN];
			}
		}
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), [&](const tbb::blocked_range2d<size_t>& range) {
			for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
				size_t iPlusNewN = i + newN;
				for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
					size_t jPlusNewN = j + newN;
					A11[i][j] = A[i][j];
					A12[i][j] = A[i][jPlusNewN];
					A21[i][j] = A[iPlusNewN][j];
					A22[i][j] = A[iPlusNewN][jPlusNewN];

					B11[i][j] = B[i][j];
					B12[i][j] = B[i][jPlusNewN];
					B21[i][j] = B[iPlusNewN][j];
					B22[i][j] = B[iPlusNewN][jPlusNewN];
				}
			}
		});
#endif


		// M2 = (A21 + A22) * B11
		Matrix M2(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp1M2(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		MatrixAdd(tmp1M2, A21, A22, newN);
		spawn(*new (allocate_child()) Strassen(M2, tmp1M2, B11, newN));

		// M3 = A11 * (B12 - B22)
		Matrix M3(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp1M3(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		MatrixSubSeq(tmp1M3, B12, B22, newN);
		spawn(*new (allocate_child()) Strassen(M3, A11, tmp1M3, newN));

		// M4 = A22 * (B21 - B11)
		Matrix M4(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp1M4(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		MatrixSubSeq(tmp1M4, B21, B11, newN);
		spawn(*new (allocate_child()) Strassen(M4, A22, tmp1M4, newN));

		// M5 = (A11 + A12) * B22128
		Matrix M5(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix tmp1M5(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		MatrixAddSeq(tmp1M5, A11, A12, newN);
		set_ref_count(5);
		spawn_and_wait_for_all(*new (allocate_child()) Strassen(M5, tmp1M5, B22, newN));

#if USE_SEQ_IN_STRASSEN
		for (size_t i = 0; i < newN; ++i) {
			size_t iPlusNewN = i + newN;
			for (size_t j = 0; j < newN; ++j) {
				size_t jPlusNewN = j + newN;
				C[i][j] = M4[i][j] - M5[i][j];
				C[i][jPlusNewN] = M3[i][j] + M5[i][j];
				C[iPlusNewN][j] = M2[i][j] + M4[i][j];
				C[iPlusNewN][jPlusNewN] = M3[i][j] - M2[i][j];
			}
		}
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), [&](const tbb::blocked_range2d<size_t>& range) {
			for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
				size_t iPlusNewN = i + newN;
				for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
					size_t jPlusNewN = j + newN;
					C[i][j] = M4[i][j] - M5[i][j];
					C[i][jPlusNewN] = M3[i][j] + M5[i][j];
					C[iPlusNewN][j] = M2[i][j] + M4[i][j];
					C[iPlusNewN][jPlusNewN] = M3[i][j] - M2[i][j];
				}
			}
		});
#endif

		M5.clear(); M5.shrink_to_fit();
		for (size_t i = 0; i < newN; ++i) {
			for (size_t j = 0; j < newN; ++j) {
				M2[i][j] = M3[i][j] = M4[i][j] = 0;
			}
		}

		// M1 = (A11 + A22) * (B11 + B22)
		// Reuse: M1 = M2 | tmp1M1 = tmp1M2 | tmp2M1 = tmp1M5
		MatrixAddSeq(tmp1M2, A11, A22, newN);
		MatrixAddSeq(tmp1M5, B11, B22, newN);
		spawn(*new (allocate_child()) Strassen(M2, tmp1M2, tmp1M5, newN));

		// M6 = (A21 - A11) * (B11 + B12)
		// Reuse: M6 = M3 | tmp1M6 = tmp1M3 | tmp2M6 = tmp2M3
		Matrix tmp2M3(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		MatrixSubSeq(tmp1M3, A21, A11, newN);
		MatrixAddSeq(tmp2M3, B11, B12, newN);
		spawn(*new (allocate_child()) Strassen(M3, tmp1M3, tmp2M3, newN));

		// M7 = (A12 - A22) * (B21 + B22)
		// Reuse: M7 = M4 | tmp1M7 = tmp1M5 | tmp2M7 = tmp2M4
		Matrix tmp2M4(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		MatrixSubSeq(tmp1M4, A12, A22, newN);
		MatrixAddSeq(tmp2M4, B21, B22, newN);
		set_ref_count(4);
		spawn_and_wait_for_all(*new (allocate_child()) Strassen(M4, tmp1M4, tmp2M4, newN));

#if USE_SEQ_IN_STRASSEN
		for (size_t i = 0; i < newN; ++i) {
			size_t iPlusNewN = i + newN;
			for (size_t j = 0; j < newN; ++j) {
				size_t jPlusNewN = j + newN;
				C[i][j] 				+= M2[i][j] + M4[i][j];
				C[iPlusNewN][jPlusNewN] += M2[i][j] + M3[i][j];
			}
		}
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), [&](const tbb::blocked_range2d<size_t>& range) {
			for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
				size_t iPlusNewN = i + newN;
				for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
					size_t jPlusNewN = j + newN;
					C[i][j] 				+= M1[i][j] + M7[i][j];
					C[iPlusNewN][jPlusNewN] += M1[i][j] + M6[i][j];
				}
			}
		});
#endif

		M2.clear(); M2.shrink_to_fit();
		M3.clear(); M3.shrink_to_fit();
		M4.clear(); M4.shrink_to_fit();

		A11.clear(); A11.shrink_to_fit();
		A12.clear(); A12.shrink_to_fit();
		A21.clear(); A21.shrink_to_fit();
		A22.clear(); A22.shrink_to_fit();

		B11.clear(); B11.shrink_to_fit();
		B12.clear(); B12.shrink_to_fit();
		B21.clear(); B21.shrink_to_fit();
		B22.clear(); B22.shrink_to_fit();

		tmp1M2.clear(); tmp1M2.shrink_to_fit();
		tmp1M3.clear(); tmp1M3.shrink_to_fit();
		tmp1M4.clear(); tmp1M4.shrink_to_fit();
		tmp1M5.clear(); tmp1M5.shrink_to_fit();

		tmp2M3.clear(); tmp2M3.shrink_to_fit();
		tmp2M4.clear(); tmp2M4.shrink_to_fit();
 	}
	return NULL;
}

/**
*  @brief  Rekursive Methode, welche zwei Matrizen nach dem
*  Strassen-Algorithmus errechnet.
*  @param  C  Matrix C (Ergebnismatrix).
*  @param  A  Matrix A.
*  @param  B  Matrix B.
*  @param  n  Matrixdimension (NxN).
*/
void strassenRecursive(Matrix& C, const Matrix& A, const Matrix& B, const size_t& n) {
	if (n <= CUT_OFF) {
		MatrixMult(C, A, B, 0, n, 0, n, 0, n);
	}
	else {
		const size_t newN = n >> 1;
		// Devide & Conquer
		Matrix A11(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix A12(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix A21(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix A22(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));

		Matrix B11(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix B12(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix B21(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));
		Matrix B22(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN)));

#if USE_SEQ_IN_STRASSEN
		for (size_t i = 0; i < newN; ++i) {
			size_t iPlusNewN = i + newN;
			for (size_t j = 0; j < newN; ++j) {
				size_t jPlusNewN = j + newN;
				A11[i][j] = A[i][j];
				A12[i][j] = A[i][jPlusNewN];
				A21[i][j] = A[iPlusNewN][j];
				A22[i][j] = A[iPlusNewN][jPlusNewN];

				B11[i][j] = B[i][j];
				B12[i][j] = B[i][jPlusNewN];
				B21[i][j] = B[iPlusNewN][j];
				B22[i][j] = B[iPlusNewN][jPlusNewN];
			}
		}
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), [&](const tbb::blocked_range2d<size_t>& range) {
			for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
				size_t iPlusNewN = i + newN;
				for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
					size_t jPlusNewN = j + newN;
					A11[i][j] = A[i][j];
					A12[i][j] = A[i][jPlusNewN];
					A21[i][j] = A[iPlusNewN][j];
					A22[i][j] = A[iPlusNewN][jPlusNewN];

					B11[i][j] = B[i][j];
					B12[i][j] = B[i][jPlusNewN];
					B21[i][j] = B[iPlusNewN][j];
					B22[i][j] = B[iPlusNewN][jPlusNewN];
				}
			}
		});
#endif

		// M2 = (A21 + A22) * B11
		Matrix tmp1(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		Matrix M2(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		MatrixAdd(tmp1, A21, A22, newN);
		strassenRecursive(M2, tmp1, B11, newN);

		// M3 = A11 * (B12 - B22)
		Matrix M3(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		MatrixSub(tmp1, B12, B22, newN);
		strassenRecursive(M3, A11, tmp1, newN);

		// M4 = A22 * (B21 - B11)
		Matrix M4(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		MatrixSub(tmp1, B21, B11, newN);
		strassenRecursive(M4, A22, tmp1, newN);

		// M5 = (A11 + A12) * B22128
		Matrix M5(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		MatrixAdd(tmp1, A11, A12, newN);
		strassenRecursive(M5, tmp1, B22, newN);

#if USE_SEQ_IN_STRASSEN
		for (size_t i = 0; i < newN; ++i) {
			size_t iPlusNewN = i + newN;
			for (size_t j = 0; j < newN; ++j) {
				size_t jPlusNewN = j + newN;
#if NO_BITWISE
				C[i][j] = M4[i][j] - M5[i][j];
				C[i][jPlusNewN] = M3[i][j] + M5[i][j];
				C[iPlusNewN][j] = M2[i][j] + M4[i][j];
				C[iPlusNewN][jPlusNewN] = M3[i][j] - M2[i][j];
#else
				C[i][j] = (M4[i][j] ^ M5[i][j]) &M5[i][j];
				C[i][jPlusNewN] = M3[i][j] ^ M5[i][j];
				C[iPlusNewN][j] = M2[i][j] ^ M4[i][j];
				C[iPlusNewN][jPlusNewN] = (M3[i][j] ^ M2[i][j]) &M2[i][j];
#endif
			}
		}
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), [&](const tbb::blocked_range2d<size_t>& range) {
			for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
				size_t iPlusNewN = i + newN;
				for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
					size_t jPlusNewN = j + newN;
					C[i][j] 				= M4[i][j] - M5[i][j];
					C[i][jPlusNewN] 		= M3[i][j] + M5[i][j];
					C[iPlusNewN][j] 		= M2[i][j] + M4[i][j];
					C[iPlusNewN][jPlusNewN] = M3[i][j] - M2[i][j];
				}
			}
		});
#endif

		M5.clear(); M5.shrink_to_fit();
		for (size_t i = 0; i < newN; ++i) {
			for (size_t j = 0; j < newN; ++j) {
				M2[i][j] = M3[i][j] = M4[i][j] = 0;
			}
		}

		// M1 = (A11 + A22) * (B11 + B22)
		Matrix tmp2(std::vector< std::vector<M_VAL_TYPE> >(newN, std::vector<M_VAL_TYPE>(newN, 0)));
		MatrixAdd(tmp1, A11, A22, newN);
		MatrixAdd(tmp2, B11, B22, newN);
		strassenRecursive(M2, tmp1, tmp2, newN);

		// M6 = (A21 - A11) * (B11 + B12)
		MatrixSub(tmp1, A21, A11, newN);
		MatrixAdd(tmp2, B11, B12, newN);
		strassenRecursive(M3, tmp1, tmp2, newN);

		// M7 = (A12 - A22) * (B21 + B22)
		MatrixSub(tmp1, A12, A22, newN);
		MatrixAdd(tmp2, B21, B22, newN);
		strassenRecursive(M4, tmp1, tmp2, newN);

#if USE_SEQ_IN_STRASSEN
		for (size_t i = 0; i < newN; ++i) {
			size_t iPlusNewN = i + newN;
			for (size_t j = 0; j < newN; ++j) {
				size_t jPlusNewN = j + newN;
#if NO_BITWISE
				C[i][j] 				+= M2[i][j] + M4[i][j];
				C[iPlusNewN][jPlusNewN] += M2[i][j] + M3[i][j];
#else
				C[i][j] 				^= M2[i][j] ^ M4[i][j];
				C[iPlusNewN][jPlusNewN] ^= M2[i][j] ^ M3[i][j];
#endif
			}
		}
#else
		tbb::parallel_for(tbb::blocked_range2d<size_t>(0, newN, 0, newN), [&](const tbb::blocked_range2d<size_t>& range) {
			for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
				size_t iPlusNewN = i + newN;
				for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
					size_t jPlusNewN = j + newN;
					C[i][j] 				^= M1[i][j] ^ M7[i][j];
					C[iPlusNewN][jPlusNewN] ^= M1[i][j] ^ M6[i][j];
				}
			}
		});
#endif

		M2.clear(); M2.shrink_to_fit();
		M3.clear(); M3.shrink_to_fit();
		M4.clear(); M4.shrink_to_fit();

		A11.clear(); A11.shrink_to_fit();
		A12.clear(); A12.shrink_to_fit();
		A21.clear(); A21.shrink_to_fit();
		A22.clear(); A22.shrink_to_fit();

		B11.clear(); B11.shrink_to_fit();
		B12.clear(); B12.shrink_to_fit();
		B21.clear(); B21.shrink_to_fit();
		B22.clear(); B22.shrink_to_fit();

		tmp1.clear(); tmp1.shrink_to_fit();
		tmp2.clear(); tmp2.shrink_to_fit();
	}
}
