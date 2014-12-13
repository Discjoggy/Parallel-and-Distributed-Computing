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
		matrixMultSeq(C, A, B, n);
	}
	else {
		const M_SIZE_TYPE newN = n >> 1;
		// Devide & Conquer
		Matrix A11(newN, InnerArray(newN));
		Matrix A12(newN, InnerArray(newN));
		Matrix A21(newN, InnerArray(newN));
		Matrix A22(newN, InnerArray(newN));

		Matrix B11(newN, InnerArray(newN));
		Matrix B12(newN, InnerArray(newN));
		Matrix B21(newN, InnerArray(newN));
		Matrix B22(newN, InnerArray(newN));

		for (M_SIZE_TYPE i = 0; i < newN; ++i) {
			M_SIZE_TYPE iPlusNewN = i + newN;
			for (M_SIZE_TYPE j = 0; j < newN; ++j) {
				M_SIZE_TYPE jPlusNewN = j + newN;
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


		// M2 = (A21 + A22) * B11
		Matrix M2(newN, InnerArray(newN));
		Matrix tmp1M2(newN, InnerArray(newN));
		matrixAddSeq(tmp1M2, A21, A22, newN);
		spawn(*new (allocate_child()) Strassen(M2, tmp1M2, B11, newN));

		// M3 = A11 * (B12 - B22)
		Matrix M3(newN, InnerArray(newN));
		Matrix tmp1M3(newN, InnerArray(newN));
		matrixSubSeq(tmp1M3, B12, B22, newN);
		spawn(*new (allocate_child()) Strassen(M3, A11, tmp1M3, newN));

		// M4 = A22 * (B21 - B11)
		Matrix M4(newN, InnerArray(newN));
		Matrix tmp1M4(newN, InnerArray(newN));
		matrixSubSeq(tmp1M4, B21, B11, newN);
		spawn(*new (allocate_child()) Strassen(M4, A22, tmp1M4, newN));

		// M5 = (A11 + A12) * B22128
		Matrix M5(newN, InnerArray(newN));
		Matrix tmp1M5(newN, InnerArray(newN));
		matrixAddSeq(tmp1M5, A11, A12, newN);
		set_ref_count(5);
		spawn_and_wait_for_all(*new (allocate_child()) Strassen(M5, tmp1M5, B22, newN));

		for (M_SIZE_TYPE i = 0; i < newN; ++i) {
			M_SIZE_TYPE iPlusNewN = i + newN;
			for (M_SIZE_TYPE j = 0; j < newN; ++j) {
				M_SIZE_TYPE jPlusNewN 	= j + newN;
				C[i][j] 				= M4[i][j] - M5[i][j];
				C[i][jPlusNewN] 		= M3[i][j] + M5[i][j];
				C[iPlusNewN][j] 		= M2[i][j] + M4[i][j];
				C[iPlusNewN][jPlusNewN] = M3[i][j] - M2[i][j];
			}
		}

		M5.clear(); M5.shrink_to_fit();
		for (M_SIZE_TYPE i = 0; i < newN; ++i) {
			for (M_SIZE_TYPE j = 0; j < newN; ++j) {
				M2[i][j] = M3[i][j] = M4[i][j] = 0;
			}
		}


		// M1 = (A11 + A22) * (B11 + B22)
		// Reuse: M1 = M2 | tmp1M1 = tmp1M2 | tmp2M1 = tmp1M5
		matrixAddSeq(tmp1M2, A11, A22, newN);
		matrixAddSeq(tmp1M5, B11, B22, newN);
		spawn(*new (allocate_child()) Strassen(M2, tmp1M2, tmp1M5, newN));

		// M6 = (A21 - A11) * (B11 + B12)
		// Reuse: M6 = M3 | tmp1M6 = tmp1M3 | tmp2M6 = tmp2M3
		Matrix tmp2M3(newN, InnerArray(newN));
		matrixSubSeq(tmp1M3, A21, A11, newN);
		matrixAddSeq(tmp2M3, B11, B12, newN);
		spawn(*new (allocate_child()) Strassen(M3, tmp1M3, tmp2M3, newN));

		// M7 = (A12 - A22) * (B21 + B22)
		// Reuse: M7 = M4 | tmp1M7 = tmp1M5 | tmp2M7 = tmp2M4
		Matrix tmp2M4(newN, InnerArray(newN));
		matrixSubSeq(tmp1M4, A12, A22, newN);
		matrixAddSeq(tmp2M4, B21, B22, newN);
		set_ref_count(4);
		spawn_and_wait_for_all(*new (allocate_child()) Strassen(M4, tmp1M4, tmp2M4, newN));

		for (M_SIZE_TYPE i = 0; i < newN; ++i) {
			M_SIZE_TYPE iPlusNewN = i + newN;
			for (M_SIZE_TYPE j = 0; j < newN; ++j) {
				C[i][j] 				+= M2[i][j] + M4[i][j];
				C[iPlusNewN][j + newN] 	+= M2[i][j] + M3[i][j];
			}
		}

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
void strassenRecursive(Matrix& C, const Matrix& A, const Matrix& B, const M_SIZE_TYPE& n) {
	if (n <= CUT_OFF) {
		matrixMultSeq(C, A, B, n);
	}
	else {
		const M_SIZE_TYPE newN = n >> 1;
		// Devide & Conquer
		Matrix A11(newN, InnerArray(newN));
		Matrix A12(newN, InnerArray(newN));
		Matrix A21(newN, InnerArray(newN));
		Matrix A22(newN, InnerArray(newN));

		Matrix B11(newN, InnerArray(newN));
		Matrix B12(newN, InnerArray(newN));
		Matrix B21(newN, InnerArray(newN));
		Matrix B22(newN, InnerArray(newN));

		for (M_SIZE_TYPE i = 0; i < newN; ++i) {
			M_SIZE_TYPE iPlusNewN = i + newN;
			for (M_SIZE_TYPE j = 0; j < newN; ++j) {
				M_SIZE_TYPE jPlusNewN = j + newN;
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


		// M2 = (A21 + A22) * B11
		Matrix tmp1(newN, InnerArray(newN));
		Matrix M2(newN, InnerArray(newN));
		matrixAddSeq(tmp1, A21, A22, newN);
		strassenRecursive(M2, tmp1, B11, newN);

		// M3 = A11 * (B12 - B22)
		Matrix M3(newN, InnerArray(newN));
		matrixSubSeq(tmp1, B12, B22, newN);
		strassenRecursive(M3, A11, tmp1, newN);

		// M4 = A22 * (B21 - B11)
		Matrix M4(newN, InnerArray(newN));
		matrixSubSeq(tmp1, B21, B11, newN);
		strassenRecursive(M4, A22, tmp1, newN);

		// M5 = (A11 + A12) * B22128
		Matrix M5(newN, InnerArray(newN));
		matrixAddSeq(tmp1, A11, A12, newN);
		strassenRecursive(M5, tmp1, B22, newN);

		for (M_SIZE_TYPE i = 0; i < newN; ++i) {
			M_SIZE_TYPE iPlusNewN = i + newN;
			for (M_SIZE_TYPE j = 0; j < newN; ++j) {
				M_SIZE_TYPE jPlusNewN = j + newN;
				C[i][j] 				= M4[i][j] - M5[i][j];
				C[i][jPlusNewN] 		= M3[i][j] + M5[i][j];
				C[iPlusNewN][j] 		= M2[i][j] + M4[i][j];
				C[iPlusNewN][jPlusNewN] = M3[i][j] - M2[i][j];
			}
		}

		M5.clear(); M5.shrink_to_fit();
		for (M_SIZE_TYPE i = 0; i < newN; ++i) {
			for (M_SIZE_TYPE j = 0; j < newN; ++j) {
				M2[i][j] = M3[i][j] = M4[i][j] = 0;
			}
		}


		// M1 = (A11 + A22) * (B11 + B22)
		Matrix tmp2(newN, InnerArray(newN));
		matrixAddSeq(tmp1, A11, A22, newN);
		matrixAddSeq(tmp2, B11, B22, newN);
		strassenRecursive(M2, tmp1, tmp2, newN);

		// M6 = (A21 - A11) * (B11 + B12)
		matrixSubSeq(tmp1, A21, A11, newN);
		matrixAddSeq(tmp2, B11, B12, newN);
		strassenRecursive(M3, tmp1, tmp2, newN);

		// M7 = (A12 - A22) * (B21 + B22)
		matrixSubSeq(tmp1, A12, A22, newN);
		matrixAddSeq(tmp2, B21, B22, newN);
		strassenRecursive(M4, tmp1, tmp2, newN);

		for (M_SIZE_TYPE i = 0; i < newN; ++i) {
			M_SIZE_TYPE iPlusNewN = i + newN;
			for (M_SIZE_TYPE j = 0; j < newN; ++j) {
				C[i][j] 				+= M2[i][j] + M4[i][j];
				C[iPlusNewN][j + newN] 	+= M2[i][j] + M3[i][j];
			}
		}

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
