//============================================================================
// Name        : Strassen.h
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
//				 Jens Overmoeller <Jens.Overmoeller@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Created on  : 15.11.2014
// Description : Algorithmus: Strassen. Parallelisierung mit tbb.
//============================================================================

#ifndef STRASSEN_H_
#define STRASSEN_H_

#include "Definitions.h"
#include "Matrix.h"

/**
*  @brief  Repraesentiert eine Klasse, welche Matrixmultiplikationen
*  mihilfe des Strassen-Algorithmusses loest.
*/
class Strassen : public tbb::task {
	Matrix& C;
	const Matrix& A;
	const Matrix& B;
	const M_SIZE_TYPE& n;

public:
	Strassen(Matrix& __C, const Matrix& __A, const Matrix& __B, const M_SIZE_TYPE& __n) : C(__C), A(__A), B(__B), n(__n) { }

	tbb::task* execute();
};

void strassenRecursive(Matrix& C, const Matrix& A, const Matrix& B, const M_SIZE_TYPE& n);

#endif
