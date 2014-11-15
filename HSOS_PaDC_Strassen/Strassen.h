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
#include <tbb/task.h>

/**
*  @brief  Repraesentiert eine Klasse, welche Matrixmultiplikationen
*  mihilfe des Strassen-Algorithmusses loest.
*/
class Strassen : public tbb::task {
	Matrix& A;
	Matrix& B;
	Matrix& C;
	size_t n;

public:
	Strassen(Matrix& __A, Matrix& __B, Matrix& __C,size_t __n);
	tbb::task* execute();
};

void strassenRecursive(Matrix& A, Matrix& B, Matrix& C, size_t n);

#endif /* STRASSEN_H_ */
