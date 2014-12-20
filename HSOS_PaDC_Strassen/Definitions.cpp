//============================================================================
// Name        : Definitions.cpp
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
//				 Jens Overmoeller <Jens.Overmoeller@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Created on  : 15.11.2014
// Description : Algorithmus: Strassen. Parallelisierung mit tbb.
//============================================================================

#include "Definitions.h"

int RUN_NAIV_SEQ 			= 0;
int RUN_NAIV_PAR 			= 0;
int RUN_STRASSEN_SEQ 		= 1;
int RUN_STRASSEN_PAR 		= 1;

M_SIZE_TYPE M_SIZE			= 4;
M_SIZE_TYPE CUT_OFF 		= 64;
unsigned int NO_THREADS		= 0;
