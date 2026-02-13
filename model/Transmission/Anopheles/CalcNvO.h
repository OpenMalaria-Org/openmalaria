/* This file is part of mcdnsa.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute
 * 
 * mcdnsa is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* This file should contain the headers of all routines that we write in the C
 * program.
 */ 

/* We also include library headers here. */ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_multiroots.h>

double CalcInitMosqEmergeRate(double* FMosqEmergeRateVector, int* daysInYearPtr,
				int* mosqRestDurationPtr, int* EIPDurationPtr, int* nHostTypesInitPtr,
				int* nMalHostTypesInitPtr, double* popSizeInitPtr, 
				double* hostAvailabilityRateInitPtr, double* mosqSeekingDeathRatePtr,
				double* mosqSeekingDurationPtr, double* mosqProbBitingPtr,
				double* mosqProbFindRestSitePtr, double* mosqProbRestingPtr,
				double* mosqProbOvipositingPtr, double* FHumanInfectivityInitVector,
				double* FEIRInitVector, double* FMosqEmergeRateInitEstimateVector);

void CalcUpsilon(gsl_matrix** Upsilon, double* PAPtr,
		double* PAiPtr, int thetap, int eta, int mt, int tau,
		int thetas, int n, int m, const double* Ni, const double* alphai,
		double muvA, double thetad, const double* PBi, const double* PCi, const double* PDi,
		double PEi, const gsl_matrix* Kvi);

void CalcSvJacobian(gsl_matrix* J, gsl_matrix** Upsilon, gsl_matrix* inv1Xtp, int eta, int mt, int thetap);

void CalcPSTS(double* sumkplusPtr, double* sumklplus, int thetas, int tau, double PA, double Pdf);

void FuncX(gsl_matrix* X, gsl_matrix** Upsilon, int t, int s, int n);

double CalcSpectralRadius(gsl_matrix* A, int n);

void CalcInv1minusA(gsl_matrix* inv1A, gsl_matrix* A, int n);

double binomial(int n, int k);

void CalcCGSLMatrixFromCArray(gsl_matrix* CMatrix, double* FArray, int nCols, int nRows);

void CalcCGSLMatrixFromFortranArray(gsl_matrix* CMatrix, double* FArray, int ColLength, int RowLength);

void CalcFortranArrayFromCGSLMatrix(gsl_matrix* CMatrix, double* FArray, int ColLength, int RowLength);

void CalcCGSLVectorFromFortranArray(gsl_vector* CVector, double* FArray, int Length);

void CalcFortranArrayFromCGSLVector(gsl_vector* CVector, double* FArray, int Length);
