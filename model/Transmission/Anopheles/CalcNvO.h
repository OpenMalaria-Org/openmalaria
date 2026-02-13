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




double TestFortranCInteractions(int* AnnualPeriodPtr, int* nsporePtr,
				   int* FTestArray, int* testArraySizePtr, 
				   double* FAMatrix, int* AMatrixColLengthPtr, int* AMatrixRowLengthPtr, 
				   double* FBMatrix, int* BMatrixColLengthPtr, int* BMatrixRowLengthPtr, 
				   double* FCMatrix, int* CMatrixColLengthPtr, int* CMatrixRowLengthPtr, 
				   double* FDMatrix, int* DMatrixColLengthPtr, int* DMatrixRowLengthPtr);

double CalcInitMosqEmergeRate(double* FMosqEmergeRateVector, int* daysInYearPtr,
				int* mosqRestDurationPtr, int* EIPDurationPtr, int* nHostTypesInitPtr,
				int* nMalHostTypesInitPtr, double* popSizeInitPtr, 
				double* hostAvailabilityRateInitPtr, double* mosqSeekingDeathRatePtr,
				double* mosqSeekingDurationPtr, double* mosqProbBitingPtr,
				double* mosqProbFindRestSitePtr, double* mosqProbRestingPtr,
				double* mosqProbOvipositingPtr, double* FHumanInfectivityInitVector,
				double* FEIRInitVector, double* FMosqEmergeRateInitEstimateVector);

void CalcUpsilonOneHost(gsl_matrix** Upsilon, double* PAPtr,
		double* PAiPtr, int thetap, int eta, int mt, int tau,
		int thetas, int n, int m, const double* Ni, const double* alphai,
		double muvA, double thetad, const double* PBi, const double* PCi, const double* PDi,
		double PEi, const gsl_matrix* Kvi, char fntestentopar[]);

int CalcSvDiff_rf(const gsl_vector* x, void* p, gsl_vector* f);

void CalcSvDiff(gsl_vector* SvDiff, gsl_vector* SvfromEIR, 
			gsl_matrix** Upsilon, gsl_vector* Nv0, gsl_matrix* inv1Xtp, 
			int eta, int mt, int thetap, char fntestentopar[]);

void CalcLambda(gsl_vector** Lambda, gsl_vector* Nv0, int eta,
				int thetap, char fntestentopar[]);

void CalcXP(gsl_vector** xp, gsl_matrix** Upsilon, 
			gsl_vector** Lambda, gsl_matrix* inv1Xtp, int eta,
			int thetap, char fntestentopar[]);

void CalcSvJacobian(gsl_matrix* J, gsl_matrix** Upsilon, gsl_matrix* inv1Xtp,
                    int eta, int mt, int thetap, char fntestentopar[]);

void CalcPSTS(double* sumkplusPtr, double* sumklplus, int thetas,
			  int tau, double PA, double Pdf);

void FuncX(gsl_matrix* X, gsl_matrix** Upsilon, int t, int s, int n);

double CalcSpectralRadius(gsl_matrix* A, int n, char fntestentopar[]);

void CalcInv1minusA(gsl_matrix* inv1A, gsl_matrix* A, int n, char fntestentopar[]);

void CalSvfromEIRdata(gsl_vector* Sv, double PAi, double PBi, double Ni, 
					  gsl_vector* Xii);

double binomial(int n, int k);

void CalcCGSLMatrixFromCArray(gsl_matrix* CMatrix, double* FArray, int nCols, int nRows);

void CalcCGSLMatrixFromFortranArray(gsl_matrix* CMatrix, double* FArray, 
				int ColLength, int RowLength);

void CalcFortranArrayFromCGSLMatrix(gsl_matrix* CMatrix, double* FArray, 
				int ColLength, int RowLength);

void CalcCGSLVectorFromFortranArray(gsl_vector* CVector, double* FArray, 
				int Length);

void CalcFortranArrayFromCGSLVector(gsl_vector* CVector, double* FArray, 
				int Length);

void PrintRootFindingStateTS(size_t iter, gsl_multiroot_fsolver* srootfind, 
							 int thetap, char fnrootfindingstate[]);

void PrintParameters(char fntestentopar[], int thetap, int tau, int thetas, 
		int n, int m, double Ni, double alphai, double muvA, 
		double thetad, double PBi, double PCi, double PDi, double PEi, 
		gsl_vector* Kvi, gsl_vector* Xii, gsl_vector* Nv0guess);

void PrintUpsilon(char fntestentopar[], gsl_matrix** Upsilon, int thetap,
		int eta, double PA, double PAi, double Pdf, gsl_vector* Pdif,
		gsl_vector* Pduf);

void PrintXP(gsl_vector** xp, int eta, int thetap, char fntestentopar[]);

void PrintLambda(gsl_vector** Lambda, int eta, char fntestentopar[]);

void PrintEigenvalues(char fntestentopar[], gsl_vector_complex* eval, int n);

void PrintMatrix(char fntestentopar[], char matrixname[], gsl_matrix* A, 
				 int RowLength, int ColLength);

void PrintVector(char fntestentopar[], char vectorname[], gsl_vector* v, int n);










