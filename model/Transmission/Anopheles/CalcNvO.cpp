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

/* Entomology model coordinator: Nakul Chitnis. */ 


/* We need to be very careful with pointers. Arrays passed in from
 * Fortran will be passed as pointers. Any changes made to those pointers
 * in C will also change the values in Fortran. We should make copies of
 * those pointers into arrays in C and work with them as arrays. From
 * Fortran we can also pass pointers to arrays of 0's for those arrays
 * that we wish to pass from C to Fortran. I think this will work. Let's
 * see what happens...
 */ 

/* We should also be careful between floats and doubles. Most of the 
 * variables in Fotran are defined as real - which would translate 
 * into floats. However, I think most of the C mathematics libraries
 * are probably for doubles so we should be careful about going back
 * and forth between these. Maybe make more variables real*8 in Fortran. 
 */ 


/* We are currently trying to create arrays and two dimensional arrays.
 * It may make more sense to just create gsl_vectors and gsl_matrices.
 * This may be a problem when we define Upsilon as a three dimensional
 * array. Maybe there is some way of dealing with that in gsl - but
 * perhaps not. Let's deal with that later.
 */ 

/* We use the naming convention that all arrays and matrices that come from 
 * Fortran and will be sent back to Fortran begin with 'F'. All vectors
 * and matrices that are created and used by C begin with 'C'.
 * Hopefully this will help to keep things less confusing
 * - although certainly not eliminate the confusion...
 */

/* In C, the first index refers to the column and the second index to the 
 * row. This is stupid, but we have no choice. However, in our attempt to
 * be as consistent as possible, we always refer to the row by 'i' and to
 * the column by 'j'. We refer to the element in the i^{th} row and j^{th} 
 * column of matrix, A, by A(j,i). This is certainly not perfect, but 
 * perhaps the best that we can do.
 */ 


/**************************************************************************
 ****************************  HEADERS ************************************
 **************************************************************************/

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

#include "CalcNvO.h"

using namespace std;


/***************************************************************************
 *********************** STRUCTURE DEFINITIONS *****************************
 ***************************************************************************/

// Structure that contains the parameters for the function used in the 
// root-finding algorithm to find the emergence rate that matches the 
// number of infectious host-seeking mosquitoes.
struct SvDiffParams
{
	gsl_vector* SvfromEIR;
	gsl_matrix** Upsilon;
	gsl_matrix* inv1Xtp;
	int eta;
	int mt;
	int thetap;
};



/***************************************************************************
 ************************ START SUBROUTINES HERE ***************************
 ***************************************************************************/


/****************************************************************************/
/* testFortranCInteractions() passes test arrays and matrices between C and
 * Fortran to see if we can get the communication to work correctly. For now
 * we aim to pass an array from Fortran to C, add 1 to all elements of the 
 * array, and read the new array in Fortran, while saving a copy of the
 * original array in C. 
 *
 * We also pass two matrices from Fortran to C:
 * A is a (defined) 2 x 3 matrix
 * B is a (defined) 3 x 2 matrix.
 * In C, we then evaluate:
 * C = AB (2 x 2)
 * D = A + B^T (2 x 3)
 *
 * We see if C returns to Fortran what we would expect.
 */

double TestFortranCInteractions(int* AnnualPeriodPtr, int* nsporePtr,
				   int* FTestArray, int* testArraySizePtr, 
				   double* FAMatrix, int* AMatrixColLengthPtr, int* AMatrixRowLengthPtr, 
				   double* FBMatrix, int* BMatrixColLengthPtr, int* BMatrixRowLengthPtr, 
				   double* FCMatrix, int* CMatrixColLengthPtr, int* CMatrixRowLengthPtr, 
				   double* FDMatrix, int* DMatrixColLengthPtr, int* DMatrixRowLengthPtr){
 

    // Initialize variables.

    int AnnualPeriod;
	int nspore;
	int i;
	int j;

	int testArraySize;
	int AMatrixRowLength;
	int AMatrixColLength;
	int BMatrixRowLength;
	int BMatrixColLength;
	int CMatrixRowLength;
	int CMatrixColLength;
	int DMatrixRowLength;
	int DMatrixColLength;
	

	double MosqEmergeRate;
	double temp;

	// Dereference pointers.
	AnnualPeriod = *AnnualPeriodPtr;
	nspore = *nsporePtr;

	testArraySize = *testArraySizePtr;
	AMatrixRowLength = *AMatrixRowLengthPtr;
	AMatrixColLength = *AMatrixColLengthPtr;
	BMatrixRowLength = *BMatrixRowLengthPtr;
	BMatrixColLength = *BMatrixColLengthPtr;
	CMatrixRowLength = *CMatrixRowLengthPtr;
	CMatrixColLength = *CMatrixColLengthPtr;
	DMatrixRowLength = *DMatrixRowLengthPtr;
	DMatrixColLength = *DMatrixColLengthPtr;


	/* We initialize copyTestArray and dynamically allocate memory to it. 
	 * This seems to work as we would like it to.
	 */
	int* copyTestArray;
	copyTestArray = (int *)malloc(testArraySize*sizeof(int));


	// We now initialize CAMatrix through CDMatrix.
	// I think we should always initialize gsl matrices to 0.
	gsl_matrix* CAMatrix = gsl_matrix_calloc(AMatrixColLength, AMatrixRowLength); 
	gsl_matrix* CBMatrix = gsl_matrix_calloc(BMatrixColLength, BMatrixRowLength); 
	gsl_matrix* CCMatrix = gsl_matrix_calloc(CMatrixColLength, CMatrixRowLength); 
	gsl_matrix* CDMatrix = gsl_matrix_calloc(DMatrixColLength, DMatrixRowLength); 

	// We allocate space for the transpose of B.
	// We're switching RowLength and ColLength - but we should be careful not
	// to get confused.
	gsl_matrix* CBTMatrix = gsl_matrix_calloc(BMatrixRowLength, BMatrixColLength);



	/* We copy testArray to copyTestArray.
	 * We add 1 to each element of testArray.
	 * We print both testArray and copyTestArray.
	 */ 
	for( i = 0; i < testArraySize; i++) {
		copyTestArray[i] = FTestArray[i];
		FTestArray[i]=FTestArray[i]+1;
		printf("FTestArray[%d] = %d\n", i, FTestArray[i]);
		printf("copyTestArray[%d] = %d\n", i, copyTestArray[i]);
	}



	/* We now try to define CAMatrix and CBMatrix as gsl matrices from the 
	 * Fortran arrays, FAMatrix and BFMatrix.
	 */ 
	CalcCGSLMatrixFromFortranArray(CAMatrix, FAMatrix, AMatrixColLength, AMatrixRowLength);
	CalcCGSLMatrixFromFortranArray(CBMatrix, FBMatrix, BMatrixColLength, BMatrixRowLength);


	/* Print out FAMatrix and CAMatrix in C to see what they look like.
	 * We add 1 to the indices so that they refer to the actual rows and columns.
	 */
	for( i=0; i<AMatrixColLength; i++) {
		for( j=0; j<AMatrixRowLength; j++) {
			printf("FAMatrix([%d],[%d]) = [%f]\n", i+1, j+1, FAMatrix[j*AMatrixColLength+i]);
			temp = gsl_matrix_get(CAMatrix, i, j);
			printf("CAMatrix([%d],[%d]) = [%f]\n", i+1, j+1, temp);
		}
	}

	/* Print out FBMatrix and CBMatrix in C to see what they look like.
	 * We add 1 to the indices so that they refer to the actual rows and columns.
	 */
	for( i=0; i<BMatrixColLength; i++) {
		for( j=0; j<BMatrixRowLength; j++) {
			printf("FBMatrix([%d],[%d]) = [%f]\n", i+1, j+1, FBMatrix[j*BMatrixColLength+i]);
			temp = gsl_matrix_get(CBMatrix, i, j);
			printf("CBMatrix([%d],[%d]) = [%f]\n", i+1, j+1, temp);
		}
	}
	/* Let's see what happens? Here comes the test...
	 * So far so good. It works as we would like it to.
	 */ 

	// Now we evaluate C = AB.

	// This routine sets C = 1*A*B + 0*C.
	// The first two terms tell it not to take the transpose of A and B.
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, CAMatrix, CBMatrix, 0.0, CCMatrix);

	// Now: D = A+B^T.

	// Evaluate transpose
	gsl_matrix_transpose_memcpy(CBTMatrix, CBMatrix);
	// First copy A into D.
	gsl_matrix_memcpy(CDMatrix, CAMatrix);
	// Add transpose(B) to D and save in D.
	gsl_matrix_add(CDMatrix, CBTMatrix);


	/* Save C and D matrices as Fortran arrays */ 
	CalcFortranArrayFromCGSLMatrix(CCMatrix, FCMatrix, CMatrixColLength, CMatrixRowLength);
	CalcFortranArrayFromCGSLMatrix(CDMatrix, FDMatrix, DMatrixColLength, DMatrixRowLength);


	/* Print out FCMatrix and CCMatrix in C to see what they look like.
	 * We add 1 to the indices so that they refer to the actual rows and columns.
	 */
	for( i=0; i<CMatrixColLength; i++) {
		for( j=0; j<CMatrixRowLength; j++) {
			printf("FCMatrix([%d],[%d]) = [%f]\n", i+1, j+1, FCMatrix[j*CMatrixColLength+i]);
			temp = gsl_matrix_get(CCMatrix, i, j);
			printf("CCMatrix([%d],[%d]) = [%f]\n", i+1, j+1, temp);
		}
	}

	/* Print out FDMatrix and CDMatrix in C to see what they look like.
	 * We add 1 to the indices so that they refer to the actual rows and columns.
	 */
	for( i=0; i<DMatrixColLength; i++) {
		for( j=0; j<DMatrixRowLength; j++) {
			printf("FDMatrix([%d],[%d]) = [%f]\n", i+1, j+1, FDMatrix[j*DMatrixColLength+i]);
			temp = gsl_matrix_get(CDMatrix, i, j);
			printf("CDMatrix([%d],[%d]) = [%f]\n", i+1, j+1, temp);
		}
	}


	//getchar();


	MosqEmergeRate = AnnualPeriod*nspore;

	// Deallocate memory for gsl matrices.
	free(copyTestArray);
	gsl_matrix_free(CAMatrix);
	gsl_matrix_free(CBMatrix);
	gsl_matrix_free(CCMatrix);
	gsl_matrix_free(CDMatrix);


	return MosqEmergeRate;
}
/***************************************************************************/






/***************************************************************************/

/* calcInitMosqEmergeRate() calculates the mosquito emergence rate given
 * all other parameters.
 *
 * We use a periodic version of the model described in "A Mathematical Model 
 * for the Dynamics of Malaria in Mosquitoes Feeding on a Heteregeneous Host
 * Population". The periodic model still needs to be written as a paper. We will
 * change these comments to refer to the approprirate paper when it is ready.
 *
 * The entomological model has a number of input parameters, including the
 * mosquito emergence rate, $N_{v0}$, and a number of output parameters, 
 * including the entomological inoculation rate, $\Xi_i$. The model produces
 * equations for $\Xi_i$ as a function of $N_{v0}$ and the other parameters.
 * However, in this function, we assume that all parameters, except $N_{v0}$ 
 * are known, and $\Xi_i$ is known. We then use these parameters, with $\Xi_i$ 
 * to calculate $N_{v0}$. The equations for $\Xi_i$ are linear in terms of 
 * $N_{v0}$ so there is a unique solution for $N_{v0}$. 
 *
 * This routine first shows the existence of a unique globally asymptotically 
 * stable periodic orbit for the system of equations describing the periodically
 * forced entomological model (for a given set of parameter values, including the
 * mosquito emergence rate). It then compares the number of infectious host-seeking
 * mosquitoes for this periodic orbit to the the number of infectious host-seeking
 * mosquitoes that would result in the given EIR. The routine then iteratively finds
 * the emergence rate that matches the given EIR.
 * 
 * However, we cannot write these equations in the form Ax=b, so we use
 * a root-finding algorithm to calculate $N_{v0}$.
 *
 * This function has a dummy return of 0.
 * 
 * FMosqEmergeRateVector is an OUT parameter.
 * All other parameters are IN parameters.
 */

double CalcInitMosqEmergeRate(double* FMosqEmergeRateVector, int* daysInYearPtr,
				int* mosqRestDurationPtr, int* EIPDurationPtr, int* nHostTypesInitPtr,
				int* nMalHostTypesInitPtr, double* popSizeInitPtr, 
				double* hostAvailabilityRateInitPtr, double* mosqSeekingDeathRatePtr,
				double* mosqSeekingDurationPtr, double* mosqProbBitingPtr,
				double* mosqProbFindRestSitePtr, double* mosqProbRestingPtr,
				double* mosqProbOvipositingPtr, double* FHumanInfectivityInitVector,
				double* FSvInitVector, double* FMosqEmergeRateInitEstimateVector){



    /* Note that from here on we use the notation from "A Mathematical Model for the
	 * Dynamics of Malaria in Mosquitoes Feeding on a Heterogeneous Host Population",
	 * and (the publication with the periodic model - yet to be written).
	 *
	 * While, this may not be the easiest notation to read for someone not familiar
	 * with the model, it will be easier to go directly from the equations in the paper
	 * to the equations, as they will be written in the code. Since the equations are not
	 * obvious in any case, anyone who wants to go through this code, will need to go 
	 * through the paper, so I think that will be ok.
	 *
	 * There are also a number of variables defined that are difficult to describe
	 * physically which we use in intermediate equations. We try to give names that
	 * we use in the papers referenced above. 
	 *
	 *  - Any complaints about this notation (or anything else in general) can be directed 
	 * to itsupport-sti@stimail.ch
	 *
	 * Once the paper on the periodic model is written/published - we should also include
	 * the equation numbers as that may help.l
	 *
	 * We may append a 'CV' or 'CM' to gsl_vectors and gsl_matrices to distinguish them
	 * if we feel it is necessary.
	 *
	 * As far as possible, we try to use gsl_vectors instead of arrays to allow more
	 * flexibility.
	 *
	 */

	int i;
	double temp;

	
	// We initialize the parameters below. Where possible, we also include the name
	// given to the parameter in Fortran. We exclude 'Init' from the name - where
	// the parameter name in the Fortran initialization contains 'Init'.

	

	// Model Parameters (input parameters to entomological model).
	// Please refer to Entomology.f for a more detailed description of these parameters.
	int thetap; // $\theta_p$: daysInYear
	int tau;	// $\tau$: mosqRestDuration
	int thetas;	// $\theta_s$: EIPDuration
	int n;		// $n$: nHostTypes
	int m;		// $m$: nMalHostTypes

	double Ni;		// $N_i$: popSize
	double alphai;	// $\alpha_i$: hostAvailabilityRate
	double muvA;	// $\mu_{vA}$: mosqSeekingDeathRate
	double thetad;	// $\theta_d$: mosqSeekingDuration
	double PBi;		// $P_{B_i}$: mosqProbBiting
	double PCi;		// $P_{C_i}$: mosqProbFindRestSite
	double PDi;		// $P_{D_i}$: mosqProbResting
	double PEi;		// $P_{E_i}$: mosqProbOvipositing

	// (NOTE that for this function Nv0 is an OUT parameter). //
	gsl_vector* Nv0;	// $N_{v0}$: mosqEmergeRate 
	gsl_vector* Kvi;	// $K_{vi}$: humanInfectivity

	// Derived Parameters
	// Probability that a mosquito survives one day of 
	// host-seeking but does not find a host. 
	// Dimensionless.
	// $P_A$ in model. Vector of length $\theta_p$.
	// For now, we assume that this is a double - we can change 
	// it later (no dependence on the phase of the period).
	double PA;		
	double* PAPtr;

	// Probability that on a given day, a mosquito finds a host 
	// of type $i$. 
	// Dimensionless.
	// $P_{A^i}$ in model. Matrix of size $n \times \theta_p$.
	// For now, we assume that this  is  a double: 
	// - no dependence on the phase of the period - or the
	// type of host. 
	double PAi;
	double* PAiPtr;

	// Spectral Radius of Xtp
	double srXtp;

	// Output Parameters (for the model)
	gsl_vector* Xii;	// $\Xi_i$: EIR


	// State variables.
	// $x_p(t)$: The periodic orbit of all eta state variables.
	gsl_vector** xp; 
	// $N_v^{(p)}(t)$. 
	// The periodic values of the total number of host-seeking mosquitoes.
	gsl_vector* Nvp;
	// $O_v^{(p)}(t)$. 
	// The periodic values of the number of infected host-seeking mosquitoes.
	gsl_vector* Ovp;
	// $S_v^{(p)}(t)$. 
	// The periodic values of the number of infectious host-seeking mosquitoes.
	gsl_vector* Svp;




	// Other Parameters
	// The initial estimate of the mosquito emergence rate. This is used
	// by the root finding algorithm to calculate Nv0.
	// Defined in Fortran as: MosqEmergeRateInitEstimate
	gsl_vector* Nv0guess;	

	// The number of infectious mosquitoes over every day of the cycle.
	// calculated from the EIR data.
	// $S_v$ (from EIR).
	gsl_vector* SvfromEIR;

	// The difference between SvfromEIR and SvfromEIR.
	// We don't know if we'll use this in the main simulation but for now
	// we use it to test CalcSvDiff().
	gsl_vector* SvDiff;

	double SvDiff1norm; //The $l^1$ norm of SvDiff.

	// The set of thetap matrices that determine the dynamics of the system
	// from one step to the next.
	// $\Upsilon(t)$ (over all time , $t \in [1, \theta_p]$).
	gsl_matrix** Upsilon;

	// The set of thetap vectors that determine the forcing of the system
	// from one step to the next.
	// $\Lambda(t)$ (over all time , $t \in [1, \theta_p]$).
	gsl_vector** Lambda;


	// Parameters that help to describe the order of the system.
	// Ask not why we call mt, mt. We use mt to index the system.
	// It is the maximum number of time steps we go back for $N_v$ and $O_v$.
	// 
	int mt;		
	int eta;	// $\eta$: The order of the system.
	int indexNv;	// Index of the total number of host-seeking mosquitoes.
	int indexOv;	// Index of the infected host-seeking mosquitoes.
	int indexSv;	// Index of the infectious host-seeking mosquitoes.

	// $X_{\theta_p}$.
	// The product of all the evolution matrices.
	// $X_{\theta_p} = X(\theta_p+1,1)$. 
	// Refer to Cushing (1995) and the paper for the periodic entomological model
	// for more information.
	gsl_matrix* Xtp;

	// $(\mathbb{I}-X_{\theta_p})^{-1}$.
	// The inverse of the identity matrix minus Xtp.
	gsl_matrix* inv1Xtp;


	// We define variables that we will need for the root-finding algorithm here.
	const gsl_multiroot_fsolver_type* Trootfind;
	gsl_multiroot_fsolver* srootfind;

	int status; 
	size_t iter=0;

	struct SvDiffParams pararootfind;
	gsl_multiroot_function frootfind;
	gsl_vector* xrootfind;

	// Maximum $l^1$ distance of error of root-finding algorithm
	double EpsAbsRF = 1e-6;	

	// Maximum number of iterations of root-finding algorithm.
	size_t maxiterRF = 1000;


	// Set Booleans.
	int ifRootFind = 1;

	int ifprintparameters = 0;
	int ifprintXtp = 0;
	int ifprintinv1Xtp = 0;
	int ifprintSv = 1;
	int ifprintSvDiff = 0;
	int ifprintfinalNv0 = 0;
	int ifprintfinalSvDiff = 0;
	int ifprintPO = 1;

	// File names
	// File where entomological parameters are printed out.
	char fnametestentopar[30] = "output_ento_para.txt";	
	char fnamerootfindoutput[30] = "output_rootfinding.txt";
	char xtpname[15] = "Xtp";
	char inv1Xtpname[15] = "inv1minusXtp";
	char SvfromEIRname[15] = "SvfromEIR";
	char SvDiffname[15] = "SvDifference";
	char finalNv0name[15] = "FinalNv0";
	char finalSvDiffname[15] = "FinalSvDiff";
	char Nvpname[15] = "NvPO";
	char Ovpname[15] = "OvPO";
	char Svpname[15] = "SvPO";




	/********************************************************************
	 ***********************  BEGIN CODE HERE ***************************
	 ********************************************************************/

	// Dereference pointers.
	thetap = *daysInYearPtr;
	tau = *mosqRestDurationPtr;
	thetas = *EIPDurationPtr;
	n = *nHostTypesInitPtr;
	m = *nMalHostTypesInitPtr;

	Ni = *popSizeInitPtr;
	alphai = *hostAvailabilityRateInitPtr;
	muvA = *mosqSeekingDeathRatePtr;
	thetad = *mosqSeekingDurationPtr;
	PBi = *mosqProbBitingPtr;
	PCi = *mosqProbFindRestSitePtr;
	PDi = *mosqProbRestingPtr;
	PEi = *mosqProbOvipositingPtr;


	// Set up the variables that we use to index the system.
	mt = thetas + tau -1;
	eta = 2*mt + tau;
	indexNv = 0;
	indexOv = mt;
	indexSv = 2*mt;


	// The set of thetap matrices that determine the dynamics of the system
	// from one step to the next, that is, the system is described by,
	// $x(t) = \Upsilon(t) x(t-1) = \Lambda(t)$.
	// $\Upsilon(t)$ is defined over time, $1 \leq t \leq \theta_p$, 
	// where $t \in \mathbb{N}$.
	Upsilon = (gsl_matrix**) malloc(thetap*sizeof(gsl_matrix*));

	// The set of thetap vectors that determine the forcing of the system
	// at every time step.
	// $\Lambda(t)$ is defined over time, $1 \leq t \leq \theta_p$, 
	// where $t \in \mathbb{N}$.
	Lambda = (gsl_vector**) malloc(thetap*sizeof(gsl_vector*));

	// The full periodic orbit.
	// $x_p(t)$.
	xp = (gsl_vector**) malloc(thetap*sizeof(gsl_vector*));

	// Allocate memory for gsl_vectors and initialize to 0.
	Nv0 = gsl_vector_calloc(thetap); 
	Kvi = gsl_vector_calloc(thetap);
	Xii = gsl_vector_calloc(thetap);
	Nv0guess = gsl_vector_calloc(thetap);
	SvfromEIR = gsl_vector_calloc(thetap);
	SvDiff = gsl_vector_calloc(thetap);
	Nvp = gsl_vector_calloc(thetap);
	Ovp = gsl_vector_calloc(thetap);
	Svp = gsl_vector_calloc(thetap);

	// Allocate memory for gsl_matrices and initialize to 0.
	Xtp = gsl_matrix_calloc(eta, eta);
	inv1Xtp = gsl_matrix_calloc(eta, eta);


	// Set Kvi and Xii from Fortran arrays.
	CalcCGSLVectorFromFortranArray(Kvi, FHumanInfectivityInitVector, thetap);
	// CalcCGSLVectorFromFortranArray(Xii, FEIRInitVector, thetap);
	CalcCGSLVectorFromFortranArray(SvfromEIR, FSvInitVector, thetap);
	CalcCGSLVectorFromFortranArray(Nv0guess, FMosqEmergeRateInitEstimateVector, thetap);

	// We now try to print these parameters to file to make sure that 
	// they show what we want them to show.
	if(ifprintparameters) {
		PrintParameters(fnametestentopar, thetap, tau, thetas, n, m, Ni, alphai,
			muvA, thetad, PBi, PCi, PDi, PEi, Kvi, Xii, Nv0guess);
	}
	// The parameter values look correct.




	// Initalize and reference pointers.
	PA = 0;
	PAi = 0;
	PAPtr = &PA;
	PAiPtr = &PAi;


	// Create matrices in Upsilon.
	// We also define PA and PAi in the same routine. 
	// For now, we treat PA and PAi as scalars since we are 
	// defining most parameters as scalars. If we do change things later, which we
	// may, then we will change the code accordingly. We will need to go through
	// a lot of changes anyway. 
	//
	CalcUpsilonOneHost(Upsilon, PAPtr, PAiPtr, thetap, eta, mt, tau, thetas, 
		n, m, Ni, alphai, muvA, thetad, PBi, PCi, PDi, PEi, Kvi, fnametestentopar);

	// Dereference PA and PAi from CalcUpsilon.
	PA = *PAPtr;
	PAi = *PAiPtr;

	// Calculate $X_{\theta_p}$.
	// Refer to Cushing (1995) and the paper for the periodic entomological model
	// for more information.
	FuncX(Xtp, Upsilon, thetap, 0, eta);

	if(ifprintXtp){
		PrintMatrix(fnametestentopar, xtpname, Xtp, eta, eta);
	}

	// We should now find the spectral radius of Xtp and show that it's less than 1.
	srXtp = CalcSpectralRadius(Xtp, eta, fnametestentopar);

	// printf("The spectral radius of Xtp = %e\n", srXtp);
	// getchar();

	// If the spectral radius of Xtp is greater than or equal to 1, then
	// we are not guaranteed the existence of a unique globally asymptotically
	// stable periodic orbit; thus it does not make sense to try to match the EIR
	// for this periodic orbit. 
	//
	// For this model, all the eigenvalues should be in the unit circle. However,
	// as we cannot show that analytically, we need to check it numerically.
	if(srXtp >= 1.0){
		// Throw an error.
	}

	// Calculate the inverse of (I-Xtp). 
	CalcInv1minusA(inv1Xtp, Xtp, eta, fnametestentopar);

	if(ifprintinv1Xtp){
		PrintMatrix(fnametestentopar, inv1Xtpname, inv1Xtp, eta, eta);
	}

	// Calculate the number of infectious host-seeking mosquitoes for the given EIR.
	// CalSvfromEIRdata(SvfromEIR, PAi, PBi, Ni, Xii);

	if(ifprintSv){
		PrintVector(fnametestentopar, SvfromEIRname, SvfromEIR, thetap);
	}



	
	
	

	// We can now work on the root-finding algorithm to calculate Nv0.
	// Then, we  can uncomment CalcLambda and CalcXP below.
	
	// We set a boolean for root-finding if we want to run the simulations quickly.
	
	if(ifRootFind){

		// But first we calculate what the periodic orbit would be with the initial
		// guess for Nv0 - to test calcXP	

		// We first test CalcSvDiff(). - It works.

		CalcSvDiff(SvDiff, SvfromEIR, Upsilon, Nv0guess, inv1Xtp, 
			eta, mt, thetap, fnametestentopar);
		if (ifprintSvDiff){
			PrintVector(fnametestentopar, SvDiffname, SvDiff, thetap);
		}
		SvDiff1norm = gsl_blas_dasum(SvDiff);
		printf("The $l^1$ norm of SvDiff is %.17e \n", SvDiff1norm);




		/************* We initialize variables for root-finding. **************/
		printf("Starting root-finding \n");

		// Parameters for root-finding function.
		// pararootfind = {SvfromEIR, Upsilon, inv1Xtp, eta, mt, thetap};
		pararootfind.SvfromEIR = SvfromEIR;
		pararootfind.Upsilon = Upsilon;
		pararootfind.inv1Xtp = inv1Xtp;
		pararootfind.eta = eta;
		pararootfind.mt = mt;
		pararootfind.thetap = thetap;

		// Set root-finding function.
		// frootfind = {&CalcSvDiff_rf, thetap, &pararootfind};
		frootfind.f = &CalcSvDiff_rf;
		frootfind.n = thetap;
		frootfind.params = &pararootfind;

		// Input vector for root-finding.
		xrootfind = gsl_vector_calloc(thetap);
		gsl_vector_memcpy(xrootfind, Nv0guess);

		// Set type of root-finding algorithm.
		Trootfind = gsl_multiroot_fsolver_hybrids;
		// Allocate memory for root-finding workspace.
		srootfind = gsl_multiroot_fsolver_alloc(Trootfind, thetap);

		printf("About to set root-finding solver \n");
		// Initialize root-finding.
		gsl_multiroot_fsolver_set(srootfind, &frootfind, xrootfind);
		printf("Set root-finding \n");

		// Print initial state (to screen):
		PrintRootFindingStateTS(iter, srootfind, thetap, fnamerootfindoutput);

		do{
			iter++;
			status = gsl_multiroot_fsolver_iterate(srootfind);
			PrintRootFindingStateTS(iter, srootfind, thetap, fnamerootfindoutput);

			// Check to see if solver is stuck
			if (status){
				break;
			}

			status = gsl_multiroot_test_residual(srootfind->f, EpsAbsRF);
		}
		while (status == GSL_CONTINUE && iter < maxiterRF);

		// Print status
		printf("status = %s \n", gsl_strerror(status)); 

		// Copy solution for Nv0 into Nv0.
		gsl_vector_memcpy(Nv0, srootfind->x);

		if(ifprintfinalNv0){
			PrintVector(fnametestentopar, finalNv0name, Nv0, thetap);
		}

		if(ifprintfinalSvDiff){
			PrintVector(fnametestentopar, finalSvDiffname, srootfind->f, thetap);
		}

		// Evaluate Lambda - as a pointer to an array of gsl_vectors.
		// Each vector is the forcing term at time, t.
		// The first term of each vector is the emergence rate at time, t.
		// The other terms are all zero.
		CalcLambda(Lambda, Nv0, eta, thetap, fnametestentopar);

		// We now evaluate the periodic orbit.
		// This function returns an array of thetap vectors, each of size eta.
		// Each vector, at time t, contains the state variables of the periodic
		// orbit of the system of equations.
		CalcXP(xp, Upsilon, Lambda, inv1Xtp, eta, thetap, fnametestentopar);



		// Retrieve the periodic orbits for Nv, Ov, and Sv.
		for (i=0; i<thetap; i++){
			temp = gsl_vector_get(xp[i], indexNv);
			gsl_vector_set(Nvp, i, temp);

			temp = gsl_vector_get(xp[i], indexOv);
			gsl_vector_set(Ovp, i, temp);

			temp = gsl_vector_get(xp[i], indexSv);
			gsl_vector_set(Svp, i, temp);
		}

		if(ifprintPO){
			PrintVector(fnametestentopar, Nvpname, Nvp, thetap);
			PrintVector(fnametestentopar, Ovpname, Ovp, thetap);
			PrintVector(fnametestentopar, Svpname, Svp, thetap);
		}


		//getchar();

	}
	else{
		gsl_vector_memcpy(Nv0, Nv0guess);
	}


	// Copy the mosquito emergence rate to the Fortran vector.
	CalcFortranArrayFromCGSLVector(Nv0, FMosqEmergeRateVector, thetap);

	// Deallocate memory for vectors and matrices.
	gsl_vector_free(Nv0);
	gsl_vector_free(Kvi);
	gsl_vector_free(Xii);
	gsl_vector_free(Nv0guess);
	gsl_vector_free(SvfromEIR);
	gsl_vector_free(SvDiff);
	gsl_vector_free(Nvp);
	gsl_vector_free(Ovp);
	gsl_vector_free(Svp);
	gsl_vector_free(xrootfind);

	gsl_matrix_free(Xtp);
	gsl_matrix_free(inv1Xtp);

	if(ifRootFind)
	{
		for (i=0; i<thetap; i++){
			gsl_matrix_free(Upsilon[i]);
			gsl_vector_free(Lambda[i]);
			gsl_vector_free(xp[i]);
		}

		gsl_multiroot_fsolver_free(srootfind);
	}

	free(Upsilon);
	free(Lambda);
	free(xp);


	return 0.0;
}
/********************************************************************/


/*******************************************************************/
/* CalcUpsilonOneHost returns a pointer to an array of thetap 
 * GSL matrices assuming there is only one host of humans..
 * Each matrix is Upsilon(t).
 *
 * $Upsilon(t)$ is the evolution of the mosquito population over one
 * time step. There are three main system variables:
 * $N_v$: The total number of host-seeking mosquitoes.
 * $O_v$: The number of infected host-seeking mosquitoes.
 * $S_v$: The number of infectious host-seeking mosquitoes.
 *
 * As the difference equations go back more than one time step, 
 * the size of the system is larger than 3.
 * For $N_v$ and $O_v$, we need to go back mt steps.
 * For $S_v$ we need to go back tau steps.
 * So the size of the system, eta = 2 mt + tau.
 * The first column of Upsilon(t) (indexed by 0 in C) corresponds to
 * $N_v(t)$ - as it depends on the other paramters at previous times.
 * The (mt+1)^th column of Upsilon(t) (indexed by mt in C) corresponds to
 * $O_v(t)$ - as it depends on the other paramters at previous times.
 * The (2mt+1)^th column of Upsilon(t) (indexed by 2mt in C) corresponds to
 * $S_v(t)$ - as it depends on the other paramters at previous times.
 * All other columns have 1 in the subdiagonal.
 *
 * For now, we write this code assuming that the parameters where we
 * are ignoring dependence on host type, or phase of the period, (and
 * have defined as doubles) will remove doubles. We do not code for
 * generality. If we make changes to these data types later, we will
 * change the code then. We code this, as is, to make it easier now, 
 * as we do not know what parameters we will change. It should 
 * hopefully, not be too difficult to change the code later (and
 * create a new general CalcUpsilon). Let's hope....
 *
 * 
 * Upsilon, PAPtr, and PAiPtr are OUT parameters.
 * All other parameters are IN parameters.
 */ 

void CalcUpsilonOneHost(gsl_matrix** Upsilon, double* PAPtr, 
		double* PAiPtr, int thetap, int eta, int mt, int tau, 
		int thetas, int n, int m, double Ni, double alphai, 
		double muvA, double thetad, double PBi, double PCi, double PDi, 
		double PEi, gsl_vector* Kvi, char fntestentopar[]){

	int i;
	int k;
	int l;
	// Prints intermediate results in calculating Upsilon.
	int ifPrintUpsilon = 0; 
	
	double PA;	// Described in CalcInitMosqEmergeRate.
	double PAi;	// Described in CalcInitMosqEmergeRate.
	// $P_{df}$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle.
	double Pdf; 
	// $P_{dif}$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle and gets infected.
	gsl_vector* Pdif;
	// $P_{duf}$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle and does not get infected.
	gsl_vector* Pduf; 
	
	// The probability of a mosquito surviving the extrinsic incubation
	// period. 
	// This is the sum of j from 0 to k_+ in (2.3c).
	double sumkplus; 
	double* sumkplusPtr;

	// This is an array of the sums from 0 to k_{l+} in (2.3c).
	// Note that sumklplus here is defined as sumlv in MATLAB.
	double* sumklplus;
	sumklplus = (double *)malloc((tau-1)*sizeof(double));





	double temp;
	// int itemp;


	// Initialize gsl_vectors.
	Pdif = gsl_vector_calloc(thetap);
	Pduf = gsl_vector_calloc(thetap);



	// We note again, that this code is written assuming there is only
	// one type of hosts. 
	
	// Refer to papers noted above for equations.
	PA = exp(-(alphai*Ni+muvA)*thetad);
	PAi = (1-PA)*(alphai*Ni)/(alphai*Ni+muvA);
	Pdf = PAi*PBi*PCi*PDi*PEi;

	// Evaluate Pdif and Pduf.
	// Note that these formulae are invalid for n>1.
	// We can generalize these to any n later 
	// - perhaps in a different function.
	
	// Pdif:
	gsl_vector_memcpy(Pdif, Kvi);
	gsl_vector_scale(Pdif, Pdf);

	// Pduf:
	gsl_vector_set_all(Pduf, 1.0);
	gsl_vector_sub(Pduf, Kvi);
	gsl_vector_scale(Pduf, Pdf);
	

	sumkplus = 0;
	sumkplusPtr = &sumkplus;

	// Calculate probabilities of mosquito surviving the extrinsic
	// incubation period.]
	// These currently do not depend on the phase of the period.
	CalcPSTS(sumkplusPtr, sumklplus, thetas, tau, PA, Pdf);
	sumkplus = *sumkplusPtr;

	/* // Print sumkplus and sumklplus
	printf("In CalcUpsilonOneHost, sumkplus = %f\n", sumkplus);
	for(l=1; l<=tau-1; l++){
		printf("sumklplus(%d) = %f \n", l, sumklplus[l-1]);
	}
	getchar();
	*/ 
	 



	// We start creating the matrices now.
	// Refer to Section 2.1 of JBD Paper for how this matrix is created.
	for (k=0; k < thetap; k++){
		Upsilon[k] = gsl_matrix_calloc(eta, eta);

		for (i=0; i<eta; i++){
			// Set 1's along the subdiagonal of all rows except the three
			// rows for the the main system variables.
			if(!((i==0) || (i==mt) || (i==(2*mt)))){
				gsl_matrix_set(Upsilon[k],i,i-1,1.0);
			}
		}

		
		// for $N_v$.
		gsl_matrix_set(Upsilon[k],0,0,PA);
		temp = Pdf + gsl_matrix_get(Upsilon[k], 0, tau-1);
		gsl_matrix_set(Upsilon[k],0,tau-1,temp);

		// for $O_v$.
		// We add thetap to i, to ensure that it's positive.
		// % is the mod function.
		temp = gsl_vector_get(Pdif,(k+thetap-tau)%thetap);		
		gsl_matrix_set(Upsilon[k],mt,tau-1,temp);
		gsl_matrix_set(Upsilon[k],mt,mt,PA);
		temp = gsl_vector_get(Pduf,(k+thetap-tau)%thetap) 
			   + gsl_matrix_get(Upsilon[k], mt, mt+tau-1);	
		gsl_matrix_set(Upsilon[k],mt,mt+tau-1,temp);

		// for $S_v$.
		temp = gsl_vector_get(Pdif,(k+thetap-thetas)%thetap)*sumkplus;
		gsl_matrix_set(Upsilon[k],2*mt,thetas-1,temp);
		gsl_matrix_set(Upsilon[k],2*mt,mt+thetas-1,-temp);
		for (l=1; l <= tau-1; l++){
			temp = gsl_vector_get(Pdif,(k+thetap-thetas-l)%thetap)*sumklplus[l-1];
			gsl_matrix_set(Upsilon[k],2*mt, thetas+l-1, temp);
			gsl_matrix_set(Upsilon[k],2*mt, mt+thetas+l-1, -temp);
		}
		gsl_matrix_set(Upsilon[k], 2*mt, 2*mt, PA);
		temp = Pdf + gsl_matrix_get(Upsilon[k], 2*mt, 2*mt+tau-1);
		gsl_matrix_set(Upsilon[k], 2*mt, 2*mt+tau-1, temp);


	}

	// We should try to print some of these matrices out to see what they  look like.
	if (ifPrintUpsilon){
		PrintUpsilon(fntestentopar, Upsilon, thetap, eta, PA, PAi, Pdf, Pdif, Pduf);
	}

	// Reference pointers.
	*PAPtr = PA;
	*PAiPtr = PAi;

	
	// Check to see if the binomial function works properly.
	/*
	temp = binomial(7,3);
	printf("%f \n",temp);
	getchar();
	*/ 

	// Deallocate memory for vectors
	gsl_vector_free(Pdif);
	gsl_vector_free(Pduf);

}
/********************************************************************/



/*******************************************************************/
/* CalcSvDiff_rf returns the difference between Sv for the periodic 
 * orbit for the given Nv0 and from the EIR data.
 * 
 * Given the input parameters to the entomological model, this routine
 * calculates the number of infectious host-seekign mosquitoes for the 
 * resulting periodic orbit. It then calculates the difference between 
 * this Sv and the periodic Sv calculated from the EIR data (which is 
 * the Sv from the periodic orbit of the system with the final 
 * calculated Nv0.
 *
 * This routine performs the same calculations as CalcSvDiff() but it
 * does so in the format required by the GSL multiroot-finding
 * algorithms.
 * 
 * f is an OUT parameter.
 * All other parameters are IN parameters.
 */  
int CalcSvDiff_rf(const gsl_vector* x, void* p, gsl_vector* f){

	// Add a static variable to keep track of how often we are in this routine.
	static int counterSvDiff = 0;

	// The $l^1$ norm of Svdiff.
	double SvDiff1norm;

	// Cast the incoming pointer, p, to point at a structure of type
	// SvDiffParams.
	struct SvDiffParams* params = (struct SvDiffParams *) p;

	// Assign parameters from params to variables defined in this routine.
	gsl_vector* SvfromEIR = (params->SvfromEIR);
	gsl_matrix** Upsilon = (params->Upsilon);
	gsl_matrix* inv1Xtp = (params->inv1Xtp);
	int eta = (params->eta);
	int mt = (params->mt);
	int thetap = (params->thetap);

	// It would be cleaner to read in the name of this file as an input
	// parameter but for now, we leave it out of the root-finding
	// algorithm and simply redefine it here.
	char fnametestentopar[30] = "output_ento_para.txt";	

	// Recreate a new Nv0 so that we're not restricted by const problems.
	gsl_vector* Nv0 = gsl_vector_calloc(thetap);
	gsl_vector_memcpy(Nv0, x);


	counterSvDiff++;
	//printf("In CalcSvDiff_rf for the %d th time \n", counterSvDiff);

	// To set f, we simply call CalcSvDiff. It's probably easier than rewriting
	// this code.
	CalcSvDiff(f, SvfromEIR, Upsilon, Nv0, inv1Xtp, eta, mt, thetap, fnametestentopar);

	SvDiff1norm = gsl_blas_dasum(f);
	
	// gsl_vector_set_all(f, 2.3);

	//printf("(%d) The $l^1$ norm of SvDiff is %.17e \n", counterSvDiff, SvDiff1norm);
	
	gsl_vector_free(Nv0);
	return GSL_SUCCESS;
}
/********************************************************************/







/*******************************************************************/
/* CalcSvDiff returns the difference between Sv for the periodic 
 * orbit for the given Nv0 and from the EIR data.
 * 
 * Given the input parameters to the entomological model, this routine
 * calculates the number of infectious host-seekign mosquitoes for the 
 * resulting periodic orbit. It then calculates the difference between 
 * this Sv and the periodic Sv calculated from the EIR data (which is 
 * the Sv from the periodic orbit of the system with the final 
 * calculated Nv0.
 * 
 * SvDiff is an OUT parameter.
 * All other parameters are IN parameters.
 */  
void CalcSvDiff(gsl_vector* SvDiff, gsl_vector* SvfromEIR, 
			gsl_matrix** Upsilon, gsl_vector* Nv0, gsl_matrix* inv1Xtp, 
			int eta, int mt, int thetap, char fntestentopar[]){


	int i;
	int indexSv;
	int ifprintSvfromNv0 = 0;

	char SvfromNv0name[15] = "SvfromNv0";

	double temp;

	// The set of thetap vectors that determine the forcing of the system
	// at every time step.
	// $\Lambda(t)$ is defined over time, $1 \leq t \leq \theta_p$, 
	// where $t \in \mathbb{N}$.
	gsl_vector** Lambda = (gsl_vector**) malloc(thetap*sizeof(gsl_vector*));

	// The full periodic orbit.
	// $x_p(t)$.
	gsl_vector** xp = (gsl_vector**) malloc(thetap*sizeof(gsl_vector*));

	// Periodic orbit of the number of infectious mosquitoes calculated for
	// the given Nv0.
	// $S_v$.
	gsl_vector* SvfromNv0 = gsl_vector_calloc(thetap);

	// Calculate the forcing term for each time in the period.
	CalcLambda(Lambda, Nv0, eta, thetap, fntestentopar);

	// Calculate the periodic orbit for the given Nv0.
	CalcXP(xp, Upsilon, Lambda, inv1Xtp, eta, thetap, fntestentopar);

	// Extract the number of infectious mosquitoes from the full periodic
	// orbit.
	indexSv = 2*mt;
	for (i=0; i<thetap; i++){
		temp = gsl_vector_get(xp[i], indexSv);
		gsl_vector_set(SvfromNv0, i, temp);
	}

	if(ifprintSvfromNv0){
		PrintVector(fntestentopar, SvfromNv0name, SvfromNv0, thetap);
	}

	// Subtract SvfromEIR from SvfromNv0
	gsl_vector_memcpy(SvDiff,SvfromNv0);
	gsl_vector_sub(SvDiff, SvfromEIR);


	for (i=0; i<thetap; i++){
		gsl_vector_free(Lambda[i]);
		gsl_vector_free(xp[i]);
	}

	free(Lambda);
	free(xp);
	gsl_vector_free(SvfromNv0);
}
/********************************************************************/





/*******************************************************************/
/* CalcLambda() returns a pointer to an array of thetap 
 * GSL vectors.
 * Each vector is Lambda(t).
 *
 * $Lambda(t)$ is the forcing of the mosquito population
 * at each time step, that is, it is the number of new
 * mosquitoes that enter the population at each time, $t$.
 *
 * We note here that Nv0 is a gsl_vector where the index, t, refers
 * to the mosquito emergence rate at time, t. Lambda[t] is a
 * gsl_vector that denotes the forcing at time t, where the index, i,
 * refers to the forcing to the i^th dimension of the system.
 *
 * God moves over the face of the waters,
 * Looking to the left and looking to the right,
 * But there is only water to see.
 * 
 * Lambda is an OUT parameter.
 * All other parameters are IN parameters.
 */ 
void CalcLambda(gsl_vector** Lambda, gsl_vector* Nv0, int eta,
				int thetap, char fntestentopar[]){

	int t;
	double temp;
	
	// Prints intermediate results in calculating Upsilon.
	int ifPrintLambda = 0; 
	
	
	for(t=0; t < thetap; t++){
		Lambda[t] = gsl_vector_calloc(eta);
		temp = gsl_vector_get(Nv0, t);
		gsl_vector_set(Lambda[t], 0, temp);
	}
	
	// We should try to print some of these vectors out to see what they  look like.
	if (ifPrintLambda){
		PrintLambda(Lambda, eta, fntestentopar);		
	}
}
/********************************************************************/




/*******************************************************************/
/* CalcXP returns a pointer to an array of thetap 
 * GSL vectors.
 * Each vector is is the periodic orbit solution to the main system
 * of equations at time, t.
 *
 * The size of each xp[t] is eta: the order of the system.
 *
 * This routine uses Theorem 2 of Cushing (1998) JDEA 3.
 *
 * We can probably improve the speed of this algorithm. We could replace
 * the vectors for Lambda[i] by simply using Nv0[i] and multiplying the first
 * column of the matrices, X[t,i] by Nv0[i]. 
 *
 * But for now we don't worry about speed and try to continue with the route
 * finding. There may be more work that we need to do to improve speed.
 * 
 * xp is an OUT parameter.
 * All other parameters are IN parameters.
 */  
void CalcXP(gsl_vector** xp, gsl_matrix** Upsilon, 
			gsl_vector** Lambda, gsl_matrix* inv1Xtp, int eta,
			int thetap, char fntestentopar[]){

	int t;
	int i;

	gsl_vector* vtemp = gsl_vector_calloc(eta);
	// gsl_vector* vtempsum = gsl_vector_calloc(eta);
	gsl_matrix* mtemp = gsl_matrix_calloc(eta, eta);   // scratch
    gsl_matrix* tail  = gsl_matrix_calloc(eta, eta);   // tail product
    gsl_matrix* mmul  = gsl_matrix_calloc(eta, eta);   // scratch for tail update

	// Initial condition for periodic orbit.
	gsl_vector* x0p = gsl_vector_calloc(eta);

	
	
	// Prints results in calculating the periodic orbit.
	int ifPrintx0p = 0;
	int ifPrintXP = 0; 

	char x0pname[15] = "x0p";
	

	//printf("Entered CalcXP() \n");

	/*
     * Evaluate the initial condition x0p.
     *
     * Original code:
     *   sum_{i=0..thetap-1} X(thetap, i+1) * Lambda[i]
     * where X(t,s) = Upsilon[s] * ... * Upsilon[t-1].
     *
     * This loop called FuncX() thetap times, each doing O(thetap) matrix products.
     * Replace with a backward pass maintaining a tail product:
     *   tail_i = prod_{k=i+1..thetap-1} Upsilon[k]
     * so X(thetap, i+1) == tail_i.
     */
    gsl_matrix_set_identity(tail);   // tail for i=thetap-1 is empty product = I
    gsl_vector_set_zero(vtemp);
    for (i = thetap - 1; i >= 0; --i) {
        // vtemp += tail * Lambda[i]
        gsl_blas_dgemv(CblasNoTrans, 1.0, tail, Lambda[i], 1.0, vtemp);

        // Update tail for next iteration (i-1): tail = Upsilon[i] * tail
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Upsilon[i], tail, 0.0, mmul);
        gsl_matrix_memcpy(tail, mmul);
    }
	gsl_blas_dgemv(CblasNoTrans, 1.0, inv1Xtp, vtemp, 0.0, x0p);

	// printf("Calculated initial condition for periodic orbit. \n");
	if (ifPrintx0p){
		PrintVector(fntestentopar, x0pname, x0p, eta);
	}

	// We evalute the full periodic orbit now. 
	// Note: to try to keep the indices consistent with our notes and MATLAB, 
	// xp[0] will refer to xp(1): because Upsilon[0] refers to Upsilon(1).
	// Thus, xp[thetap-1] = x0p. We can check this to make sure.

	// for(t=0; t<thetap; t++){
	// 	if(t==100 || t==200 || t==300){
	// 		printf("t=%d \n", t);
	// 	}
	// 	xp[t] = gsl_vector_calloc(eta);
	// 	// gsl_vector_set_zero(vtemp);
	// 	// gsl_vector_set_zero(vtempsum);
	// 	FuncX(mtemp, Upsilon, t+1, 0, eta);
	// 	gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, x0p, 1.0, xp[t]);
	// 	for(i=0; i<=t; i++){
	// 		// printf("t=%d i=%d \n", t, i);
	// 		FuncX(mtemp, Upsilon, t+1, i+1, eta);
	// 		gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, Lambda[i], 1.0, xp[t]);
	// 	}
	// }

	// gsl_matrix* X = gsl_matrix_calloc(eta, eta); 
	// gsl_matrix* temp = gsl_matrix_calloc(eta, eta); 
	// gsl_matrix_set_identity(X);

	xp[0] = gsl_vector_calloc(eta);
	FuncX(mtemp, Upsilon, 0, 0, eta);
	gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, x0p, 1.0, xp[0]);

	// x(t+1) = Upsilon(t)*x(t) + Lambda(t)
	for(t=1; t<thetap; t++){
		xp[t] = gsl_vector_calloc(eta);
		gsl_blas_dgemv(CblasNoTrans, 1.0, Upsilon[t-1], xp[t-1], 1.0, xp[t]);
		gsl_vector_add(xp[t], Lambda[t-1]);
	}

	//gsl_matrix_free(temp);

	// printf("Calculated periodic orbit. \n");

	// We should try to print some of these vectors out to see what they  look like.
	if (ifPrintXP){
		PrintXP(xp, eta, thetap, fntestentopar);		
	}

	/*
	t=0
		X(0,1)
		i=0
			X(1,1)
	t=1
		X(0,2)
		i=0
			X(1,2)
		i=1
			X(2,2)
	t=2
		X(0,3) 
	t
		X(t+1) = X(t) * Upsilon(t)
		xp[t] = X(t+1) * Lambda[t]

	*/


	gsl_vector_free(vtemp);
	// gsl_vector_free(vtempsum);
	gsl_vector_free(x0p);
	gsl_matrix_free(mtemp);
	gsl_matrix_free(tail);
    gsl_matrix_free(mmul);
}
/********************************************************************/



/*******************************************************************/
/* CalcPSTS() calculates probabilities of surviving the extrinsic
 * incubation period (or part of). The returned variables are the sums
 * to $k_+$ and $k_{l+}$ (including the binomial coefficients and 
 * probabilities in (2.3c) of the paper. 
 *
 * Currently, this returns scalar values because neither $P_A$, nor
 * $P_{df}$, depend on the phase of the period.
 *
 * Note that sumklplus here is defined as sumlv in MATLAB.
 * 
 * sumkplusPtr and sumklplus are OUT parameters.
 * All other parameters are IN parameter.
 */ 

void CalcPSTS(double* sumkplusPtr, double* sumklplus, int thetas,
			  int tau, double PA, double Pdf){

	int j;
	int l;
	int kplus;	// $k_+$ in model.
	int klplus;	// $k_{l+}$ in model.
	// int itemp;

	double sumkplus;
	double thetasd;
	double taud;
	double temp;
	double tempbin;
	double temppap;
	double temppdfp;


	taud = (double)tau;
	thetasd = (double) thetas;


	// klplus = (int *)malloc((tau-1)*sizeof(int)); Define temporarily.
	kplus = (int) ((thetasd/taud)-1.); // = floor(thetas/tau)-1;

	// Evaluate sumkplus
	sumkplus = 0.;
	for (j=0; j <= kplus; j++){
		tempbin = binomial(thetas-(j+1)*tau+j,j);
		temppap = pow(PA,thetas-(j+1)*tau);
		temppdfp = pow(Pdf,j);
		temp = tempbin*temppap*temppdfp;
		sumkplus=sumkplus+temp;
	}
	*sumkplusPtr = sumkplus;

	/* // Check results for kplus. 
	itemp = (int) (7./2.) - 1.;
	printf("kplus = %d\n", kplus);
	printf("itemp = %d\n", itemp);
	getchar();
	*/

	/* // Print sumkplus
	printf("In CalcPSTS(), sumkplus = %f \n", sumkplus);
	getchar();
	*/ 

	// Evaluate sumklplus
	for (l=1; l <= tau-1; l++){
		klplus = (int) (((thetasd+l)/taud) - 2); // = floor((thetas+l)/tau)-2;
		sumklplus[l-1] = 0;
		// printf("For l = %d, klplus = %d \n", l, klplus);

		for(j=0; j<=klplus; j++){
			tempbin = binomial(thetas+l-(j+2)*tau+j,j);
			temppap = pow(PA,thetas+l-(j+2)*tau);
			temppdfp = pow(Pdf,j+1);
			temp = tempbin*temppap*temppdfp;
			sumklplus[l-1] = sumklplus[l-1]+temp;
			// printf("For j = %d, tempsum = %f \n", j, temp);
		}
		// printf("sumklplus(%d) = %f \n", l, sumklplus[l-1]);
	}
}
/********************************************************************/



/*******************************************************************/
/* FuncX() calculates X(t,s).
 *
 * Note that we have to be careful with indices here. 
 * Cushing (1995) has indices starting at 0 and ending at $\theta_p -1$.
 * In our notes, and in MATLAB, the indices start at 1 and end at $\theta_p$.
 *
 *       X(t,s) = \Upsilon(t-1)*...*Upsilon(s) for t \geq s+1
 *              = I                            for t = s.
 *
 * Here, FuncX() is defined for s>=0 and t>=1.
 * 
 * X is an OUT parameter.
 * All other parameters are IN parameters.
 */ 

void FuncX(gsl_matrix* X, gsl_matrix** Upsilon, int t, int s, int eta){

	int i;

	gsl_matrix* temp = gsl_matrix_calloc(eta, eta); 

	gsl_matrix_set_identity(X);

	for (i=s; i<t; i++){
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Upsilon[i], X, 0.0, temp);
		gsl_matrix_memcpy(X, temp);
	}
	gsl_matrix_free(temp);
}
/********************************************************************/




/*******************************************************************/
/* CalcSpectralRadius() calculates the spectral radius of a given matrix.
 *
 * Given an n by n, real, nonsymmetric matrix, A, 
 * this routine calcultes its spectral radius,
 * that is, the eigenvalue with the largest absolute value.
 * 
 * A, n, and fntestentopar are IN parameters.
 */ 

double CalcSpectralRadius(gsl_matrix* A, int n, char fntestentopar[]){

	int ifprinteval = 0;	// Flag to print eigenvalues.

	int i;

	double sr;	// sprectral radius

	double temp;
	gsl_complex ztemp;

	gsl_vector* abseval = gsl_vector_calloc(n);	// Vector of the absolute values of eigenvalues.
	gsl_matrix* B = gsl_matrix_calloc(n, n); // Use to keep A safe.
	gsl_vector_complex* eval = gsl_vector_complex_calloc(n); // Vector of eigenvalues.
	// Allocate memory for workspace to evaluate the eigenvalues.
	gsl_eigen_nonsymm_workspace* w = gsl_eigen_nonsymm_alloc(n); 



	// Copy A into B to keep it safe.
	gsl_matrix_memcpy(B, A);

	// Calculate eigenvalues of B:
	gsl_eigen_nonsymm(B, eval, w);

	if (ifprinteval) {
		PrintEigenvalues(fntestentopar, eval, n);
	}

	// Calculate the absolute values of the eigenvalues.
	for(i=0; i<n; i++){
		ztemp = gsl_vector_complex_get(eval, i);
		temp = gsl_complex_abs(ztemp);
		gsl_vector_set(abseval, i, temp);
	}
	
	// Find the largest eigenvalue.
	sr = gsl_vector_max(abseval);


	// Free memory.
	gsl_matrix_free(B);
	gsl_vector_complex_free(eval);
	gsl_eigen_nonsymm_free(w);
	gsl_vector_free(abseval);

	return sr;
}
/********************************************************************/




/*******************************************************************/
/* CalcInv1minusA() calculates the inverse of (I-A) where A is a 
 * given matrix.
 *
 * Given an n by n, real matrix, A, 
 * this routine calcultes the inverse of (I-A) where I is the 
 * n by n identity matrix.
 * 
 * A, n, and fntestentopar are IN parameters.
 * inv1A is an OUT parameter.
 */ 

void CalcInv1minusA(gsl_matrix* inv1A, gsl_matrix* A, int n, char fntestentopar[]){
	
	int ifprintinv = 0;	// Flag to print inverse.

	char invname[15] = "inv1minusA"; // Name of matrix (when printing to file).

	// Data types required to compute inverse.
	gsl_matrix* B = gsl_matrix_calloc(n, n); // We calculate (I-A) in B.
	int signum;
	gsl_permutation* p = gsl_permutation_alloc(n);


	gsl_matrix_set_identity(B); // B = I.
	gsl_matrix_sub(B, A);	// B = I-A.

	// Calculate LU decomposition of (I-A).
	gsl_linalg_LU_decomp(B, p, &signum);

	// Use LU decomposition to calculate inverse.
	gsl_linalg_LU_invert(B, p, inv1A);


	if (ifprintinv) {
		PrintMatrix(fntestentopar, invname, inv1A, n, n);
	}

	

	// Free memory.
	gsl_matrix_free(B);
	gsl_permutation_free(p);


}
/********************************************************************/



/*******************************************************************/
/* CalcSvfromEIRdata() calculates Sv, given the EIR.
 *
 * Given EIR, and the parameters that determine host-biting,
 * this routine calculates the number of infectious host-seeking
 * mosquitoes, Sv.
 *
 * The EIR is assumed to be periodic so the resulting vector for Sv
 * is also periodic.
 *
 * The other parameters are constant. 
 *
 * Once the periodic paper is written, we should add a refernece to the 
 * equation that we use.
 * 
 * PAi, PBi, Ni and Xii are IN parameters.
 * Sv is an OUT parameter.
 */ 
// #include <iostream>
// void CalSvfromEIRdata(gsl_vector* Sv, double PAi, double PBi, double Ni, 
// 					  gsl_vector* Xii){

// 	double temp;
	
// 	// Sv(t) = Xii(t)*(Ni/(PAi*PBi))
// 	temp = Ni/(PAi*PBi) * 2.0; // DEBUG, remove * 2!!!!!!!!!!!
// 	gsl_vector_memcpy(Sv, Xii);
// 	gsl_vector_scale(Sv, temp);	
	
// 	std::cout << endl;
// 	std::cout << "Root-Finding: svFromEIR:" << std::endl;
// 	std::cout << gsl_vector_get(Xii, 0) << " * " << temp << " = " << gsl_vector_get(Sv, 0) << std::endl;

// }
/********************************************************************/





/*******************************************************************/
/* binomial() calculates the binomial coefficient given two integers.
 * 
 * Note that we do not check for errors.
 * 
 * All parameters are IN parameters.
 */ 

double binomial(int n, int k){
	
	unsigned int nunsigned;
	unsigned int kunsigned;
	unsigned int nminusk;

	double bc;	// Binomial coefficient

	nminusk = (unsigned int) (n-k);
	nunsigned = (unsigned int) n;
	kunsigned = (unsigned int) k;

	bc = gsl_sf_fact(nunsigned)/(gsl_sf_fact(kunsigned)*gsl_sf_fact(nminusk));

	return bc;
}
/********************************************************************/




/*******************************************************************/
/* CalcCGSLMatrixFromFortranArray() returns a GSL matrix, defined
 * according to C convention from an array passed into C from 
 * Fortran.
 * 
 * We now define a function to return a matrix, defined according to 
 * C convention from a matrix defined in Fortran, but passed as  
 * one-dimensional array. We assume that the matrix (in both C and 
 * Fortran) consists of doubles. We would need to rewrite this function
 * for other data types. We also need to ensure that the relevant
 * matrices in Fortran consists of real*8.
 *
 * We assume that the array and matrix are defined appropriately, 
 * that is, they have the correct dimensions. We do not check for
 * errors resulting from differences in sizes.
 *
 * Wanda was the name of the fish, or so thought the little boy. But
 * the fish had no name.
 *
 * CMatrix is an OUT parameter.
 * FArray is an IN parameter.
 */ 


void CalcCGSLMatrixFromFortranArray(gsl_matrix* CMatrix, double* FArray, 
				int ColLength, int RowLength){
	/* Note that ColLength is the number of rows
	 *       and RowLength is the number of columns.
	 */

	int i; // We use i to refer to the row.
	int j; // We use j to refer to the column.
	double temp; // Temporary value of i,j element.


	for (i=0; i<ColLength; i++){
		for (j=0; j<RowLength; j++){
			temp = FArray[i+j*ColLength];
			gsl_matrix_set(CMatrix, i, j, temp);
		}
	}
}
/*********************************************************************/




/*******************************************************************/
/* CalcCFortranArrayfromCGSLMatrix() returns an array defined 
 * according to Fortran matrix convention from a GSL matrix.
 *
 * This function is currently only defined for doubles. We will
 * probably need to rewrite this if we use it for anything else.
 *
 * We assume that the array and matrix are defined appropriately, 
 * that is, they have the correct dimensions. We do not check for
 * errors resulting from differences in sizes.
 *
 * FArray is an OUT parameter.
 * CMatrix is an IN parameter.
 */ 


void CalcFortranArrayFromCGSLMatrix(gsl_matrix* CMatrix, double* FArray, 
				int ColLength, int RowLength){
	/* Note that ColLength is the number of rows
	 *       and RowLength is the number of columns.
	 */
					
	int i;
	int j;
	double temp;

	for (i=0; i<ColLength; i++){
		for (j=0; j<RowLength; j++){
			temp = gsl_matrix_get(CMatrix, i, j);
			FArray[i+j*ColLength] = temp;
		}
	}
}
/********************************************************************/






/*******************************************************************/
/* CalcCGSLVectorFromFortranArray() returns a GSL vector, defined
 * according to C convention from an array passed into C from 
 * Fortran.
 * 
 * This function is currently only defined for doubles. We will
 * probably need to rewrite this if we use it for anything else.
 *
 * We assume that the array and vector are defined appropriately, 
 * that is, they have the correct dimensions. We do not check for
 * errors resulting from differences in sizes.
 *
 * CVector is an OUT parameter.
 * FArray is an IN parameter.
 */ 


void CalcCGSLVectorFromFortranArray(gsl_vector* CVector, double* FArray, 
				int Length){
	int i; 
	double temp; // Temporary value of i^{th} element.


	for (i=0; i<Length; i++){
		temp = FArray[i];
		gsl_vector_set(CVector, i, temp);
	}
}

/*********************************************************************/




/*******************************************************************/
/* CalcCFortranArrayfromCGSLVector() returns an array defined 
 * according to Fortran matrix convention from a GSL vector.
 *
 * This function is currently only defined for doubles. We will
 * probably need to rewrite this if we use it for anything else.
 *
 * We assume that the array and vector are defined appropriately, 
 * that is, they have the correct dimensions. We do not check for
 * errors resulting from differences in sizes.
 *
 * FArray is an OUT parameter.
 * CVector is an IN parameter.
 */ 


void CalcFortranArrayFromCGSLVector(gsl_vector* CVector, double* FArray, 
				int Length){
	int i; 
	double temp; // Temporary value of i^{th} element.


	for (i=0; i<Length; i++){
		temp = gsl_vector_get(CVector, i);
		FArray[i] = temp;
		
	}
}
/********************************************************************/



/********************************************************************/
/* PrintRootFindingStateTS() prints the current status of the root-
 * finding algorithm to the screen and to the given file.
 *
 * There are numerous quantities that we could print to see how the
 * root-finding algorithm is doing. It is not reasonable to print
 * all thetap terms, so for now, we print out the value of Nv0[0]
 * to see one of the values of the emergence rate, and the $l^1$
 * norm of $f$.
 *
 * Note that we print to screen and to a file.
 *
 * All parameters are IN parameters.
 */

void PrintRootFindingStateTS(size_t iter, gsl_multiroot_fsolver* srootfind, 
							 int thetap, char fnrootfindingstate[]){
    
	double svdiffsum;
	double Nv0_0;

	FILE* fpp = fopen(fnrootfindingstate, "a");

	// Calculate the $l^1$ norm of f.
	svdiffsum = gsl_blas_dasum(srootfind->f);

	// Get the 0th element of Nv0.
	Nv0_0 = gsl_vector_get(srootfind->x, 0);

	// Print to screen:
	printf("iter = %5lu Nv0(1) = % .3f ||f||_1 = % .3f \n", iter, Nv0_0, svdiffsum);
	fprintf(fpp, "iter = %5lu Nv0(1) = % .3f ||f||_1 = % .3f \n", iter, Nv0_0, svdiffsum);
	fclose(fpp);
}




/********************************************************************/
/* PrintParameters() prints the input parameters to a given file. 
 * We currently use this to make sure that the inputs we have in C
 * are what we expect from what we've sent from Fortran. 
 *
 * We may transform/copy this into a new function that does more.
 * 
 * All parameters are IN parameters.
 */

void PrintParameters(char fntestentopar[], int thetap, int tau, int thetas, 
		int n, int m, double Ni, double alphai, double muvA, 
		double thetad, double PBi, double PCi, double PDi, double PEi, 
		gsl_vector* Kvi, gsl_vector* Xii, gsl_vector* Nv0guess){

	int i;
	double temp;

	FILE* fpp = fopen(fntestentopar, "a");

	fprintf(fpp, "thetap = %d\n", thetap);
	fprintf(fpp, "tau = %d\n", tau);
	fprintf(fpp, "thetas = %d\n", thetas);
	fprintf(fpp, "n = %d\n", n);
	fprintf(fpp, "m = %d\n", m);

	fprintf(fpp, "Ni = %f\n", Ni);
	fprintf(fpp, "alphai = %f\n", alphai);
	fprintf(fpp, "muvA = %f\n", muvA);
	fprintf(fpp, "thetad = %f\n", thetad);
	fprintf(fpp, "PBi = %f\n", PBi);
	fprintf(fpp, "PCi = %f\n", PCi);
	fprintf(fpp, "PDi = %f\n", PDi);
	fprintf(fpp, "PEi = %f\n", PEi);

	
	fprintf(fpp, "Kvi = \n");
	gsl_vector_fprintf(fpp, Kvi, "%f");

	fprintf(fpp, "Xii = \n");
	gsl_vector_fprintf(fpp, Xii, "%f");

	fprintf(fpp, "Nv0guess = \n");
	gsl_vector_fprintf(fpp, Nv0guess, "%f");
	

	// Let's do this properly.
	
	for (i=0; i<thetap; i++){
		temp = gsl_vector_get(Kvi, i);
		fprintf(fpp, "Kvi(%d) = %f \n", i+1, temp);
	}

	for (i=0; i<thetap; i++){
		temp = gsl_vector_get(Xii, i);
		fprintf(fpp, "Xii(%d) = %f \n", i+1, temp);
	}

	for (i=0; i<thetap; i++){
		temp = gsl_vector_get(Nv0guess, i);
		fprintf(fpp, "Nv0guess(%d) = %f \n", i+1, temp);
	}
	

	fclose(fpp);

	

}
/********************************************************************/



/********************************************************************/
/* PrintUpsilon() prints the intermediate results while calculating 
 * Upsilon.
 * 
 * All parameters are IN parameters.
 */

void PrintUpsilon(char fntestentopar[], gsl_matrix** Upsilon, int thetap,
		int eta, double PA, double PAi, double Pdf, gsl_vector* Pdif,
		gsl_vector* Pduf){

	int i;
	int j;
	int k;
	double temp;

	FILE* fpp = fopen(fntestentopar, "a");

	fprintf(fpp, "PA = %f\n", PA);
	fprintf(fpp, "PAi = %f\n", PAi);
	fprintf(fpp, "Pdf = %f\n", Pdf);


	/*
	for (i=0; i<thetap; i++){
		temp = gsl_vector_get(Pdif, i);
		fprintf(fpp, "Pdif(%d) = %f \n", i+1, temp);
	}

	for (i=0; i<thetap; i++){
		temp = gsl_vector_get(Pduf, i);
		fprintf(fpp, "Pduf(%d) = %f \n", i+1, temp);
	}
	*/ 



	// Print some Upsilon[k].
	k=0;

	fprintf(fpp, "Upsilon[%d] = \n", k);
	for (i=0; i < eta; i++){
		for (j=0; j < eta; j++){
			temp = gsl_matrix_get(Upsilon[k], i, j);
			fprintf(fpp, "%f ", temp);
		}
		fprintf(fpp, "\n");
	}
	// gsl_matrix_fprintf(fpp, Upsilon[k], "%f");


	k = 364;

	fprintf(fpp, "Upsilon[%d] = \n", k);
	for (i=0; i < eta; i++){
		for (j=0; j < eta; j++){
			temp = gsl_matrix_get(Upsilon[k], i, j);
			fprintf(fpp, "%f ", temp);
		}
		fprintf(fpp, "\n");
	}



	fclose(fpp);

	

}
/********************************************************************/




/********************************************************************/
/* PrintXP() prints out values of XP, the periodic orbit.
 * 
 * All parameters are IN parameters.
 */
void PrintXP(gsl_vector** xp, int eta, int thetap, char fntestentopar[]){

	int t;
	char xpname[15] = "xp";
	char xpvecname[15];
	char timestring[15];
	// double temp;

	FILE* fpp = fopen(fntestentopar, "a");

	// Print all xp[t]:
	for(t=0; t<thetap; t++){
		sprintf(timestring, "%d", t+1);
		strcpy(xpvecname, xpname);
		strcat(xpvecname, "(");
		strcat(xpvecname, timestring);
		strcat(xpvecname, ")");
		PrintVector(fntestentopar, xpvecname, xp[t], eta);
	}

	fclose(fpp);
}
/********************************************************************/



/********************************************************************/
/* PrintLambda() prints some values of Lambda.
 * 
 * All parameters are IN parameters.
 */
void PrintLambda(gsl_vector** Lambda, int eta, char fntestentopar[]){

	int t;
	// double temp;

	FILE* fpp = fopen(fntestentopar, "a");

	// Print some Lambda[t].
	t=0;

	fprintf(fpp, "Lambda[%d] = \n", t);
	gsl_vector_fprintf(fpp, Lambda[t], "%f");

	t=139;

	fprintf(fpp, "Lambda[%d] = \n", t);
	gsl_vector_fprintf(fpp, Lambda[t], "%f");


	t=363;

	fprintf(fpp, "Lambda[%d] = \n", t);
	gsl_vector_fprintf(fpp, Lambda[t], "%f");


	fclose(fpp);

	

}
/********************************************************************/



/********************************************************************/
/* PrintEigenvalues() prints eigenvalues to the given file.
 * 
 * All parameters are IN parameters.
 */

void PrintEigenvalues(char fntestentopar[],gsl_vector_complex* eval, int n){

	// int i;
	// int j;
	// double temp;

	FILE* fpp = fopen(fntestentopar, "a");

	fprintf(fpp, "Eigenvalues = \n");
	/*
	for (i=0; i < eta; i++){
		for (j=0; j < eta; j++){
			temp = gsl_matrix_get(Xtp, i, j);
			fprintf(fpp, "%e ", temp);
		}
		fprintf(fpp, "\n");
	}
	*/
	gsl_vector_complex_fprintf(fpp, eval, "%e");
	fclose(fpp);
}
/********************************************************************/





/********************************************************************/
/* PrintMatrix() prints the given matrix to the given file.
 * 
 * All parameters are IN parameters.
 */

void PrintMatrix(char fntestentopar[], char matrixname[], gsl_matrix* A, 
				 int RowLength, int ColLength){

	int i;
	int j;
	double temp;

	FILE* fpp = fopen(fntestentopar, "a");

	fprintf(fpp, "%s = \n", matrixname);
	for (i=0; i < ColLength; i++){
		for (j=0; j < RowLength; j++){
			temp = gsl_matrix_get(A, i, j);
			fprintf(fpp, "%e ", temp);
		}
		fprintf(fpp, "\n");
	}

	fclose(fpp);
}
/********************************************************************/




/********************************************************************/
/* PrintVector() prints the given vector to the given file.
 * 
 * All parameters are IN parameters.
 */

void PrintVector(char fntestentopar[], char vectorname[], gsl_vector* v, int n){

	int i;
	double temp;

	FILE* fpp = fopen(fntestentopar, "a");

	for (i=0; i < n; i++){
		temp = gsl_vector_get(v, i);
		fprintf(fpp, "%s(%d) = %f \n", vectorname, i+1, temp);
		
	}
	
	fclose(fpp);
}
/********************************************************************/



