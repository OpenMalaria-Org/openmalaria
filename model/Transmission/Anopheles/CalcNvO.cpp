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

	const double* Ni;		// $N_i$: popSize				(length n)
	const double* alphai;	// $\alpha_i$: hostAvailabilityRate	(length n)
	double muvA;	// $\mu_{vA}$: mosqSeekingDeathRate
	double thetad;	// $\theta_d$: mosqSeekingDuration
	const double* PBi;		// $P_{B_i}$: mosqProbBiting		(length n)
	const double* PCi;		// $P_{C_i}$: mosqProbFindRestSite	(length n)
	const double* PDi;		// $P_{D_i}$: mosqProbResting		(length n)
	double PEi;		// $P_{E_i}$: mosqProbOvipositing

	// (NOTE that for this function Nv0 is an OUT parameter). //
	gsl_vector* Nv0;	// $N_{v0}$: mosqEmergeRate 
	gsl_matrix* Kvi;	// $K_{vi}$: humanInfectivity (now n x thetap)

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
	int mt;		
	int eta;	// $\eta$: The order of the system.

	// $X_{\theta_p}$.
	// The product of all the evolution matrices.
	// $X_{\theta_p} = X(\theta_p+1,1)$. 
	// Refer to Cushing (1995) and the paper for the periodic entomological model
	// for more information.
	gsl_matrix* Xtp;

	// $(\mathbb{I}-X_{\theta_p})^{-1}$.
	// The inverse of the identity matrix minus Xtp.
	gsl_matrix* inv1Xtp;
	
	int status;

	// Dereference pointers.
	thetap = *daysInYearPtr;
	tau = *mosqRestDurationPtr;
	thetas = *EIPDurationPtr;
	n = *nHostTypesInitPtr;
	m = *nMalHostTypesInitPtr;
	
	muvA = *mosqSeekingDeathRatePtr;
	thetad = *mosqSeekingDurationPtr;

	Ni = popSizeInitPtr;
	alphai = hostAvailabilityRateInitPtr;
	PBi = mosqProbBitingPtr;
	PCi = mosqProbFindRestSitePtr;
	PDi = mosqProbRestingPtr;
	PEi = *mosqProbOvipositingPtr;

	// Set up the variables that we use to index the system.
	mt = thetas + tau -1;
	eta = 2*mt + tau;

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
	Nv0guess = gsl_vector_calloc(thetap);
	SvfromEIR = gsl_vector_calloc(thetap);
	// SvDiff = gsl_vector_calloc(thetap);
	Nvp = gsl_vector_calloc(thetap);
	Ovp = gsl_vector_calloc(thetap);
	Svp = gsl_vector_calloc(thetap);

	// Allocate memory for gsl_matrices and initialize to 0.
	Kvi = gsl_matrix_calloc(n, thetap);
	Xtp = gsl_matrix_calloc(eta, eta);
	inv1Xtp = gsl_matrix_calloc(eta, eta);

	// Set Kvi, Sv and Nv0guess
	CalcCGSLMatrixFromCArray(Kvi, FHumanInfectivityInitVector, n, thetap);
	CalcCGSLVectorFromFortranArray(SvfromEIR, FSvInitVector, thetap);
	CalcCGSLVectorFromFortranArray(Nv0guess, FMosqEmergeRateInitEstimateVector, thetap);

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
	CalcUpsilon(Upsilon, PAPtr, PAiPtr, thetap, eta, mt, tau, thetas, 
		n, m, Ni, alphai, muvA, thetad, PBi, PCi, PDi, PEi, Kvi);

	// Dereference PA and PAi from CalcUpsilon.
	PA = *PAPtr;
	PAi = *PAiPtr;

	// Calculate $X_{\theta_p}$.
	// Refer to Cushing (1995) and the paper for the periodic entomological model
	// for more information.
	FuncX(Xtp, Upsilon, thetap, 0, eta);

	// We should now find the spectral radius of Xtp and show that it's less than 1.
	srXtp = CalcSpectralRadius(Xtp, eta);

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
	CalcInv1minusA(inv1Xtp, Xtp, eta);

	/************* Analytic solve: build Jacobian and solve linear system. **************/
	printf("Solving linear system with analytic Jacobian (Nv0 -> Sv)\n");

	gsl_matrix* J = gsl_matrix_calloc(thetap, thetap);
	CalcSvJacobian(J, Upsilon, inv1Xtp, eta, mt, thetap);

	// Solve J * Nv0 = SvfromEIR
	gsl_matrix* JLU = gsl_matrix_alloc(thetap, thetap);
	gsl_matrix_memcpy(JLU, J);
	gsl_permutation* perm = gsl_permutation_alloc(thetap);
	int signum = 0;
	status = gsl_linalg_LU_decomp(JLU, perm, &signum);
	if (status) {
		printf("LU_decomp failed: %s\n", gsl_strerror(status));
	}
	gsl_linalg_LU_solve(JLU, perm, SvfromEIR, Nv0);

	printf("Post-solve\n");

	gsl_permutation_free(perm);
	gsl_matrix_free(JLU);
	gsl_matrix_free(J);

	// Copy the mosquito emergence rate to the Fortran vector.
	CalcFortranArrayFromCGSLVector(Nv0, FMosqEmergeRateVector, thetap);

	// Deallocate memory for vectors and matrices.
	gsl_vector_free(Nv0);
	gsl_vector_free(Nv0guess);
	gsl_vector_free(SvfromEIR);
	gsl_vector_free(Nvp);
	gsl_vector_free(Ovp);
	gsl_vector_free(Svp);

	gsl_matrix_free(Kvi);
	gsl_matrix_free(Xtp);
	gsl_matrix_free(inv1Xtp);

	for (i=0; i<thetap; i++)
		gsl_matrix_free(Upsilon[i]);

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
void CalcUpsilon(gsl_matrix** Upsilon, double* PAPtr,
		double* PAiPtr, int thetap, int eta, int mt, int tau,
		int thetas, int n, int m, const double* Ni, const double* alphai,
		double muvA, double thetad, const double* PBi, const double* PCi, const double* PDi,
		double PEi, const gsl_matrix* Kvi){

	int i;
	int k;
	int l;
	// Prints intermediate results in calculating Upsilon.
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

	// Initialize gsl_vectors.
	Pdif = gsl_vector_calloc(thetap);
	Pduf = gsl_vector_calloc(thetap);

	// ---- Generalization to n host types (including non-human) ----

	// PA = exp(- (sum_i alpha_i N_i + mu_vA) * theta_d)
	// (computed in the same "sum then exp" style as the original code)
	temp = 0.0;
	for (i = 0; i < n; i++){
		temp += alphai[i] * Ni[i];
	}
	PA = exp(-(temp + muvA) * thetad);

	// Keep PAi as a scalar for minimal changes:
	// Probability of encountering one of the first m host types.
	PAi = 0.0;
	for (i = 0; i < m; i++){
		PAi += (1.0 - PA) * (alphai[i] * Ni[i]) / (temp + muvA);
	}

	// Pdf = sum_i PAi_i * PBi_i * PCi_i * PDi_i * PEi
	Pdf = 0.0;
	for (i = 0; i < n; i++){
		const double PAi_i = (1.0 - PA) * (alphai[i] * Ni[i]) / (temp + muvA);
		Pdf += PAi_i * PBi[i] * PCi[i] * PDi[i] * PEi;
	}

	// Evaluate Pdif and Pduf.
	// Pdif(k) = sum_i [ PAi_i * PBi_i * PCi_i * PDi_i * PEi * Kvi(i,k) ]
	// Pduf(k) = Pdf - Pdif(k)
	for (k = 0; k < thetap; k++){
		double pdif_k = 0.0;

		for (i = 0; i < n; i++){
			const double PAi_i = (1.0 - PA) * (alphai[i] * Ni[i]) / (temp + muvA);
			const double wi = PAi_i * PBi[i] * PCi[i] * PDi[i] * PEi;

			// NOTE: this file consistently uses gsl_matrix_get(M, col, row).
			// Here: col = host type i, row = time k.
			pdif_k += wi * gsl_matrix_get(Kvi, i, k);
		}

		// numerical safety
		if (pdif_k < 0.0) pdif_k = 0.0;
		if (pdif_k > Pdf) pdif_k = Pdf;

		gsl_vector_set(Pdif, k, pdif_k);
		gsl_vector_set(Pduf, k, Pdf - pdif_k);
	}

	sumkplus = 0;
	sumkplusPtr = &sumkplus;

	// Calculate probabilities of mosquito surviving the extrinsic
	// incubation period.]
	// These currently do not depend on the phase of the period.
	CalcPSTS(sumkplusPtr, sumklplus, thetas, tau, PA, Pdf);
	sumkplus = *sumkplusPtr;

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

	// Reference pointers.
	*PAPtr = PA;
	*PAiPtr = PAi;

	// Deallocate memory for vectors
	gsl_vector_free(Pdif);
	gsl_vector_free(Pduf);

	free(sumklplus);
}
/********************************************************************/


/*******************************************************************/
/* CalcSvJacobian builds the analytic Jacobian J for the map Nv0 -> Sv.
 *
 * J is a thetap x thetap matrix where:
 *   J(t,s) = d Sv(t) / d Nv0(s)
 * and Sv(t) is extracted from the periodic orbit state xp[t] at indexSv = 2*mt.
 *
 * This exploits linearity of the periodic-orbit construction:
 *  - Lambda[s] = e0 * Nv0(s)
 *  - x0 = (I - Xtp)^{-1} * sum_s X(thetap, s+1) * Lambda[s]
 *  - xp recursion: xp[t] = Upsilon[t-1] * xp[t-1] + Lambda[t-1]
 */
void CalcSvJacobian(gsl_matrix* J, gsl_matrix** Upsilon, gsl_matrix* inv1Xtp,
					int eta, int mt, int thetap){

	const int indexSv = 2 * mt;

	// D(t) is eta x thetap: D(t)[:,s] = d xp(t)[*] / d Nv0(s)
	gsl_matrix* Dcur = gsl_matrix_calloc(eta, thetap);
	gsl_matrix* Dnext = gsl_matrix_calloc(eta, thetap);

	// Scratch objects.
	gsl_matrix* tail = gsl_matrix_calloc(eta, eta);
	gsl_matrix* mmul = gsl_matrix_calloc(eta, eta);
	gsl_vector* col0 = gsl_vector_calloc(eta);
	gsl_vector* tmp = gsl_vector_calloc(eta);

	// 1) Build D(0): derivative of x0p wrt each Nv0(s).
	//    We mirror CalcXP's backward tail product so that tail at step s is
	//    exactly the X(thetap, s+1) used in the forcing sum.
	gsl_matrix_set_identity(tail);
	for (int s = thetap - 1; s >= 0; --s) {
		// col0 = tail * e0, but e0 selects the first column of tail.
		for (int r = 0; r < eta; ++r) {
			gsl_vector_set(col0, r, gsl_matrix_get(tail, r, 0));
		}

		// tmp = inv1Xtp * col0
		gsl_blas_dgemv(CblasNoTrans, 1.0, inv1Xtp, col0, 0.0, tmp);

		// Set column s of Dcur to tmp
		gsl_vector_view Dcol = gsl_matrix_column(Dcur, s);
		gsl_vector_memcpy(&Dcol.vector, tmp);

		// Update tail = Upsilon[s] * tail
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Upsilon[s], tail, 0.0, mmul);
		gsl_matrix_memcpy(tail, mmul);
	}

	// Fill Jacobian row for t=0
	for (int s = 0; s < thetap; ++s) {
		gsl_matrix_set(J, 0, s, gsl_matrix_get(Dcur, indexSv, s));
	}

	// 2) Forward sensitivity recursion matching CalcXP:
	//      xp[t] = Upsilon[t-1] * xp[t-1] + Lambda[t-1]
	//    with Lambda[t-1] = e0 * Nv0(t-1)
	for (int t = 1; t < thetap; ++t) {
		// Dnext = Upsilon[t-1] * Dcur
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Upsilon[t - 1], Dcur, 0.0, Dnext);

		// Add e0 to column (t-1): Dnext(0, t-1) += 1
		gsl_matrix_set(Dnext, 0, t - 1, gsl_matrix_get(Dnext, 0, t - 1) + 1.0);

		// Jacobian row t = row indexSv of Dnext
		for (int s = 0; s < thetap; ++s) {
			gsl_matrix_set(J, t, s, gsl_matrix_get(Dnext, indexSv, s));
		}

		// Swap Dcur and Dnext
		gsl_matrix* tmpM = Dcur;
		Dcur = Dnext;
		Dnext = tmpM;
	}

	gsl_vector_free(tmp);
	gsl_vector_free(col0);
	gsl_matrix_free(mmul);
	gsl_matrix_free(tail);
	gsl_matrix_free(Dnext);
	gsl_matrix_free(Dcur);
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

	// Evaluate sumklplus
	for (l=1; l <= tau-1; l++){
		klplus = (int) (((thetasd+l)/taud) - 2); // = floor((thetas+l)/tau)-2;
		sumklplus[l-1] = 0;

		for(j=0; j<=klplus; j++){
			tempbin = binomial(thetas+l-(j+2)*tau+j,j);
			temppap = pow(PA,thetas+l-(j+2)*tau);
			temppdfp = pow(Pdf,j+1);
			temp = tempbin*temppap*temppdfp;
			sumklplus[l-1] = sumklplus[l-1]+temp;
		}
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

double CalcSpectralRadius(gsl_matrix* A, int n){
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

void CalcInv1minusA(gsl_matrix* inv1A, gsl_matrix* A, int n){
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

	// Free memory.
	gsl_matrix_free(B);
	gsl_permutation_free(p);
}
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

void CalcCGSLMatrixFromCArray(gsl_matrix* CMatrix, double* FArray,
				int nCols, int nRows){
	int i;
	int k;
	double temp; // Temporary value

	// We assume host-major layout:
	// FArray[i*nRows + k] = element (col=i, row=k)

	for (i=0; i<nCols; i++){
		for (k=0; k<nRows; k++){
			temp = FArray[i*nRows + k];
			gsl_matrix_set(CMatrix, i, k, temp);
		}
	}
}



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


