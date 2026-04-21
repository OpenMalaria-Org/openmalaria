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
#include "util/errors.h"

// Regularized smoothness solve for Nv0:
//   x = argmin ||J x - S||_2^2 + lambda * ||D2 x||_2^2
// where D2 is the *periodic* second-difference operator.
// Solves (J^T J + lambda * (D2^T D2)) x = J^T S via Cholesky (SPD).
//
// Requires GSL: <gsl/gsl_blas.h>, <gsl/gsl_linalg.h>

#include <vector>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

static inline int modp(int i, int n) {
    int r = i % n;
    return (r < 0) ? r + n : r;
}

// Adds lambda * I to A in-place.
static void add_ridge_I(gsl_matrix* A, double lambda) {
    const int n = static_cast<int>(A->size1);
    for (int i = 0; i < n; ++i) {
        gsl_matrix_set(A, i, i, gsl_matrix_get(A, i, i) + lambda);
    }
}

// Adds lambda * (D2^T D2) to A in-place, with periodic boundary conditions.
// (D2^T D2) has 5-point stencil: [ +1, -4, +6, -4, +1 ] on offsets [-2,-1,0,+1,+2].
static void add_periodic_D2tD2(gsl_matrix* A, double lambda) {
    const int n = static_cast<int>(A->size1);
    // A is assumed square n×n
    for (int i = 0; i < n; ++i) {
        const int im2 = modp(i - 2, n);
        const int im1 = modp(i - 1, n);
        const int ip1 = modp(i + 1, n);
        const int ip2 = modp(i + 2, n);

        // Diagonal
        gsl_matrix_set(A, i, i,   gsl_matrix_get(A, i, i)   + lambda * 6.0);
        // First off-diagonals
        gsl_matrix_set(A, i, im1, gsl_matrix_get(A, i, im1) + lambda * (-4.0));
        gsl_matrix_set(A, i, ip1, gsl_matrix_get(A, i, ip1) + lambda * (-4.0));
        // Second off-diagonals
        gsl_matrix_set(A, i, im2, gsl_matrix_get(A, i, im2) + lambda * 1.0);
        gsl_matrix_set(A, i, ip2, gsl_matrix_get(A, i, ip2) + lambda * 1.0);
    }
}

// Solves for x (length n):
//   (J^T J + lambdaRidge * I + lambdaD2 * D2^T D2) x = J^T S
// J is n×n (your thetap×thetap Jacobian).
// S is length n (SvfromEIR).
// x must be allocated length n.
static int solve_Nv0_smooth_reg(
    gsl_vector* x,
    const gsl_matrix* J,
    const gsl_vector* S,
    double lambdaD2,
    double lambdaRidge
) {
    const int n = static_cast<int>(J->size1);
    if (J->size2 != (size_t)n || S->size != (size_t)n || x->size != (size_t)n) {
        return GSL_EBADLEN;
    }
    if (lambdaD2 < 0.0 || lambdaRidge < 0.0) {
        return GSL_EDOM;
    }

    // A = J^T J
    gsl_matrix* A = gsl_matrix_calloc(n, n);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, J, J, 0.0, A);

    // Add ridge: lambdaRidge * I
    if (lambdaRidge > 0.0) add_ridge_I(A, lambdaRidge);

    // Add smoothness: lambdaD2 * (D2^T D2)
    if (lambdaD2 > 0.0) add_periodic_D2tD2(A, lambdaD2);

    // rhs = J^T S
    gsl_vector* rhs = gsl_vector_calloc(n);
    gsl_blas_dgemv(CblasTrans, 1.0, J, S, 0.0, rhs);

    // Solve A x = rhs
    int status = 0;

    // Try Cholesky first
    gsl_matrix* Achol = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(Achol, A);
    status = gsl_linalg_cholesky_decomp(Achol);
    if (status == 0) {
        gsl_vector_memcpy(x, rhs);
        status = gsl_linalg_cholesky_solve(Achol, rhs, x);
    } else {
        // Fallback: LU
        gsl_matrix* ALU = gsl_matrix_alloc(n, n);
        gsl_matrix_memcpy(ALU, A);
        gsl_permutation* p = gsl_permutation_alloc(n);
        int signum = 0;
        status = gsl_linalg_LU_decomp(ALU, p, &signum);
        if (status == 0) {
            status = gsl_linalg_LU_solve(ALU, p, rhs, x);
        }
        gsl_permutation_free(p);
        gsl_matrix_free(ALU);
    }

    gsl_matrix_free(Achol);
    gsl_matrix_free(A);
    gsl_vector_free(rhs);
    return status;
}

// Drop-in: choose lambdaSmooth automatically from J using SVD.
// Idea:
//   You solve (J^T J + lambda * D2^T D2) x = J^T S
//   J^T J has eigenvalues sigma_i^2 (sigma_i singular values of J).
//   D2^T D2 (periodic second-difference) has eigenvalues mu_k in [0, 16].
// We pick an "effective cutoff" sigma_c = sigma_max / kappa_target,
// and choose lambda so that lambda * mu_max ~= sigma_c^2, i.e.
//   lambda = (sigma_max / kappa_target)^2 / mu_max  with mu_max = 16.
//
// This yields a stable default without a sweep.
// Tune kappa_target (1e4..1e6 typical); larger => weaker regularization.
// Frobenius norm squared: ||J||_F^2
static double frob_norm_sq(const gsl_matrix* M)
{
    double sum = 0.0;
    for (size_t i = 0; i < M->size1; ++i) {
        for (size_t j = 0; j < M->size2; ++j) {
            const double v = gsl_matrix_get(M, i, j);
            sum += v * v;
        }
    }
    return sum;
}

// Compute lambdaSmooth = cSmooth * ||J||_F^2 / (6n)
static double compute_lambdaSmooth_from_J(
    const gsl_matrix* J,
    double cSmooth,
    bool print_diag = true)
{
    const double n = static_cast<double>(J->size1); // assume square n×n
    const double JF2 = frob_norm_sq(J);
    const double D2F2 = 6.0 * n;   // periodic second difference operator

    double lambda = cSmooth * (JF2 / D2F2);

    if (print_diag) {
        printf("Scale-based lambda selection:\n");
        printf("  ||J||_F^2 = %.6e\n", JF2);
        printf("  n = %.0f\n", n);
        printf("  ||D2||_F^2 = %.6e\n", D2F2);
        printf("  cSmooth = %.6g\n", cSmooth);
        printf("  lambdaSmooth = %.6e\n", lambda);
    }

    return lambda;
}

static double compute_lambdaRidge_from_JtJ_trace(const gsl_matrix* J, double cRidge) {
    const size_t n = J ? J->size1 : 0;
    if (!J || J->size2 != n) return 0.0;
    if (cRidge <= 0.0) return 0.0;

    // trace(J^T J) = sum_{i,j} J_ij^2 = ||J||_F^2
    double frob2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            const double v = gsl_matrix_get(J, i, j);
            frob2 += v * v;
        }
    }

    const double mean_diag_JtJ = frob2 / double(n); // (trace / n)
    return cRidge * mean_diag_JtJ;
}

struct SvDiffParams
{
	const gsl_vector* SvfromEIR;
	const std::vector<gsl_matrix*>* Upsilon;
	gsl_matrix* inv1Xtp;
	int eta;
	int mt;
	int thetap;
};

static void CalcLambda(gsl_vector** Lambda, gsl_vector* Nv0, int eta, int thetap) {
	for (int t = 0; t < thetap; ++t) {
		Lambda[t] = gsl_vector_calloc(eta);
		gsl_vector_set(Lambda[t], 0, gsl_vector_get(Nv0, t));
	}
}

[[maybe_unused]] static void CalcXP_pre_performance_patch(gsl_vector** xp,
		const std::vector<gsl_matrix*> &Upsilon, gsl_vector** Lambda,
		gsl_matrix* inv1Xtp, int eta, int thetap) {
	gsl_vector* vtemp = gsl_vector_calloc(eta);
	gsl_matrix* mtemp = gsl_matrix_calloc(eta, eta);
	gsl_vector* x0p = gsl_vector_calloc(eta);

	// Pre-7b6a2b91 version: recompute X(thetap, i+1) with FuncX() for each i.
	for (int i = 0; i < thetap; ++i) {
		FuncX(mtemp, Upsilon, thetap, i + 1, eta);
		gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, Lambda[i], 1.0, vtemp);
	}
	gsl_blas_dgemv(CblasNoTrans, 1.0, inv1Xtp, vtemp, 0.0, x0p);

	xp[0] = gsl_vector_calloc(eta);
	FuncX(mtemp, Upsilon, 0, 0, eta);
	gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, x0p, 1.0, xp[0]);

	for (int t = 1; t < thetap; ++t) {
		xp[t] = gsl_vector_calloc(eta);
		gsl_blas_dgemv(CblasNoTrans, 1.0, Upsilon[t - 1], xp[t - 1], 1.0, xp[t]);
		gsl_vector_add(xp[t], Lambda[t - 1]);
	}

	gsl_vector_free(vtemp);
	gsl_vector_free(x0p);
	gsl_matrix_free(mtemp);
}

static void CalcXP_optimized(gsl_vector** xp, const std::vector<gsl_matrix*> &Upsilon,
		gsl_vector** Lambda, gsl_matrix* inv1Xtp, int eta, int thetap) {
	gsl_vector* vtemp = gsl_vector_calloc(eta);
	gsl_matrix* mtemp = gsl_matrix_calloc(eta, eta);
	gsl_matrix* tail = gsl_matrix_calloc(eta, eta);
	gsl_matrix* mmul = gsl_matrix_calloc(eta, eta);
	gsl_vector* x0p = gsl_vector_calloc(eta);

	// Performance patch from 7b6a2b91: compute x0p with a backward tail-product
	// pass instead of calling FuncX() repeatedly inside the forcing sum.
	gsl_matrix_set_identity(tail);
	gsl_vector_set_zero(vtemp);
	for (int i = thetap - 1; i >= 0; --i) {
		gsl_blas_dgemv(CblasNoTrans, 1.0, tail, Lambda[i], 1.0, vtemp);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tail, Upsilon[i], 0.0, mmul);
		gsl_matrix_memcpy(tail, mmul);
	}
	gsl_blas_dgemv(CblasNoTrans, 1.0, inv1Xtp, vtemp, 0.0, x0p);

	xp[0] = gsl_vector_calloc(eta);
	FuncX(mtemp, Upsilon, 0, 0, eta);
	gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, x0p, 1.0, xp[0]);

	for (int t = 1; t < thetap; ++t) {
		xp[t] = gsl_vector_calloc(eta);
		gsl_blas_dgemv(CblasNoTrans, 1.0, Upsilon[t - 1], xp[t - 1], 1.0, xp[t]);
		gsl_vector_add(xp[t], Lambda[t - 1]);
	}

	gsl_vector_free(vtemp);
	gsl_vector_free(x0p);
	gsl_matrix_free(mtemp);
	gsl_matrix_free(tail);
	gsl_matrix_free(mmul);
}

static void CalcXP(gsl_vector** xp, const std::vector<gsl_matrix*> &Upsilon,
		gsl_vector** Lambda, gsl_matrix* inv1Xtp, int eta, int thetap) {
	// Uncomment the next line and comment out the optimized call below to use
	// the exact pre-performance-patch CalcXP implementation for comparison.
	//CalcXP_pre_performance_patch(xp, Upsilon, Lambda, inv1Xtp, eta, thetap);
	CalcXP_optimized(xp, Upsilon, Lambda, inv1Xtp, eta, thetap);
}

static void CalcSvDiff(gsl_vector* SvDiff, const gsl_vector* SvfromEIR,
		const std::vector<gsl_matrix*> &Upsilon, gsl_vector* Nv0, gsl_matrix* inv1Xtp,
		int eta, int mt, int thetap) {
	gsl_vector** Lambda = (gsl_vector**) malloc(thetap * sizeof(gsl_vector*));
	gsl_vector** xp = (gsl_vector**) malloc(thetap * sizeof(gsl_vector*));
	gsl_vector* SvfromNv0 = gsl_vector_calloc(thetap);

	CalcLambda(Lambda, Nv0, eta, thetap);
	CalcXP(xp, Upsilon, Lambda, inv1Xtp, eta, thetap);

	const int indexSv = 2 * mt;
	for (int i = 0; i < thetap; ++i) {
		gsl_vector_set(SvfromNv0, i, gsl_vector_get(xp[i], indexSv));
	}

	gsl_vector_memcpy(SvDiff, SvfromNv0);
	gsl_vector_sub(SvDiff, SvfromEIR);

	for (int i = 0; i < thetap; ++i) {
		gsl_vector_free(Lambda[i]);
		gsl_vector_free(xp[i]);
	}

	free(Lambda);
	free(xp);
	gsl_vector_free(SvfromNv0);
}

static int CalcSvDiff_rf(const gsl_vector* x, void* p, gsl_vector* f) {
	struct SvDiffParams* params = (struct SvDiffParams*) p;

	gsl_vector* Nv0 = gsl_vector_calloc(params->thetap);
	gsl_vector_memcpy(Nv0, x);

	CalcSvDiff(f, params->SvfromEIR, *params->Upsilon, Nv0, params->inv1Xtp,
		params->eta, params->mt, params->thetap);

	gsl_vector_free(Nv0);
	return GSL_SUCCESS;
}

[[maybe_unused]] static int solve_Nv0_rootfind(gsl_vector* Nv0, const gsl_vector* SvfromEIR,
		const std::vector<gsl_matrix*> &Upsilon, gsl_matrix* inv1Xtp,
		int eta, int mt, int thetap) {
	const double EpsAbsRF = 1e-7;
	const size_t maxiterRF = 1000;

	printf("Solving with GSL multiroot (Sv(Nv0) - SvfromEIR)... ");

	struct SvDiffParams pararootfind = {
		SvfromEIR,
		&Upsilon,
		inv1Xtp,
		eta,
		mt,
		thetap
	};

	gsl_multiroot_function frootfind;
	frootfind.f = &CalcSvDiff_rf;
	frootfind.n = thetap;
	frootfind.params = &pararootfind;

	gsl_vector* SvDiff = gsl_vector_calloc(thetap);
	CalcSvDiff(SvDiff, SvfromEIR, Upsilon, Nv0, inv1Xtp, eta, mt, thetap);
	printf("The $l^1$ norm of SvDiff is %.17e \n", gsl_blas_dasum(SvDiff));
	gsl_vector_free(SvDiff);

	gsl_vector* xrootfind = gsl_vector_alloc(thetap);
	gsl_vector_memcpy(xrootfind, Nv0);

	const gsl_multiroot_fsolver_type* Trootfind = gsl_multiroot_fsolver_hybrids;
	gsl_multiroot_fsolver* srootfind = gsl_multiroot_fsolver_alloc(Trootfind, thetap);

	printf("Starting root-finding \n");
	printf("About to set root-finding solver \n");
	int status = gsl_multiroot_fsolver_set(srootfind, &frootfind, xrootfind);
	printf("Set root-finding \n");
	size_t iter = 0;
	printf("iter = %5lu Nv0(1) = % .3f ||f||_1 = % .3f \n",
		iter, gsl_vector_get(srootfind->x, 0), gsl_blas_dasum(srootfind->f));
	if (status == GSL_SUCCESS) {
		do {
			++iter;
			status = gsl_multiroot_fsolver_iterate(srootfind);
			printf("iter = %5lu Nv0(1) = % .3f ||f||_1 = % .3f \n",
				iter, gsl_vector_get(srootfind->x, 0), gsl_blas_dasum(srootfind->f));
			if (status != GSL_SUCCESS) {
				break;
			}
			status = gsl_multiroot_test_residual(srootfind->f, EpsAbsRF);
		} while (status == GSL_CONTINUE && iter < maxiterRF);
	}

	gsl_vector_memcpy(Nv0, srootfind->x);

	if (status == GSL_CONTINUE) {
		status = GSL_EMAXITER;
	}
	if (status != GSL_SUCCESS) {
		printf("root finding failed: %s\n", gsl_strerror(status));
	} else {
		printf("Done\n");
	}

	gsl_multiroot_fsolver_free(srootfind);
	gsl_vector_free(xrootfind);
	return status;
}

[[maybe_unused]] static int solve_Nv0_linear(gsl_vector* Nv0, const gsl_vector* SvfromEIR,
		const std::vector<gsl_matrix*> &Upsilon, gsl_matrix* inv1Xtp,
		int eta, int mt, int thetap) {
	printf("Solving linear system with analytic Jacobian (Nv0 -> Sv)... ");

	gsl_matrix* J = gsl_matrix_calloc(thetap, thetap);
	CalcSvJacobian(J, Upsilon, inv1Xtp, eta, mt, thetap);

	gsl_matrix* JLU = gsl_matrix_alloc(thetap, thetap);
	gsl_matrix_memcpy(JLU, J);
	gsl_permutation* perm = gsl_permutation_alloc(thetap);
	int signum = 0;

	int status = gsl_linalg_LU_decomp(JLU, perm, &signum);
	if (status == GSL_SUCCESS) {
		status = gsl_linalg_LU_solve(JLU, perm, SvfromEIR, Nv0);
	}
	if (status) {
		printf("Linear LU Nv0 solve failed: %s\n", gsl_strerror(status));
	} else {
		printf("Done\n");
	}

	gsl_permutation_free(perm);
	gsl_matrix_free(JLU);
	gsl_matrix_free(J);
	return status;
}

[[maybe_unused]] static int solve_Nv0_linear_smooth(gsl_vector* Nv0, const gsl_vector* SvfromEIR,
		const std::vector<gsl_matrix*> &Upsilon, gsl_matrix* inv1Xtp,
		int eta, int mt, int thetap) {
	printf("Solving linear system with analytic Jacobian (Nv0 -> Sv)... ");

	gsl_matrix* J = gsl_matrix_calloc(thetap, thetap);
	CalcSvJacobian(J, Upsilon, inv1Xtp, eta, mt, thetap);

	const double cSmooth = 1e5;
	const double lambdaSmooth = compute_lambdaSmooth_from_J(J, cSmooth);

	const double cRidge = 1.5e-1;
	const double lambdaRidge = compute_lambdaRidge_from_JtJ_trace(J, cRidge);

	int status = solve_Nv0_smooth_reg(Nv0, J, SvfromEIR, lambdaSmooth, lambdaRidge);
	if (status) {
		printf("Smooth regularized Nv0 solve failed: %s\n", gsl_strerror(status));
	} else {
		printf("Done\n");
	}

	gsl_matrix_free(J);
	return status;
}

using namespace std;

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
 * mosqEmergeRateVector is an OUT parameter.
 * All other parameters are IN parameters.
 */
double CalcInitMosqEmergeRate(
	std::vector<double> &mosqEmergeRateVector, // output vector
	int thetap, // $\theta_p$: daysInYear
	int tau, // $\tau$: mosqRestDuration
	int thetas, // $\theta_s$: EIPDuration
	int n, // $n$: nHostTypes
	int m, // $m$: nMalHostTypes
	const std::vector<double> &Ni, // $N_i$: popSize (length n)
	const std::vector<double> &alphai, // $\alpha_i$: hostAvailabilityRate	(length n)
	double muvA, // $\mu_{vA}: mosqSeekingDeathRate
	double thetad, // $\theta_d$: mosqSeekingDuration
	const std::vector<double> &PBi, // $P_{B_i}$: mosqProbBiting (length n)
	const std::vector<double> &PCi, // $P_{C_i}$: mosqProbFindRestSite (length n)
	const std::vector<double> &PDi, // $P_{D_i}$: mosqProbResting (length n)
	const std::vector<double> &PEi, // $P_{E_i}$: mosqProbOvipositing -- assumed not affected by nhh
	const std::vector<double> KviInit, // Kvi (size n * thetap)
	const std::vector<double> SvInit // Sv (length n)
){
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

	gsl_vector_view Nv0_view = gsl_vector_view_array(mosqEmergeRateVector.data(), thetap);
	gsl_vector* Nv0 = &Nv0_view.vector; // $N_{v0}$: mosqEmergeRate 

	const auto Kvi_view = gsl_matrix_const_view_array(KviInit.data(), thetap, n);
	const gsl_matrix* Kvi = &Kvi_view.matrix;

	// State variables.
	// $x_p(t)$: The periodic orbit of all eta state variables.
	gsl_vector** xp = (gsl_vector**) malloc(thetap*sizeof(gsl_vector*)); 
	// $N_v^{(p)}(t)$. 
	// The periodic values of the total number of host-seeking mosquitoes.
	gsl_vector* Nvp = gsl_vector_calloc(thetap);
	// $O_v^{(p)}(t)$. 
	// The periodic values of the number of infected host-seeking mosquitoes.
	gsl_vector* Ovp = gsl_vector_calloc(thetap);
	// $S_v^{(p)}(t)$. 
	// The periodic values of the number of infectious host-seeking mosquitoes.
	gsl_vector* Svp = gsl_vector_calloc(thetap);

	// The number of infectious mosquitoes over every day of the cycle.
	// calculated from the EIR data.
	// $S_v$ (from EIR).
	gsl_vector_const_view SvfromEIR_view = gsl_vector_const_view_array(SvInit.data(), thetap);
	const gsl_vector* SvfromEIR = &SvfromEIR_view.vector;

	// The set of thetap matrices that determine the dynamics of the system
	// from one step to the next.
	// $\Upsilon(t)$ (over all time , $t \in [1, \theta_p]$).
	std::vector<gsl_matrix*> Upsilon(thetap, nullptr);

	// The set of thetap vectors that determine the forcing of the system
	// from one step to the next.
	// $\Lambda(t)$ (over all time , $t \in [1, \theta_p]$).
	gsl_vector** Lambda = (gsl_vector**) malloc(thetap*sizeof(gsl_vector*));

	// Parameters that help to describe the order of the system.
	// It is the maximum number of time steps we go back for $N_v$ and $O_v$. 
	int mt = thetas + tau -1;		
	int eta = 2*mt + tau;	// $\eta$: The order of the system.

	// $X_{\theta_p}$.
	// The product of all the evolution matrices.
	// $X_{\theta_p} = X(\theta_p+1,1)$. 
	// Refer to Cushing (1995) and the paper for the periodic entomological model
	// for more information.
	gsl_matrix* Xtp = gsl_matrix_calloc(eta, eta);

	// $(\mathbb{I}-X_{\theta_p})^{-1}$.
	// The inverse of the identity matrix minus Xtp.
	gsl_matrix* inv1Xtp = gsl_matrix_calloc(eta, eta);

	// Allocate memory for gsl_matrices and initialize to 0.
	Xtp = gsl_matrix_calloc(eta, eta);
	inv1Xtp = gsl_matrix_calloc(eta, eta);

	// Probability that a mosquito survives one day of host-seeking but does not find a host. 
	// Dimensionless.
	// $P_A$ in model. Vector of length $\theta_p$.
	double PA = 0.0;		
	
	// Probability that on a given day, a mosquito finds a host of type $i$. 
	// Dimensionless.
	// $P_{A^i}$ in model. Matrix of size $n \times \theta_p$.
	// - no dependence on the phase of the period - or the type of host. 
	double PAi = 0.0;

	// Create matrices in Upsilon.
	// We also define PA and PAi in the same routine. 
	// For now, we treat PA and PAi as scalars since we are 
	// defining most parameters as scalars. If we do change things later, which we
	// may, then we will change the code accordingly. We will need to go through
	// a lot of changes anyway. 
	CalcUpsilon(Upsilon, PA, PAi, thetap, eta, mt, tau, thetas, 
		n, m, Ni.data(), alphai.data(), muvA, thetad, PBi.data(), PCi.data(), PDi.data(), PEi.data(), Kvi);

	// Calculate $X_{\theta_p}$.
	// Refer to Cushing (1995) and the paper for the periodic entomological model
	// for more information.
	FuncX(Xtp, Upsilon, thetap, 0, eta);

	// We should now find the spectral radius of Xtp and show that it's less than 1.
	double srXtp = CalcSpectralRadius(Xtp, eta);

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

	// Uncomment exactly one solve below.
	int st = solve_Nv0_rootfind(Nv0, SvfromEIR, Upsilon, inv1Xtp, eta, mt, thetap);
	//int st = solve_Nv0_linear(Nv0, SvfromEIR, Upsilon, inv1Xtp, eta, mt, thetap);
	//int st = solve_Nv0_linear_smooth(Nv0, SvfromEIR, Upsilon, inv1Xtp, eta, mt, thetap);
	(void) st;

	bool hasNegativeNv0 = false;
	for (int i = 0; i < thetap; ++i) {
		if (gsl_vector_get(Nv0, i) < 0.0) {
			hasNegativeNv0 = true;
			break;
		}
	}

	// Deallocate memory for vectors and matrices.
	gsl_vector_free(Nvp);
	gsl_vector_free(Ovp);
	gsl_vector_free(Svp);

	gsl_matrix_free(Xtp);
	gsl_matrix_free(inv1Xtp);

	for (int i=0; i<thetap; i++)
		gsl_matrix_free(Upsilon[i]);

	free(Lambda);
	free(xp);

	if (hasNegativeNv0) {
		throw ::OM::util::base_exception(
			"Emergence rate contains negative values. Seasonality is probably too steep even after smoothing, check EIR input.",
			::OM::util::Error::VectorFitting);
	}

	return 0.0;
}
/********************************************************************/


/*******************************************************************/
/* CalcUpsilon returns a pointer to an array of thetap 
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
void CalcUpsilon(std::vector<gsl_matrix*> &Upsilon, double &PA,
		double &PAi, int thetap, int eta, int mt, int tau,
		int thetas, int n, int m, const double* Ni, const double* alphai,
		double muvA, double thetad, const double* PBi, const double* PCi, const double* PDi,
		const double* PEi, const gsl_matrix* Kvi) {

	// $P_{dif}$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle and gets infected.
	gsl_vector* Pdif = gsl_vector_calloc(thetap);

	// $P_{duf}$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle and does not get infected.
	gsl_vector* Pduf = gsl_vector_calloc(thetap);

	// PA = exp(- (sum_i alpha_i N_i + mu_vA) * theta_d)
	// (computed in the same "sum then exp" style as the original code)
	double temp = 0.0;
	for (int i = 0; i < n; i++)
		temp += alphai[i] * Ni[i];
	
	PA = exp(-(temp + muvA) * thetad);

	// Keep PAi as a scalar for minimal changes:
	// Probability of encountering one of the first m host types.
	PAi = 0.0;
	for (int i = 0; i < m; i++)
		PAi += (1.0 - PA) * (alphai[i] * Ni[i]) / (temp + muvA);

	// $P_{df}$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle.
	double Pdf = 0.0;; 

	// Pdf = sum_i PAi_i * PBi_i * PCi_i * PDi_i * PEi
	for (int i = 0; i < n; i++){
		const double PAi_i = (1.0 - PA) * (alphai[i] * Ni[i]) / (temp + muvA);
		Pdf += PAi_i * PBi[i] * PCi[i] * PDi[i] * PEi[i];
	}

	// Precompute weights w[i] (independent of k)
	std::vector<double> w(n);
	const double denom = (temp + muvA);

	for (int i = 0; i < n; ++i) {
		const double PAi_i = (1.0 - PA) * (alphai[i] * Ni[i]) / denom;
		w[i] = PAi_i * PBi[i] * PCi[i] * PDi[i] * PEi[i];
	}

	// Evaluate Pdif and Pduf.
	// Pdif(k) = sum_i [ w_i * Kvi(k,i) ]
	// Pduf(k) = Pdf - Pdif(k)
	for (int k = 0; k < thetap; ++k) {
		const double* Krow = gsl_matrix_const_ptr(Kvi, k, 0); // row k, col 0
		double pdif_k = 0.0;

		// If only human groups can have non-zero Kvi, iterate only them:
		for (int i = 0; i < m; ++i) {
			pdif_k += w[i] * Krow[i];
		}
		// If you want fully general (allow Kvi for any host), use i < n.

		if (pdif_k < 0.0) pdif_k = 0.0;
		if (pdif_k > Pdf)  pdif_k = Pdf;

		gsl_vector_set(Pdif, k, pdif_k);
		gsl_vector_set(Pduf, k, Pdf - pdif_k);
	}

	// The probability of a mosquito surviving the extrinsic incubation
	// period. 
	// This is the sum of j from 0 to k_+ in (2.3c).
	double sumkplus = 0; 

	// This is an array of the sums from 0 to k_{l+} in (2.3c).
	// Note that sumklplus here is defined as sumlv in MATLAB.
	double* sumklplus = (double *)malloc((tau-1)*sizeof(double));

	// Calculate probabilities of mosquito surviving the extrinsic
	// incubation period.]
	// These currently do not depend on the phase of the period.
	CalcPSTS(sumkplus, sumklplus, thetas, tau, PA, Pdf);

	// We start creating the matrices now.
	// Refer to Section 2.1 of JBD Paper for how this matrix is created.
	for (int k=0; k < thetap; k++){
		Upsilon[k] = gsl_matrix_calloc(eta, eta);

		for (int i=0; i<eta; i++){
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
		for (int l=1; l <= tau-1; l++){
			temp = gsl_vector_get(Pdif,(k+thetap-thetas-l)%thetap)*sumklplus[l-1];
			gsl_matrix_set(Upsilon[k],2*mt, thetas+l-1, temp);
			gsl_matrix_set(Upsilon[k],2*mt, mt+thetas+l-1, -temp);
		}
		gsl_matrix_set(Upsilon[k], 2*mt, 2*mt, PA);
		temp = Pdf + gsl_matrix_get(Upsilon[k], 2*mt, 2*mt+tau-1);
		gsl_matrix_set(Upsilon[k], 2*mt, 2*mt+tau-1, temp);
	}

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
void CalcSvJacobian(gsl_matrix* J, const std::vector<gsl_matrix*> &Upsilon, gsl_matrix* inv1Xtp,
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

		// Update tail = tail * Upsilon[s]
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tail, Upsilon[s], 0.0, mmul);
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
void CalcPSTS(double &sumkplus, double* sumklplus, int thetas,
			  int tau, double PA, double Pdf){
	int kplus = thetas / tau - 1; // $k_+$ in model.

	// Evaluate sumkplus
	sumkplus = 0.;
	for (int j=0; j <= kplus; j++){
		double tempbin = gsl_sf_choose(thetas - (j + 1) * tau + j, j); // binomial
		double temppap = pow(PA,thetas-(j+1)*tau);
		double temppdfp = pow(Pdf,j);
		double temp = tempbin*temppap*temppdfp;
		sumkplus=sumkplus+temp;
	}

	// Evaluate sumklplus
	for (int l=1; l <= tau-1; l++){
		int klplus = (thetas+l) / tau - 2; // $k_{l+}$ in model.
		sumklplus[l-1] = 0;

		for(int j=0; j<=klplus; j++){
			double tempbin = gsl_sf_choose(thetas + l - (j + 2) * tau + j, j); // binomial
			double temppap = pow(PA,thetas+l-(j+2)*tau);
			double temppdfp = pow(Pdf,j+1);
			double temp = tempbin*temppap*temppdfp;
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
void FuncX(gsl_matrix* X, const std::vector<gsl_matrix*> &Upsilon, int t, int s, int eta){
	gsl_matrix* temp = gsl_matrix_calloc(eta, eta); 

	gsl_matrix_set_identity(X);

	for (int i=s; i<t; i++){
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
	for(int i=0; i<n; i++){
		gsl_complex ztemp = gsl_vector_complex_get(eval, i);
		double temp = gsl_complex_abs(ztemp);
		gsl_vector_set(abseval, i, temp);
	}
	
	// Find the largest eigenvalue.
	double sr = gsl_vector_max(abseval); // sprectral radius

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
	gsl_permutation* p = gsl_permutation_alloc(n);

	gsl_matrix_set_identity(B); // B = I.
	gsl_matrix_sub(B, A);	// B = I-A.

	int signum;
	// Calculate LU decomposition of (I-A).
	gsl_linalg_LU_decomp(B, p, &signum);

	// Use LU decomposition to calculate inverse.
	gsl_linalg_LU_invert(B, p, inv1A);	

	// Free memory.
	gsl_matrix_free(B);
	gsl_permutation_free(p);
}
/********************************************************************/
