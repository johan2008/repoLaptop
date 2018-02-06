#ifdef WITH_GSL
	#include <gsl/gsl_math.h>
	#include <gsl/gsl_complex.h>
	#include <gsl/gsl_complex_math.h>
	#include <gsl/gsl_matrix_complex_double.h>
	#include <gsl/gsl_fft_complex.h>
	#include <gsl/gsl_fft_real.h>
	#include <gsl/gsl_blas.h>
	#include <gsl/gsl_permutation.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_fit.h>
	#include <gsl/gsl_eigen.h>
#endif

#ifdef WITH_CXSPARSE
	#include <Include/cs.h>
#endif

class LinAlgComplex;

#ifdef OS__MAC_OSX__
	//http://developer.apple.com/library/mac/#documentation/Accelerate/Reference/BLAS_Ref/Reference/reference.html
//	#include <Accelerate/Accelerate.h>
	#include <Accelerate/../Frameworks/vecLib.framework/Versions/A/Headers/cblas.h>
#else
	#include <cblas.h>
#endif
extern "C" void dgetrf_(int* n, int * n, double* a, int* lda, int * ipvt, int * info);
extern "C" void dgetri_(int* n, double* a, int* lda, int * ipvt, 
						double * work, int * lwork, int * info);
extern "C" void   dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, 
						  double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, 
						  int* lwork, int* info);
extern "C" void   dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, 
						 double *b, int *ldb, int *info);

extern "C" void dgels_(char* trans, int* m, int* n, int* nrhs, double* a, int* lda,
					   double* b, int* ldb, double* work, int* lwork, int* info );

extern "C" void dpbtrf_(char uplo, int n, int kd, double  *a,  int  lda, int *info);
extern "C" void dpbtrs_(char uplo, int n, int kd, int nrhs, double * a, int lda, 
						double * b, int ldb, int * info);
extern "C" void dpftrf_(char transfr, char uplo, int n, double  *a,  int *info);
extern "C" void dpftrs_(char transfr, char uplo, int n, int nrhs, double * A, 
						double * B, int ldb, int *info);

extern "C" void dpotrf_(char uplo, int n, double  *a,  int lda, int *info);
extern "C" void dpotrs_(char uplo, int n, int nrhs, double * A, int lda, 
						double * B, int ldb, int *info);

// NOTE: there are several lapack options for this function, see benchmark:
//			http://www.netlib.org/lapack/lug/node71.html
extern "C" void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
					   double* w, double* work, int* lwork, int* info );

extern "C" void zgesv_(int *n, int *nrhs, LinAlgComplex *a, int *lda, int *ipiv, 
					   LinAlgComplex *b, int *ldb, int *info);

extern "C" void zgels_(char* trans, int* m, int* n, int* nrhs, LinAlgComplex* a,
					   int* lda, LinAlgComplex* b, int* ldb, LinAlgComplex* work, 
					   int* lwork, int* info );
extern "C" void zgesvd_( char* jobu, char* jobvt, int* m, int* n, LinAlgComplex* a,
						int* lda, double* s, LinAlgComplex * u, int* ldu, LinAlgComplex* vt, 
						int* ldvt, LinAlgComplex* work, int* lwork, double* rwork, int* info );

extern "C" void zgetrf_(int* n, int * n, LinAlgComplex* a, int* lda, int * ipvt, int * info);
extern "C" void zgetri_(int* n, LinAlgComplex* a, int* lda, int * ipvt, 
						LinAlgComplex * work, int * lwork, int * info);



#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <string.h>

#include "VkFunctions.h"

#ifndef __LINALG_H
#define __LINALG_H

class LinAlg
{
	public:
		enum SVDMethod
		{
			SVD_AUTO,
			SVD_REGULAR, 
			SVD_GOLUB_REINSCH,		// use if M>>N
			SVD_JACOBI				// use if M>>N and higher accuracy is desired
		};
		enum InverseMethod
		{
			INV_AUTO,
			INV_LU,
			INV_SVD,
			INV_CHOLESKY
		};
};

#endif

