#include <math.h>
#include <stdio.h>
#include "BasicNumerical/BasicNumerical.h"
double LinAlgMatrixComplex::s_fEpsilonForLogExp = 0.000001;
double LinAlgMatrixComplex::s_fQForLogExp = 1000;

std::map<int, std::map<int, LinAlgComplex *> > LinAlgMatrixComplex::m_LapackComplexWorkspaces;
std::map<int, std::map<int, double *> > LinAlgMatrixComplex::m_LapackWorkspaces;
std::map<int, std::map<int, int *> > LinAlgMatrixComplex::m_LapackIntWorkspaces;

double * LinAlgMatrixComplex::GetLapackWorkspace(int arraySize, int id)
{
	std::map<int, std::map<int, double *> >::iterator iter = m_LapackWorkspaces.find(id);
	if (iter == m_LapackWorkspaces.end())
	{
		double * newVal = new double[arraySize];
		m_LapackWorkspaces[id][arraySize] = newVal;
		return newVal;
	}
	else
	{
		std::map<int, double *>::iterator iter2 = iter->second.find(arraySize);
		if (iter2 == iter->second.end())
		{
			double * newVal = new double[arraySize];
			iter->second[arraySize] = newVal;
			return newVal;
		}
		return iter2->second;	
	}
	assert(false);
	return NULL;
}

LinAlgComplex * LinAlgMatrixComplex::GetLapackComplexWorkspace(int arraySize, int id)
{
	std::map<int, std::map<int, LinAlgComplex *> >::iterator iter = m_LapackComplexWorkspaces.find(id);
	if (iter == m_LapackComplexWorkspaces.end())
	{
		LinAlgComplex * newVal = new LinAlgComplex[arraySize];
		m_LapackComplexWorkspaces[id][arraySize] = newVal;
		return newVal;
	}
	else
	{
		std::map<int, LinAlgComplex *>::iterator iter2 = iter->second.find(arraySize);
		if (iter2 == iter->second.end())
		{
			LinAlgComplex * newVal = new LinAlgComplex[arraySize];
			iter->second[arraySize] = newVal;
			return newVal;
		}
		return iter2->second;	
	}
	assert(false);
	return NULL;
}

int * LinAlgMatrixComplex::GetLapackIntWorkspace(int arraySize, int id)
{
	std::map<int, std::map<int, int *> >::iterator iter = m_LapackIntWorkspaces.find(id);
	if (iter == m_LapackIntWorkspaces.end())
	{
		int * newVal = new int[arraySize];
		m_LapackIntWorkspaces[id][arraySize] = newVal;
		return newVal;
	}
	else
	{
		std::map<int, int *>::iterator iter2 = iter->second.find(arraySize);
		if (iter2 == iter->second.end())
		{
			int * newVal = new int[arraySize];
			iter->second[arraySize] = newVal;
			return newVal;
		}
		return iter2->second;	
	}
	assert(false);
	return NULL;
}

LinAlgMatrixComplex::LinAlgMatrixComplex(int m, int n)
{
	this->m = m;
	this->n = n;
	m_aComplexData = new LinAlgComplex[m*n];
}

LinAlgMatrixComplex::LinAlgMatrixComplex(const LinAlgMatrixComplex & otherMatrix)
{
	m = otherMatrix.Rows();
	n = otherMatrix.Cols();	
	m_aComplexData = new LinAlgComplex[m*n];
	
	AddMatrix(otherMatrix);
}

LinAlgMatrixComplex::~LinAlgMatrixComplex()
{
	delete [] m_aComplexData;
}

void LinAlgMatrixComplex::PrintMatrix(double absMax)
{
	for (int i=0; i<Rows(); i++)
	{
		std::cout<<"[\t";		
		for (int j=0; j<Cols(); j++)
		{
			LinAlgComplex &val = (*this)(i,j);
			if (val.Power()<absMax)
				std::cout<<"0*\t";
			else
			{
				if (val.i>0)
					std::cout<<val.r<<"+"<<val.i<<"i"<<"\t";
				else if (val.i<0)
					std::cout<<val.r<<val.i<<"i"<<"\t";
				else
					std::cout<<val.r<<"\t";
			}
				
		}
		std::cout<<"]\n";
	}
	std::cout<<std::endl;
}

void LinAlgMatrixComplex::PrintMatlab()
{
	std::cout<<"[";			
	for (int i=0; i<Rows(); i++)
	{	
		std::cout<<"[";
		for (int j=0; j<Cols(); j++)
		{
			LinAlgComplex val = (*this)(i,j);
			std::string endChar=", ";
			if (j==Cols()-1)
				endChar="";
			if (val.i>0)
				std::cout<<val.r<<"+"<<val.i<<"i"<<endChar.c_str();
			else if (val.i<0)
				std::cout<<val.r<<val.i<<"i"<<endChar.c_str();
			else
				std::cout<<val.r<<endChar.c_str();
		}
		std::cout<<"];";
	}
	std::cout<<"]"<<std::flush;	
}

double LinAlgMatrixComplex::GetMax() const
{
	double maxVal = GetValue(0, 0).Power();
	for (int i=0; i<Rows(); i++)
	for (int j=0; j<Cols(); j++)
		maxVal = vkMax(maxVal, GetValue(i, j).Power());
				
	return maxVal;
}

double LinAlgMatrixComplex::GetMin() const
{
	double minVal = GetValue(0, 0).Power();
	for (int i=0; i<Rows(); i++)
	for (int j=0; j<Cols(); j++)
		minVal = vkMin(minVal, GetValue(i,j).Power());
	
	return minVal;	
}

int LinAlgMatrixComplex::Rows() const
{
	return m;
}

int LinAlgMatrixComplex::Cols() const
{
	return n;
}
	
void LinAlgMatrixComplex::PointwiseMultiply(LinAlgMatrixComplex & m2, LinAlgMatrixComplex & trg)
{
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
		trg(i,j) = m2(i,j) * (*this)(i,j);
}

void LinAlgMatrixComplex::Multiply(LinAlgMatrixComplex & m2, LinAlgMatrixComplex & trg)
{
	//http://www.prism.gatech.edu/~ndantam3/cblas-doc/doc/html/cblas_8h.html
	assert(m2.m==n);	// inner dimensions agree
	assert(m==trg.m && m2.n==trg.n);	// outer dimensions agree with output

	LinAlgComplex alpha(1., 0);
	LinAlgComplex beta(0., 0.);
	cblas_zgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, m, m2.n, n, &alpha, 
				 m_aComplexData, m, m2.m_aComplexData, m2.m, &beta, trg.m_aComplexData, trg.m);
}

void LinAlgMatrixComplex::Multiply(LinAlgVectorComplex & v, LinAlgVectorComplex & trg)
{
	assert(n==v.m && m==trg.m);	// matrix - vector dimensions agree
	
	LinAlgComplex alpha(1., 0);
	LinAlgComplex beta(0., 0.);
	
	cblas_zgemv (CblasColMajor, CblasNoTrans, m, n, &alpha, m_aComplexData, m, 
				 v.m_aComplexData, 1, &beta, trg.m_aComplexData, 1);	
}

void LinAlgMatrixComplex::MultiplyByReal(double val)
{
	for (int i=0; i<Rows(); i++)
	for (int j=0; j<Cols(); j++)
		(*this)(i,j) *= val;
}

void LinAlgMatrixComplex::AddMatrix(const LinAlgMatrixComplex & m2)
{
	assert(Rows()==m2.Rows());
	assert(Cols()==m2.Cols());	
	for (int i=0; i<Rows(); i++)
	for (int j=0; j<Cols(); j++)
		(*this)(i,j) += m2.GetValue(i,j);
}

void LinAlgMatrixComplex::Transpose()
{
	assert(false);
	assert(Rows()==Cols());	//if matrix is not square, use Transpose(trg) instead

	for (int i=0; i<m; i++)
	for (int j=0; j<i; j++)
	{
		LinAlgComplex t = (*this)(i,j);
		(*this)(i,j) = (*this)(j,i);
		(*this)(j,i) = t;
	}
}

void LinAlgMatrixComplex::Transpose(LinAlgMatrixComplex & trg)
{
	assert(false);
	assert(m==trg.n && n==trg.m);
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
		trg(j,i)=(*this)(i,j);
}

LinAlgComplex & LinAlgMatrixComplex::operator()(int i, int j)
{
	assert(i<m && i>=0 && j<n && j>=0);
	return m_aComplexData[j*m + i];
}

LinAlgComplex LinAlgMatrixComplex::GetValue(int i, int j) const
{
	return (*((LinAlgMatrixComplex*)this))(i,j);
}

void LinAlgMatrixComplex::SetValue(int i, int j, const LinAlgComplex & value)
{
	(*this)(i,j) = value;
}


LinAlgMatrixComplex & LinAlgMatrixComplex::operator=(LinAlgMatrixComplex & src)
{
	assert(Rows()==src.Rows() && Cols()==src.Cols());
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
		(*this)(i,j) = src(i,j);
	return *this;	
}

void LinAlgMatrixComplex::Conj(LinAlgMatrixComplex & trg)
{
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
		trg(i,j) = (*this)(i,j).Conjugate();
}

	
void LinAlgMatrixComplex::Inverse()
{
	int nb=2;
	int lwork = m*n*nb;	
	assert(m==n);
	LinAlgComplex * workspace = GetLapackComplexWorkspace(lwork, 0);	
	int * pivotID = GetLapackIntWorkspace(m, 0);
	
	int info;
	int N = m;
	//std::cout<<"A=";	PrintMatlab();	std::cout<<std::endl;
	zgetrf_(&N, &N, m_aComplexData, &N, pivotID, &info);
	assert(info==0);
	zgetri_(&N, m_aComplexData, &N, pivotID, workspace, &lwork, &info );
	assert(info==0);
	//std::cout<<"Ainv=";	PrintMatlab();	std::cout<<std::endl;
	//std::cout<<"========"<<std::endl;
}


void LinAlgMatrixComplex::Solve(LinAlgVectorComplex & b, LinAlgVectorComplex &x)
{
	if (m==n)			// exact solution
		SolveLUDestructive(b, x);
	else if (m > n)		// overcontrained - least squares fit
		LeastSquaresQRDestructive(b, x);
	else
		assert(false);
}


void LinAlgMatrixComplex::SolveLUDestructive(LinAlgVectorComplex & b, LinAlgVectorComplex &xVec)
{
	// Allocate temporary memory
	assert(m==n && b.m==m && xVec.m==m);
	int *ipiv = GetLapackIntWorkspace(n, 0);
	LinAlgComplex * a = m_aComplexData;
	LinAlgComplex * x = xVec.m_aComplexData;
	for (int i = 0; i < n; i++) 
		xVec(i) = b(i);
	
	// Call LAPACK function
	int info; 
	int nrhs = 1; 
	int N = this->n;
	zgesv_(&N, &nrhs, a, &N, ipiv, x, &N, &info);
	
	// Check return status
	if (info != 0) 
	{
		fprintf(stderr, "[ERROR] Error solving system of equations: %d\n", info);
		assert(false);
	}
}


void LinAlgMatrixComplex::LeastSquaresQRDestructive(LinAlgVectorComplex & b, LinAlgVectorComplex &x)
{
	
	assert(m==b.m && n==x.m);
	
	int nrhs = 1;
	int M = m, N = n; 
	int lwork = -1;
	LinAlgComplex wkopt;
	LinAlgComplex * work;
	int info;
	
	zgels_("No transpose", &M, &N, &nrhs, m_aComplexData, &M, b.m_aComplexData, &M, &wkopt, &lwork, &info );
	lwork = (int)(wkopt.r);
	work = GetLapackComplexWorkspace(lwork, 0);
	zgels_( "No transpose", &M, &N, &nrhs, m_aComplexData, &M, b.m_aComplexData, &M, work, &lwork, &info );
	for (int i=0; i<n; i++)
		x(i) = b(i);	
}

double LinAlgMatrixComplex::Norm()
{
	int m = Rows();
	int n = Cols();	
	
	double norm=0;
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			norm += pow((*this)(i, j).Power(), 2.);
	return sqrt(norm);
}

void LinAlgMatrixComplex::Log(LinAlgMatrixComplex & X)
{
	assert(false);	//TODO: verify
	int m = Rows();
	int n = Cols();	
	
	assert(X.Rows()==Rows() && 
		   X.Cols()==Cols() &&
		   Rows()==Cols());
	int k = 0;

	LinAlgMatrixComplex A(m, n);	
	// X = A-I
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
	{
		X.SetValue(i, j, GetValue(i, j) - LinAlgComplex((i==j) ? 1. : 0.));
		A.SetValue(i, j, GetValue(i, j));
	}

	// while ||A-I|| > .5
	while (X.Norm() > .5)
//	for (int iter=0; iter < 20; iter++)
	{
		// A = sqrt(A)
		A.Sqrt(X);
		for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			A.SetValue(i, j, X.GetValue(i, j));
		k++;
		
		// X = A - I
		for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			X.SetValue(i, j, A.GetValue(i, j) - ((i==j) ? 1. : 0.));
		
//		std::cout<<"Norm="<<X.Norm()<<std::endl;
		assert (k < 50);
	}
	
	// A = I - A
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
		A.SetValue(i, j, LinAlgComplex((i==j) ? 1. : 0.) - A.GetValue(i, j));
	
	// Z = A; X = X = A; i_alg = 1
	double i_alg=1.;
	LinAlgMatrixComplex Z(m, n);
	LinAlgMatrixComplex Zsave(m, n);	
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
	{
		Z.SetValue(i, j, A.GetValue(i, j));
		X.SetValue(i, j, A.GetValue(i, j));		
	}
	
	while (Z.Norm() > s_fEpsilonForLogExp)
	{
		// Z = Z * A
		Z.Multiply(A, Zsave);
		for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			Z.SetValue(i, j, Zsave.GetValue(i, j));
		
		// X = X + Z / i
		i_alg+=1.;		
		Zsave.MultiplyByReal(1. / i_alg);
		X.AddMatrix(Zsave);
	}
	X.MultiplyByReal(pow(2., k));
	X.MultiplyByReal(-1.);	// WHY?
}

void LinAlgMatrixComplex::Exp(LinAlgMatrixComplex & X)
{
	assert(false);	//TODO: verify
	int m = Rows();
	int n = Cols();	
	
	assert(X.Rows()==Rows() && 
		   X.Cols()==Cols() &&
		   Rows()==Cols());

	// A = 2^-j A
	int j_APow = vkMax(0, 1 + floor(log2(Norm())));
	LinAlgMatrixComplex A(m,n);	
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
		A.SetValue(i, j, GetValue(i, j));
	
	A.MultiplyByReal(pow(2., -(double)j_APow));
	
	// D = I	N = I	X = I	 c = 1
	double c=1.;
	LinAlgMatrixComplex D(m,n);		
	LinAlgMatrixComplex N(m,n);		
	LinAlgMatrixComplex AX(m,n);		
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)	
	{
		D.SetValue(i, j, (i==j) ? 1. : 0.);
		N.SetValue(i, j, (i==j) ? 1. : 0.);
		X.SetValue(i, j, (i==j) ? 1. : 0.);		
	}
		
	for (int k=1; k<=s_fQForLogExp; k++)
	{
		c = c * (s_fQForLogExp - k + 1.) / (k * (2 * s_fQForLogExp - k + 1));
		
		// X = AX;	
		A.Multiply(X, AX);
		for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)	
			X.SetValue(i, j, AX.GetValue(i, j));
		
		// N = N + cX	
		AX.MultiplyByReal(c);
		N.AddMatrix(AX);
		
		// D = D + (-1)^k cX
		AX.MultiplyByReal(pow(-1., k));
		D.AddMatrix(AX);
	}
	
	// X = D^-1 N
	D.Inverse();
	D.Multiply(N, X);
	
	// X = X^(2^j)
	for (int k=0; k<j_APow; k++)
	{
		X.Multiply(X, AX);
		for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)	
			X.SetValue(i, j, AX.GetValue(i, j));		
	}
}

void LinAlgMatrixComplex::Sqrt(LinAlgMatrixComplex & X)
{
	assert(false);	//TODO: verify
	int m = Rows();
	int n = Cols();	
	assert(X.Rows()==Rows() && 
		   X.Cols()==Cols() &&
		   Rows()==Cols());

	// X=A, Y=I
	LinAlgMatrixComplex A(m,n);	
	LinAlgMatrixComplex Y(m,n);
	LinAlgMatrixComplex XX(m,n);
	LinAlgMatrixComplex iX(m,n);	
	LinAlgMatrixComplex iY(m,n);		
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
	{
		Y.SetValue(i, j, (i==j) ? 1. : 0);
		X.SetValue(i, j, GetValue(i, j));
		A.SetValue(i, j, GetValue(i, j));
	}
	
	// |XX-A|
	X.Multiply(X, XX);
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			XX.SetValue(i, j, XX.GetValue(i, j) - A.GetValue(i, j));
	
	while(XX.Norm() > s_fEpsilonForLogExp)
	{
		// iX = X^-1	iY = Y^-1
		for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
		{
			iX.SetValue(i, j, X.GetValue(i,j));
			iY.SetValue(i, j, Y.GetValue(i,j));
		}
		
		iX.Inverse();
		iY.Inverse();
		
		// X = .5 (X + iY)		
		X.AddMatrix(iY);
		X.MultiplyByReal(.5);
		
		// Y = .5 (Y+iX)
		Y.AddMatrix(iX);
		Y.MultiplyByReal(.5);
		
		// XX - A
		X.Multiply(X, XX);
		for (int i=0; i<m; i++)
			for (int j=0; j<n; j++)
				XX.SetValue(i, j, XX.GetValue(i, j) - A.GetValue(i, j));		
	}
}

void LinAlgMatrixComplex::SVDDecomposeDestructive(LinAlgMatrixComplex & U,
												  LinAlgVectorComplex & S, LinAlgMatrixComplex & Vt)
{
	int M = m;
	int N = n;

	assert(U.m==m && U.n==m);
	assert(S.m==vkMin(m, n));
	assert(Vt.m==n && Vt.n==n);
	
	int lda = m, ldu = m, ldvt = n;
	int info;
	
	LinAlgComplex wkopt;
	LinAlgComplex * work;
	int lwork=-1;
	double * rwork = GetLapackWorkspace(5*vkMin(m,n), 0);
	double * s = GetLapackWorkspace(m, 1);
	
	LinAlgComplex * a = m_aComplexData;
//	LinAlgComplex * a = GetLapackComplexWorkspace(m*n, 1);
//	for (int i=0; i<m*n; i++)
//		a[i] = m_aComplexData[i];
	
//	std::cout<<"A=";	(*this).PrintMatlab();		std::cout<<std::endl;
	
	zgesvd_( "All", "All", &M, &N, a, &lda, s, 
			U.m_aComplexData, &ldu, Vt.m_aComplexData, &ldvt, &wkopt, &lwork, rwork, &info );
	lwork = (int)wkopt.r;

	work = GetLapackComplexWorkspace(lwork, 0);
	zgesvd_("All", "All", &M, &N, a, &lda, s,
			U.m_aComplexData, &ldu, Vt.m_aComplexData, &ldvt, work, &lwork, rwork, &info );	
	
	// fill s
	for (int i=0; i<S.Rows(); i++)
		S(i) = s[i];	
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//////////////////////////// GSL STUFF //////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

#ifdef WITH_GSL
void LinAlgMatrixComplex::WriteRow(int row, gsl_complex_packed_array data)
{
	for (int i=0; i<Cols(); i++)
		gsl_matrix_complex_set(m_pData, row, i, gsl_complex_rect(data[2*i + 0], data[2*i + 1]));
}

void LinAlgMatrixComplex::WriteCol(int col, gsl_complex_packed_array data)
{
	for (int i=0; i<Rows(); i++)
		gsl_matrix_complex_set(m_pData, i, col, gsl_complex_rect(data[2*i + 0], data[2*i + 1]));	
}

void LinAlgMatrixComplex::ReadRow(int row, gsl_complex_packed_array data)
{
	for (int j=0; j<Cols(); j++)
	{
		gsl_complex val = (*this)(row, j);
		data[2*j+0] = val.dat[0];
		data[2*j+1] = val.dat[1];
	}		
}

void LinAlgMatrixComplex::ReadCol(int col, gsl_complex_packed_array data)
{
	for (int j=0; j<Rows(); j++)
	{
		gsl_complex val = (*this)(j, col);
		data[2*j+0] = val.dat[0];
		data[2*j+1] = val.dat[1];
	}	
}

void OTHERSTUFF()
{
	gsl_matrix_complex_free(m_pData);
	if (m_pWorkMatrixMxN!=NULL)
	gsl_matrix_complex_free(m_pWorkMatrixMxN);
	if (m_pWorkMatrixNxN!=NULL)
	gsl_matrix_complex_free(m_pWorkMatrixNxN);
	if (m_pWorkMatrixNxM!=NULL)
	gsl_matrix_complex_free(m_pWorkMatrixNxM);
	if (m_pWorkMatrixNxN_additional!=NULL)
	gsl_matrix_complex_free(m_pWorkMatrixNxN_additional);

	m_pWorkMatrixMxN = NULL;
	m_pWorkMatrixNxN = NULL;
	m_pWorkMatrixNxM = NULL;
	m_pWorkMatrixNxN_additional = NULL;	

	m_pData = gsl_matrix_complex_alloc (m, n);
	gsl_matrix_complex_set_zero(m_pData);
	
	//	gsl_matrix_complex_set_zero(trg.m_pData);
	//	gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gsl_complex_rect (1, 0), 
	//					m_pData, m2.m_pData, gsl_complex_rect (0, 0), trg.m_pData);	
	//	//gsl_linalg_complex_cholesky_decomp(m_pData);
	//			gsl_complex & trgVal = (*this)(i,j);
	//			trgVal.dat[0] *= val;
	//			trgVal.dat[1] *= val;			
	//gsl_matrix_complex_transpose_memcpy(trg.m_pData, m_pData);
}

//// PRIVATE
void LinAlgMatrixComplex::PrintMatlab(gsl_matrix_complex * data)
{
	std::cout<<"[";			
	for (int i=0; i<m; i++)
	{	
		std::cout<<"[";
		for (int j=0; j<n; j++)
		{
			//gsl_complex val = *(gsl_matrix_complex_ptr (data, i, j));
			LinAlgComplex val = (*this)(i,j);
			std::string endChar=", ";
			if (j==n-1)
				endChar="";
			if (val.i>0)
				std::cout<<val.r<<"+"<<val.i<<"i"<<endChar.c_str();
			else if (val.i<0)
				std::cout<<val.r<<val.i<<"i"<<endChar.c_str();
			else
				std::cout<<val.r<<endChar.c_str();
		}
		std::cout<<"];";
	}
	std::cout<<"];"<<std::endl;	
}
#endif

