#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "BasicNumerical/BasicNumerical.h"
std::map<int, std::map<int, double *> > LinAlgMatrixReal::m_LapackWorkspaces;
std::map<int, std::map<int, int *> > LinAlgMatrixReal::m_LapackIntWorkspaces;

int * LinAlgMatrixReal::GetLapackIntWorkspace(int arraySize, int id)
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

double * LinAlgMatrixReal::GetLapackWorkspace(int arraySize, int id)
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

LinAlgMatrixReal::LinAlgMatrixReal()
{
	this->m = 0;
	this->n = 0;
	m_aRealData = NULL;
	m_bFreeData = false;
}

LinAlgMatrixReal::LinAlgMatrixReal(int m, int n, double initValue)
{
	this->m = m;
	this->n = n;
	m_aRealData = new double[m * n];
	m_bFreeData = true;
	for (int i=0; i<this->m; i++)
	{
		for (int j=0; j<this->n; j++)
		{
			(*this)(i, j) = initValue;
		}
	}
}

LinAlgMatrixReal::LinAlgMatrixReal(LinAlgMatrixReal & srcMatrix)
{
	m = srcMatrix.m;
	n = srcMatrix.n;
	
	m_bFreeData = true;	
	m_aRealData = new double[m * n];
	for (int i=0; i<srcMatrix.m; i++)
		for (int j=0; j<srcMatrix.n; j++)	
			(*this)(i,j) = srcMatrix(i,j);
}

LinAlgMatrixReal::~LinAlgMatrixReal()
{
	if (m_bFreeData)
		delete [] m_aRealData;
}

void LinAlgMatrixReal::PrintMatlab()
{
	std::cout<<"["<<std::flush;
	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
			std::cout<<(*this)(i,j)<<" "<<std::flush;
		if (i!=m-1)
			std::cout<<";";
		else
			std::cout<<"]";
	}
	
}

void LinAlgMatrixReal::PrintMatrix(double absMax)
{
	for (int i=0; i<m; i++)
	{
		std::cout<<"[\t";		
		for (int j=0; j<n; j++)
		{
			if (vkAbs((*this)(i,j))<absMax)
				std::cout<<"0*\t";
			else
				std::cout<<(*this)(i,j)<<"\t";
		}
		std::cout<<"]\n";
	}
	std::cout<<std::endl;
}

double LinAlgMatrixReal::GetMax(int * saveRow, int * saveCol)
{
	double maxVal = (*this)(0,0);
	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
		{
			if (maxVal <= (*this)(i, j))
			{
				maxVal = (*this)(i, j);
				if (saveRow!=NULL)
					*saveRow = i;
				if (saveCol!=NULL)
					*saveCol = j;
			}
			
		}
	}
	return maxVal;
}

double LinAlgMatrixReal::GetMin(int * saveRow, int * saveCol)
{
	double minVal = (*this)(0, 0);
	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
		{
			if (minVal >= (*this)(i, j))
			{
				minVal = (*this)(i, j);
				if (saveRow!=NULL)
					*saveRow = i;
				if (saveCol!=NULL)
					*saveCol = j;
			}
		}	
	}
	return minVal;	
}

int LinAlgMatrixReal::Rows()
{
	return m;
}

int LinAlgMatrixReal::Cols()
{
	return n;
}


LinAlgMatrixReal & LinAlgMatrixReal::operator=(LinAlgMatrixReal & src)
{
	assert(m==src.m && n==src.n);
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			(*this)(i,j) = src(i,j);
	return *this;	
}


void LinAlgMatrixReal::Add(LinAlgMatrixReal & otherMatrix)
{
	assert(m==otherMatrix.m && n==otherMatrix.n);
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			(*this)(i,j) += otherMatrix(i,j);
}

LinAlgMatrixReal LinAlgMatrixReal::operator*(LinAlgMatrixReal & m2) 
{
	LinAlgMatrixReal C(m, m2.n);
	Multiply(m2, C);
	return C;
}

void LinAlgMatrixReal::SetValues(double value)
{
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
		(*this)(i,j) = value;
}

void LinAlgMatrixReal::PointwiseMultiply(LinAlgMatrixReal & m2, LinAlgMatrixReal & trg)
{
	assert(m==m2.m && n==m2.n && trg.m==m2.m && trg.n==m2.n);
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
		trg(i,j) = m2(i,j) * (*this)(i,j);
}


void LinAlgMatrixReal::Multiply(LinAlgMatrixReal & m2, LinAlgMatrixReal & trg)
{
	assert(trg.m==m && trg.n==m2.n && n==m2.m);
	cblas_dgemm (CblasColMajor,  CblasNoTrans, CblasNoTrans, m, m2.n, n,  
				 1.0, m_aRealData, m, m2.m_aRealData, m2.m, 0.0, trg.m_aRealData, trg.m);

}

void LinAlgMatrixReal::Multiply(LinAlgVectorReal & v, LinAlgVectorReal & vTrg)
{
	assert(v.m==n && vTrg.m==m);
	cblas_dgemv (CblasColMajor, CblasNoTrans, m, n, 1.0, m_aRealData, m, 
				 v.m_aRealData, 1, 0.0, vTrg.m_aRealData, 1);
}

void LinAlgMatrixReal::Transpose(LinAlgMatrixReal & trg)
{
	assert(m==trg.n && n==trg.m);
	for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
		trg(j,i) = (*this)(i,j);
}

void LinAlgMatrixReal::Inverse()
{
	int nb=2;
	int lwork = m*n*nb;	
	assert(m==n);
	double * workspace = GetLapackWorkspace(lwork, 0);	
	int * pivotID = GetLapackIntWorkspace(m, 0);
	
	int info;
	int N = m;
	dgetrf_(&N, &N, m_aRealData, &N, pivotID, &info);
	assert(info==0);
	dgetri_(&N, m_aRealData, &N, pivotID, workspace, &lwork, &info );
	assert(info==0);
}


void LinAlgMatrixReal::Solve(LinAlgVectorReal & b, LinAlgVectorReal &x)
{
	if (m==n)			// exact solution
		SolveLUDestructive(b, x);
	else if (m > n)		// overcontrained - least squares fit
		LeastSquaresQRDestructive(b, x);
	else
		assert(false);
}

void LinAlgMatrixReal::SVDDecomposeDestructive(LinAlgMatrixReal & U, LinAlgVectorReal & S, 
											   LinAlgMatrixReal & Vt)
{
	// example: http://software.intel.com/sites/products/documentation/hpc/mkl/lapack/mkl_lapack_examples/dgesvd_ex.c.htm
	int info, lwork;
	
	assert(U.Rows()==m);		assert(U.Cols()==m);
	assert(Vt.Rows()==n);		assert(Vt.Cols()==n);
	assert(S.Rows()==vkMin(m, n));
	
	double * a = m_aRealData;
	double * s = S.m_aRealData;
	double * u = U.m_aRealData;
	double * vt = Vt.m_aRealData;
		
	double wkopt;
	double* work;
	
	lwork = -1;
	dgesvd_( "All", "All", &m, &n, a, &m, s, u, &m, vt, &n, &wkopt, &lwork, &info);
	lwork = (int)wkopt;
	work = GetLapackWorkspace(lwork, 4);
	dgesvd_( "All", "All", &m, &n, a, &m, s, u, &m, vt, &n, work, &lwork, &info);	
	assert(info==0);	// check that algorithm converged
}


void LinAlgMatrixReal::SolveLUDestructive(LinAlgVectorReal & b, LinAlgVectorReal &xVec)
{
	// Allocate temporary memory
	assert(m==n && b.m==m && xVec.m==m);
	int *ipiv = GetLapackIntWorkspace(n, 0);
	double * a = m_aRealData;
	double * x = xVec.m_aRealData;
	for (int i = 0; i < n; i++) 
		xVec(i) = b(i);
	    
	// Call LAPACK function
	int info; 
	int nrhs = 1; 
	int n = this->n;
	dgesv_(&n, &nrhs, a, &n, ipiv, x, &n, &info);
	
	// Check return status
	if (info != 0) 
	{
		fprintf(stderr, "[ERROR] Error solving system of equations: %d\n", info);
		assert(false);
	}
}


void LinAlgMatrixReal::LeastSquaresQRDestructive(LinAlgVectorReal & b, LinAlgVectorReal &x)
{
	assert(m>n);	// if underconstrained - can use same as below, but store output in x
	assert(m==b.m && n==x.m);
	
	int nrhs = 1;
	int M = m, N = n; 
	int lwork = -1;
	double wkopt;
	double * work;
	int info;
	
	dgels_("No transpose", &M, &N, &nrhs, m_aRealData, &M, b.m_aRealData, &M, &wkopt, &lwork, &info );
	lwork = (int)(wkopt);
	work = GetLapackWorkspace(lwork, 0);
	dgels_( "No transpose", &M, &N, &nrhs, m_aRealData, &M, b.m_aRealData, &M, work, &lwork, &info );
	for (int i=0; i<n; i++)
		x(i) = b(i);
}

void LinAlgMatrixReal::DecomposeCholesky()
{
	assert(false);	// cannot fix - segfault
	assert(m==n);
	std::cout<<"DecomposeCholesky: m="<<m<<" x "<<n<<" dpotrf"<<std::endl;
	int info=0;
//	int kd=n-1;
//	dpbtrf('L', n, kd, m_aRealData,  n, &info);
//	dpftrf('N', 'L', n, m_aRealData, &info);
//	dpotrf('U', n, m_aRealData,  n, &info);
	std::cout<<"DecomposeCholesky: DONE!"<<std::endl;	
	assert(info==0);
}

void LinAlgMatrixReal::SolveDecomposedCholesky(LinAlgVectorReal & b, LinAlgVectorReal &x)
{
	assert(false);	// cannot decompose
	assert(m==n);
	for (int i = 0; i < n; i++) 
		x(i) = b(i);

	int info=0;
//	dpbtrs('U', n, n-1, 1, m_aRealData, n, x.m_aRealData, 1, &info);
//	dpftrs('N', 'L', n, 1, m_aRealData, x.m_aRealData, 1, &info);
	dpotrs_('L', n, 1, m_aRealData,  n, x.m_aRealData, 1, &info);
	assert(info==0);
}


double LinAlgMatrixReal::Norm(double p)
{
	double sum = 0;
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			sum += pow((*this)(i,j), p);
	return sum;
}

void LinAlgMatrixReal::MultiplyBy(double value)
{
	assert(value!=0);
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			(*this)(i,j) *= value;
}

void LinAlgMatrixReal::AddVal(int i, int j, double val)
{
	(*this)(i,j) += val;
}

double & LinAlgMatrixReal::operator()(int i, int j)
{
	assert(i<m && i>=0 && j<n && j>=0);
	return m_aRealData[j*m+i];
}

bool LinAlgMatrixReal::WriteMatrixBin(const char * outputFile)
{	
	FILE * file = fopen(outputFile, "wb");
	assert (file!=NULL);
	fwrite(m_aRealData, sizeof(double), m*n, file);
	fclose(file);
	return true;
}

bool LinAlgMatrixReal::ReadMatrixBin(const char * inputFile)
{
	FILE * file = fopen(inputFile, "rb");
	if (file==NULL)
		return false;
	fread(m_aRealData, sizeof(double), m*n, file);
	fclose(file);
	return true;
}

bool LinAlgMatrixReal::WriteMatrixASCII(const char * outputFile)
{
	std::ofstream textStream(outputFile);
	assert(textStream.is_open());
	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
		{
			if (j<n-1)
				textStream<<(*this)(i,j)<<" ";
			else
				textStream<<(*this)(i,j)<<"\n";
		}
	}
	return true;
}

bool LinAlgMatrixReal::ReadMatrixASCII(const char * inputFile)
{
	std::ifstream textStream(inputFile);
	if (!textStream.is_open())
		return false;
	
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			textStream>>(*this)(i,j);
	return true;
}

///////////// Eigen decomposition ///////////
void LinAlgMatrixReal::EigendecompositionSymmetric(LinAlgVectorReal & eigenValues, 
												   LinAlgMatrixReal & eigenVectors)
{	
	assert(m==n && eigenValues.m==m && eigenVectors.m==m && eigenVectors.n==n);
	int N = n;
	int lwork, info;
	double wkopt;
	double* work;
	
	eigenVectors = *this;
	
	double * a = eigenVectors.m_aRealData;
	double * w = eigenValues.m_aRealData;
	lwork = -1;
	dsyev_( "Vectors", "Upper", &N, a, &N, w, &wkopt, &lwork, &info );
	lwork = (int)wkopt;
	work = GetLapackWorkspace(lwork, 0);
	dsyev_( "Vectors", "Upper", &N, a, &N, w, work, &lwork, &info );
	if (info!=0)
	{
		std::cout<<"[ERROR] LAPACK: dsyevd, info="<<info<<std::endl;
		assert(info==0);
	}
	
	// sort eigenvalues by absolute value
	int * mapIDToPosition = GetLapackIntWorkspace(N, 0);
	int * mapPositionToID = GetLapackIntWorkspace(N, 1);
	std::vector<Sortable> eigs;
	for (int i=0; i<N; i++)
	{
		eigs.push_back(Sortable(vkAbs(eigenValues(i)), NULL, i));
		mapIDToPosition[i] = i;
		mapPositionToID[i] = i;
	}
	std::sort(eigs.begin(), eigs.end());
	for (int i=0; i<N; i++)
	{
		int putHereID = eigs[N-i-1].id;
		int putHerePosition = mapIDToPosition[putHereID];
	
		// swap values		
		double temp = eigenValues(putHerePosition);
		eigenValues(putHerePosition) = eigenValues(i);
		eigenValues(i) = temp;
		
		// swap columns
		for (int j=0; j<N; j++)
		{
			temp = eigenVectors(j, putHerePosition);
			eigenVectors(j, putHerePosition) = eigenVectors(j, i);
			eigenVectors(j, i) = temp;
		}
		
		mapIDToPosition[putHereID] = i;
		mapIDToPosition[mapPositionToID[i]] = putHerePosition;
		mapPositionToID[putHerePosition] = mapPositionToID[i];
		mapPositionToID[i] = putHereID;
	}
}

void LinAlgMatrixReal::RowNormalize()
{
	for (int i=0; i<m; i++)
	{
		double norm = 0;
		for (int j=0; j<n; j++)
			norm += (*this)(i,j);
		if (norm!=0)
			for (int j=0; j<n; j++)
				(*this)(i,j)/=norm;
	}
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//////////////////////////// GSL STUFF //////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



#ifdef WITH_GSL
// ALL GSL STUFF
void tempstuff
{
	//	m_pData = gsl_matrix_alloc (m, n);
	//	gsl_matrix_set_all(m_pData, initValue);
		//gsl_matrix_free (m_pData);	
	//double real = gsl_matrix_complex_get(src.m_pData, i, j).dat[0];
	//double imaginary = gsl_matrix_complex_get(src.m_pData, i, j).dat[1];
	//gsl_matrix_set(m_pData, i, j, sqrt(real*real+imaginary*imaginary));
	

}

// ///////// STATIC WORKSPACE //////////
std::vector<std::map<int, gsl_vector *>  > LinAlgMatrixReal::s_WorspaceVector;
std::vector<std::map<int, std::map<int, gsl_matrix*> > > LinAlgMatrixReal::s_WorspaceMatrix;	
std::map<int, gsl_eigen_symmv_workspace *> LinAlgMatrixReal::m_EigenNonSymmWorkspace;

gsl_vector * LinAlgMatrixReal::GetWorkspaceVector(int N, int unique)
{
	while (unique>=(int)s_WorspaceVector.size())
		s_WorspaceVector.push_back(std::map<int, gsl_vector *>());
		
	if (s_WorspaceVector[unique].find(N)!=s_WorspaceVector[unique].end())
		return s_WorspaceVector[unique][N];
	
	s_WorspaceVector[unique][N] = gsl_vector_alloc(N);
	return s_WorspaceVector[unique][N];
}

gsl_matrix * LinAlgMatrixReal::GetWorkspaceMatrix(int N, int M, int unique)
{
	while (unique>=(int)s_WorspaceMatrix.size())
		s_WorspaceMatrix.push_back(std::map<int, std::map<int, gsl_matrix *> >());
	
	if (s_WorspaceMatrix[unique].find(N)!=s_WorspaceMatrix[unique].end() 
		&& s_WorspaceMatrix[unique][N].find(M)!=s_WorspaceMatrix[unique][N].end())
		return s_WorspaceMatrix[unique][N][M];
	
	if (s_WorspaceMatrix[unique].find(N)==s_WorspaceMatrix[unique].end())
		s_WorspaceMatrix[unique][N] = std::map<int, gsl_matrix*>();
	
	s_WorspaceMatrix[unique][N][M] = gsl_matrix_alloc(N, M);
	return s_WorspaceMatrix[unique][N][M];
}

// FFT stuff
// need it for FFT - need to modify
void LinAlgMatrixReal::WriteRow(int row, gsl_complex_packed_array data)
{
	for (int i=0; i<n; i++)
		gsl_matrix_set(m_pData, row, i, data[2*i + 0]);	
}

void LinAlgMatrixReal::WriteCol(int col, gsl_complex_packed_array data)
{
	for (int i=0; i<m; i++)
		gsl_matrix_set(m_pData, i, col, data[2*i + 0]);		
}

void LinAlgMatrixReal::ReadRow(int row, gsl_complex_packed_array data)
{
	double * rowData = gsl_matrix_ptr(m_pData, row, 0);
	gsl_fft_real_unpack (rowData, data, 1, n);
}

void LinAlgMatrixReal::ReadCol(int col, gsl_complex_packed_array data)
{
	for (int j=0; j<m; j++)
	{
		data[2*j+0] = (*this)(j, col);
		data[2*j+1] = 0;
	}		
}

#endif




