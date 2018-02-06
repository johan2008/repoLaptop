
#include "LinAlg.h"
class LinAlgMatrixComplex;
class LinAlgVectorReal;

#ifndef __LINALG_MATRIX_REAL_H
#define __LINALG_MATRIX_REAL_H

class LinAlgMatrixReal
{
	public:
		friend class LinAlgVectorReal;
		friend class LinAlgMatrixComplex;	
	
	// Initialization	
		LinAlgMatrixReal(int m, int n, double initValue = 0);
		LinAlgMatrixReal(LinAlgMatrixReal & srcMatrix);		
		virtual ~LinAlgMatrixReal();
	
	// Matrix operations, changing matrix
		LinAlgMatrixReal operator * (LinAlgMatrixReal & m2);	
		virtual void SetValues(double value);
		virtual void MultiplyBy(double value);
		virtual void Add(LinAlgMatrixReal & otherMatrix);	
		LinAlgMatrixReal & operator=(LinAlgMatrixReal & src);
		void AddVal(int i, int j, double val);
		double & operator()(int i, int j);	
	
	// Matrix operations, save result to a target
		virtual void Transpose(LinAlgMatrixReal & trg);	
		virtual void PointwiseMultiply(LinAlgMatrixReal & m2, LinAlgMatrixReal & trg);
		virtual void Multiply(LinAlgMatrixReal & m2, LinAlgMatrixReal & trg);
		virtual void Multiply(LinAlgVectorReal & v, LinAlgVectorReal & vTrg);
	
	// System of Linear Equations
		virtual void Inverse();
		virtual void Solve(LinAlgVectorReal & b, LinAlgVectorReal &x);
		virtual void SolveLUDestructive(LinAlgVectorReal & b, LinAlgVectorReal &x);
		virtual void LeastSquaresQRDestructive(LinAlgVectorReal & b, LinAlgVectorReal &x);
		
		virtual void DecomposeCholesky();
		virtual void SolveDecomposedCholesky(LinAlgVectorReal & b, LinAlgVectorReal &x);
	
	// eigendecomposition
		virtual void EigendecompositionSymmetric(LinAlgVectorReal & eigenValues, 
												 LinAlgMatrixReal & eigenVectors);
		virtual void RowNormalize();
	
	// svd decomposition
		virtual void SVDDecomposeDestructive(LinAlgMatrixReal & U, 
											 LinAlgVectorReal & S, LinAlgMatrixReal & Vt);	

	// Getters	
		virtual double GetMax(int * saveRow=NULL, int * saveCol=NULL);
		virtual double GetMin(int * saveRow=NULL, int * saveCol=NULL);
		virtual double Norm(double p);		//pth norm		
		virtual int Rows();
		virtual int Cols();
	
	// Visualizing / Type Conversion matrix
		virtual void PrintMatrix(double absMax=0);
		virtual void PrintMatlab();	
		virtual bool WriteMatrixBin(const char * outputFile);
		virtual bool ReadMatrixBin(const char * inputFile);
		virtual bool WriteMatrixASCII(const char * outputFile);
		virtual bool ReadMatrixASCII(const char * inputFile);
		
	protected:		
		int m, n;
		double * m_aRealData;
		bool m_bFreeData;
	
		static std::map<int, std::map<int, double *> > m_LapackWorkspaces;
		static std::map<int, std::map<int, int *> > m_LapackIntWorkspaces;
	
		static double * GetLapackWorkspace(int arraySize, int id);
		static int * GetLapackIntWorkspace(int arraySize, int id);
	
#ifdef WITH_GSL
		gsl_matrix * m_pGSLData;
		// Native operations 	
		virtual void WriteRow(int row, gsl_complex_packed_array data);
		virtual void WriteCol(int col, gsl_complex_packed_array data);
		virtual void ReadRow(int row, gsl_complex_packed_array data);
		virtual void ReadCol(int col, gsl_complex_packed_array data);
		
		static gsl_vector * GetWorkspaceVector(int N, int unique=0);
		static gsl_matrix * GetWorkspaceMatrix(int N, int M, int unique=0);
		
		static std::vector<std::map<int, gsl_vector*> > s_WorspaceVector;
		static std::vector<std::map<int, std::map<int, gsl_matrix*> > > s_WorspaceMatrix;	
		static std::map<int, gsl_eigen_symmv_workspace *> m_EigenNonSymmWorkspace;	
#endif
	
};

#endif
