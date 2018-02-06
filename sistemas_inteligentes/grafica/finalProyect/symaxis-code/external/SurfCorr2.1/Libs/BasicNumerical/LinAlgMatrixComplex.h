
#include "LinAlg.h"
#include "LinAlgComplex.h"

class LinAlgMatrixReal;
class LinAlgVectorReal;
class LinAlgVectorComplex;

#ifndef __LINALG_MATRIX_COMPLEX_H
#define __LINALG_MATRIX_COMPLEX_H

class LinAlgMatrixComplex
{
	friend class FFTManager;
	friend class LinAlgMatrixReal;
	friend class LinAlgVectorReal;
	friend class LinAlgVectorComplex;
	
	public:
	// Initialization
		LinAlgMatrixComplex(int m, int n);
		LinAlgMatrixComplex(const LinAlgMatrixComplex & otherMatrix);	
		virtual ~LinAlgMatrixComplex();

	// Visualizing Matrix
		virtual void PrintMatrix(double absMax = 0 );
		virtual void PrintMatlab();
	
	// Getters/setters
		LinAlgComplex & operator()(int i, int j);		
		LinAlgComplex GetValue(int i, int j) const;
		void SetValue(int i, int j, const LinAlgComplex & value); 
		void SetValue(int i, int j, double realValue, double imgValue); 
	
	// Matrix operations, changing matrix
		virtual void Transpose();
		virtual void MultiplyByReal(double val);
		virtual void AddMatrix(const LinAlgMatrixComplex & m2);	
		LinAlgMatrixComplex & operator=(LinAlgMatrixComplex & src);	
		virtual double Norm();
	
	// Matrix operations, save result to a target
		virtual void PointwiseMultiply(LinAlgMatrixComplex & m2, LinAlgMatrixComplex & trg);
		virtual void Multiply(LinAlgMatrixComplex & m2, LinAlgMatrixComplex & trg);
		virtual void Multiply(LinAlgVectorComplex & v, LinAlgVectorComplex & trg);
		virtual void Conj(LinAlgMatrixComplex & trg);		
		virtual void Transpose(LinAlgMatrixComplex & trg);
		virtual void Log(LinAlgMatrixComplex & trg);
		virtual void Sqrt(LinAlgMatrixComplex & trg);	
		virtual void Exp(LinAlgMatrixComplex & trg);
	
	// System of Linear Equations
		virtual void Inverse();
		virtual void Solve(LinAlgVectorComplex & b, LinAlgVectorComplex &x);
		virtual void SolveLUDestructive(LinAlgVectorComplex & b, LinAlgVectorComplex &x);
		virtual void LeastSquaresQRDestructive(LinAlgVectorComplex & b, LinAlgVectorComplex &x);

	// singular value decomposition
		virtual void SVDDecomposeDestructive(LinAlgMatrixComplex & U, 
											 LinAlgVectorComplex & S, LinAlgMatrixComplex & Vt);

	// Getters
		virtual double GetMax() const;
		virtual double GetMin() const;	
		virtual int Rows() const;
		virtual int Cols() const;
			
	protected:
		LinAlgComplex * m_aComplexData;
		int m, n;
		static void SortEigenvalues(LinAlgVectorComplex & eigenValues, 
									LinAlgMatrixComplex & eigenVectors);
		
		static std::map<int, std::map<int, LinAlgComplex *> > m_LapackComplexWorkspaces;
		static std::map<int, std::map<int, double *> > m_LapackWorkspaces;
		static std::map<int, std::map<int, int *> > m_LapackIntWorkspaces;
		
		static LinAlgComplex * GetLapackComplexWorkspace(int arraySize, int id);
		static double * GetLapackWorkspace(int arraySize, int id);
		static int * GetLapackIntWorkspace(int arraySize, int id);
	
	private:
		//static void PrintMatlab(gsl_matrix_complex * data);
		
		static double s_fEpsilonForLogExp;
		static double s_fQForLogExp;
	
#ifdef WITH_GSL
		// used internally
		gsl_complex & operator()(int i, int j);		
	
		// Used by FFTManager
		virtual void WriteRow(int row, gsl_complex_packed_array data);
		virtual void WriteCol(int col, gsl_complex_packed_array data);
		virtual void ReadRow(int row, gsl_complex_packed_array data);
		virtual void ReadCol(int col, gsl_complex_packed_array data);
	
		// Used by LinAlgMatrixComplex
		gsl_matrix_complex * m_pGSLData;
	
		gsl_matrix_complex * m_pWorkMatrixMxN;
		gsl_matrix_complex * m_pWorkMatrixNxN;
		gsl_matrix_complex * m_pWorkMatrixNxN_additional;	
		gsl_matrix_complex * m_pWorkMatrixNxM;	
#endif
	
};


#endif

