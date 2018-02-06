class LinAlgMatrixComplex;
class LinAlgVectorReal;

#ifndef __LINALG_MATRIX_SPARSE_REAL_H
#define __LINALG_MATRIX_SPARSE_REAL_H

class LinAlgMatrixSparseReal
{
	public:
		enum SupportedMatrixFormats
		{
			SUPPORT_CX_SPARSE = 1,	// solve system of linear equations
			SUPPORT_LIBSVD = 2,		// eigendecomposition
			SUPPORT_ARPACK = 4,		// proper eigendecomposition
			SUPPORT_MATLAB = 8,		// very slow - but should be the most accurate
			SUPPORT_ALL_DYNAMIC_CAST = 16,
			SUPPORT_ALL = 0xFFFFFFFF
		};
		friend class LinAlgVectorReal;
		friend class LinAlgMatrixComplex;	
		friend class LinAlgMatrixReal;
	
	// Initialization	
		LinAlgMatrixSparseReal(int m, int n, int formats, int maxNonZero = -1);
		virtual ~LinAlgMatrixSparseReal();
	
	// Matrix operations, changing matrix
		LinAlgMatrixSparseReal operator * (LinAlgMatrixSparseReal & m2) const;	
		void AddVal(int i, int j, double val);
		double GetVal(int i, int j, double val);
		void RowNormalize();
	
	// System of Linear Equations
		virtual void Solve(const LinAlgVectorReal & b, LinAlgVectorReal &x);
		virtual void EigendecompositionSymmetric(LinAlgVectorReal & eigenValues, 
												 LinAlgMatrixReal & eigenVectors, 
												 int dimensions);
		virtual void EigendecompositionSymmetricSVDLib(LinAlgVectorReal & eigenValues, 
													   LinAlgMatrixReal & eigenVectors, 
													   int dimensions);
		virtual void EigendecompositionSymmetricMatlab(LinAlgVectorReal & eigenValues, 
													   LinAlgMatrixReal & eigenVectors, 
													   int dimensions, const char * uniquePrefix = "none");
	
		virtual void WriteMatrixASCII(const char * rows, 
									  const char * cols,
									  const char * vals);

		virtual void WriteMatrixBin(const char * filename);	
		virtual bool ReadMatrixBin(const char * filename);
	
	protected:	
		bool IsFormatSupported(SupportedMatrixFormats format);
		int m_Formats;
		int m, n;
		bool m_bFreeData;
		
		cs_di * m_pData;
	
		SMat m_pSVDLibData;
		int m_iCurrSVDLibDataCol;
		int m_iCurrSVDLibDataRow;
		int m_iCurrentSVDLibDataN;
	
	// my format for sparse
		std::map<int, double> * m_ArrayOfColsMapRows;
		void FillSVDLibData();
		void FillCXSparseFormat();
};

#endif
