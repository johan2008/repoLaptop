#include "LinAlg.h"
#include "LinAlgComplex.h"

#ifndef __LINALG_VECTOR_COMPLEX_H
#define __LINALG_VECTOR_COMPLEX_H

class LinAlgVectorComplex
{
	friend class FFTManager;
	friend class LinAlgMatrixReal;
	friend class LinAlgMatrixComplex;	
	friend class LinAlgVectorReal;
	
	public:
	// Initialization
		LinAlgVectorComplex(int m);
		virtual ~LinAlgVectorComplex();
	
	// Output
		void PrintMatlab();
	
	// Manipulation	
		LinAlgComplex & operator()(int i);
		LinAlgComplex GetValue(int i) const;
		void SetValue(int i, const LinAlgComplex & value);
	
	// Getters
		int Rows() const;
		
	protected:
		int m;
		LinAlgComplex * m_aComplexData;
	
#ifdef WITH_GSL
		//gsl_complex & operator()(int i);
		gsl_vector_complex * m_pGSLData;
#endif
};


#endif
