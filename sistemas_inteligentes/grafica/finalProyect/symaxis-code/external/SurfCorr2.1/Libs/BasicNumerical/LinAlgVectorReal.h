#include "LinAlg.h"


#ifndef __LINALG_VECTOR_REAL_H
#define __LINALG_VECTOR_REAL_H

class LinAlgVectorReal
{
	friend class FFTManager;
	friend class LinAlgMatrixReal;
	
	public:
	// Initialization
		LinAlgVectorReal(int m, double initVal = 0);
		LinAlgVectorReal(std::vector<double> & inVector);
		virtual ~LinAlgVectorReal();
	
	// Type Conversion
		std::vector<double> ToVector();

	// Manipulation
		double & operator()(int i);
		double GetValue(int i) const;
		void SetValue(int i, double value);
		LinAlgVectorReal & operator = (LinAlgVectorReal & vec);
		void Multiply(double val);
		void Add(LinAlgVectorReal & v2);
		void SetAll(double val);	
	
	// Visualization
		void PrintMatlab(double absMin = 0) const;	
		bool WriteVectorBin(const char * outputFile);
		bool ReadVectorBin(const char * inputFile);
		bool WriteVectorASCII(const char * outputFile);
		bool ReadVectorASCII(const char * inputFile);
	
	// Getters
		double Norm(double p) const;
		int Rows() const;

	protected:
		bool m_bFreeData;
		int m;

		double * m_aRealData;
		
#ifdef WITH_GSL
		gsl_vector * m_pGSLData;
#endif
};

#endif
