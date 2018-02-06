#ifdef WITH_GSL

#include "LinAlg.h"

#include "LinAlgMatrixComplex.h"
#include "LinAlgMatrixReal.h"

#ifndef __CALCULUS_H
#define __CALCULUS_H

class FFTManager
{
	public:	
		struct FFTWorkspace
		{			
			int n;
			gsl_fft_complex_wavetable * m_pTrigonometricTable;
			gsl_fft_complex_workspace * m_pWorkspace;			
			gsl_complex_packed_array m_RowArray;
		};
		void FFT2D(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg);
		void IFFT2D(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg);
		void FFT1DRow(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg);
		void IFFT1DRow(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg);
		void FFT1DCol(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg);
		void IFFT1DCol(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg);

	
		FFTWorkspace Get1DWorkspace(int n);
		void Clear1DWorkspace(int n);
	
		//testing
		static void TestFFT2D();
	
		std::map<int,  FFTWorkspace> m_WorkspaceTable;
};

class NumericalMethods
{
	public:
		static void IntegratedLagrangianWghts(std::vector<double> & smplLoc, std::vector<double> & w, int m);
		static void TestIntegration();
		static void TestIntegration2D();
	
	protected:
		static double IntegrateLagrangianInterpolant(int id, std::vector<double> & smplLoc, 
													 int m, int minID, int maxID, int interval);
};

#endif

#endif // if WITH_GSL
