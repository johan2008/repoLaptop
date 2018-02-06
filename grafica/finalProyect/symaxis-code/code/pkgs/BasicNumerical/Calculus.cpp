#ifdef WITH_GSL
#include "BasicNumerical/BasicNumerical.h"

//TODO: assert that matrices are 2D (if necessary)
void FFTManager::FFT2D(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg)
{
	//MATLAB: fft2D(f) = fft(fft(f).').' where .' is transpose	
	//F1 = fft(f') - do fft on each row of f
	//trg = fft(F1') - do fft on each row of F1

	FFTWorkspace workspace = Get1DWorkspace(src.Cols());
	
	// find F1=fft(f')
	for (int i=0; i<src.Rows(); i++)//for each row of src
	{
		//do forward fourier transform
		src.ReadRow(i, workspace.m_RowArray);
		gsl_fft_complex_forward(workspace.m_RowArray, 1, src.Cols(), 
								workspace.m_pTrigonometricTable, 
								workspace.m_pWorkspace);
		trg.WriteRow(i, workspace.m_RowArray);
	}
	
	workspace = Get1DWorkspace(src.Rows());
	// find fft(F1')
	for (int i=0; i<src.Cols(); i++) //for each col of F1
	{
		//do forward fourier transform
		trg.ReadCol(i, workspace.m_RowArray);
		gsl_fft_complex_forward(workspace.m_RowArray, 1, src.Rows(), 
								workspace.m_pTrigonometricTable, 
								workspace.m_pWorkspace);
		trg.WriteCol(i, workspace.m_RowArray);
	}
}

void FFTManager::IFFT2D(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg)
{
	//MATLAB: ifft2D(F) = ifft(ifft(F).').'
	FFTWorkspace workspace = Get1DWorkspace(src.Rows());
	
	for (int i=0; i<src.Cols(); i++) //for each col
	{
		src.ReadCol(i, workspace.m_RowArray);
		gsl_fft_complex_inverse(workspace.m_RowArray, 1, src.Rows(), 
								workspace.m_pTrigonometricTable, 
								workspace.m_pWorkspace);			
		trg.WriteCol(i, workspace.m_RowArray);
	}
	
//	std::cout<<"IFFT: for each col: "<<std::endl;
//	trg.PrintMatrix();
	
	workspace = Get1DWorkspace(src.Cols());
	for (int i=0; i<src.Rows(); i++) //for each row
	{
		trg.ReadRow(i, workspace.m_RowArray);
		gsl_fft_complex_inverse(workspace.m_RowArray, 1, src.Cols(), 
								workspace.m_pTrigonometricTable, 
								workspace.m_pWorkspace);			
		trg.WriteRow(i, workspace.m_RowArray);
	}	
//	std::cout<<"IFFT: for each row: "<<std::endl;
//	trg.PrintMatrix();
	
}

void FFTManager::FFT1DRow(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg)
{
	FFTWorkspace workspace = Get1DWorkspace(src.Cols());
	for (int i=0; i<src.Rows(); i++)//for each row of src
	{
		//do forward fourier transform
		src.ReadRow(i, workspace.m_RowArray);
		gsl_fft_complex_forward(workspace.m_RowArray, 1, src.Cols(), 
								workspace.m_pTrigonometricTable, 
								workspace.m_pWorkspace);
		trg.WriteRow(i, workspace.m_RowArray);
	}
}

void FFTManager::IFFT1DRow(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg)
{
	FFTWorkspace workspace = Get1DWorkspace(src.Cols());
	for (int i=0; i<src.Rows(); i++) //for each row
	{
		src.ReadRow(i, workspace.m_RowArray);
		gsl_fft_complex_inverse(workspace.m_RowArray, 1, src.Cols(), 
								workspace.m_pTrigonometricTable, 
								workspace.m_pWorkspace);			
		trg.WriteRow(i, workspace.m_RowArray);
	}	
}

void FFTManager::FFT1DCol(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg)
{
	FFTWorkspace workspace = Get1DWorkspace(src.Rows());
	// find fft(F1')
	for (int i=0; i<src.Cols(); i++) //for each col of F1
	{
		//do forward fourier transform
		src.ReadCol(i, workspace.m_RowArray);
		gsl_fft_complex_forward(workspace.m_RowArray, 1, src.Rows(), 
								workspace.m_pTrigonometricTable, 
								workspace.m_pWorkspace);
		trg.WriteCol(i, workspace.m_RowArray);
	}	
}

void FFTManager::IFFT1DCol(LinAlgMatrixComplex & src, LinAlgMatrixComplex & trg)
{
	FFTWorkspace workspace = Get1DWorkspace(src.Rows());
	for (int i=0; i<src.Cols(); i++) //for each col
	{
		src.ReadCol(i, workspace.m_RowArray);
		gsl_fft_complex_inverse(workspace.m_RowArray, 1, src.Rows(), 
								workspace.m_pTrigonometricTable, 
								workspace.m_pWorkspace);			
		trg.WriteCol(i, workspace.m_RowArray);
	}	
}



FFTManager::FFTWorkspace FFTManager::Get1DWorkspace(int n)
{
	if (m_WorkspaceTable.find(n)!=m_WorkspaceTable.end())
	{
		return m_WorkspaceTable[n];
	}
	else
	{
		FFTWorkspace workspace;
		workspace.n = n;
		workspace.m_RowArray = new double[n*2];
		workspace.m_pTrigonometricTable = gsl_fft_complex_wavetable_alloc(n);
		workspace.m_pWorkspace = gsl_fft_complex_workspace_alloc(n);
		m_WorkspaceTable[n] = workspace;
		return workspace;
	}
}

void FFTManager::Clear1DWorkspace(int n)
{
	FFTWorkspace & removeMe = m_WorkspaceTable[n];
	gsl_fft_complex_wavetable_free(removeMe.m_pTrigonometricTable);
	gsl_fft_complex_workspace_free(removeMe.m_pWorkspace);
	m_WorkspaceTable.erase(m_WorkspaceTable.find(n));
}

void FFTManager::TestFFT2D()
{
	FFTManager mng;
	LinAlgMatrixComplex f(4, 4);
	int baseInd = 0;
   	f(baseInd+0, baseInd+0).dat[0]= 4;	f(baseInd+0, baseInd+1).dat[0]= 0;	
		f(baseInd+0, baseInd+2).dat[0]= 0;	f(baseInd+0, baseInd+3).dat[0]= 0;			
	f(baseInd+1, baseInd+0).dat[0]= 0; f(baseInd+1, baseInd+1).dat[0]= 4;	
		f(baseInd+1, baseInd+2).dat[0]= 0;	f(baseInd+1, baseInd+3).dat[0]= 0;	
	
	f(baseInd+2, baseInd+0).dat[0]= 0;	f(baseInd+2, baseInd+1).dat[0]= 0;	
		f(baseInd+2, baseInd+2).dat[0]= 4;	f(baseInd+2, baseInd+3).dat[0]= 0;			
	f(baseInd+3, baseInd+0).dat[0]= 0;	f(baseInd+3, baseInd+1).dat[0]= 0;	
		f(baseInd+3, baseInd+2).dat[0]= 0;	f(baseInd+3, baseInd+3).dat[0]= 4;
	
	LinAlgMatrixComplex fftF(f.Rows(), f.Cols());
	LinAlgMatrixComplex ifftF(f.Rows(), f.Cols());
	mng.FFT2D(f, fftF);
//	std::cout<<"f="<<std::endl;
//	f.PrintMatrix();
//	std::cout<<"fft(f)="<<std::endl;
//	fftF.PrintMatrix();
	mng.IFFT2D(fftF, ifftF);
	std::cout<<"ifft(fft(f))="<<std::endl;
	ifftF.PrintMatrix();
}



void NumericalMethods::IntegratedLagrangianWghts(std::vector<double> & smplLoc, 
												 std::vector<double> & w, int m)
{
	assert(m*2+2 <= (int)smplLoc.size());
	// initialize weights to 0
	for (int i=0; i<(int)smplLoc.size(); i++)
		w.push_back(0.);
	
	int N = (int)smplLoc.size();
	// for each interval
	for (int interval = 0; interval < N-1; interval++)
	{
		// find neighborhood
		int minID = interval-m;
		int maxID = interval+1+m;
		if (interval-m<0)
		{
			minID = 0;
			maxID = 2*m+1;
		}
		else if (interval+m+1>=N)
		{
			minID = N - 2*m-2;
			maxID = N-1;
		}
		
//		std::cout<<"Considering interval "<<interval<<" : ["<<smplLoc[interval]<<", "<<smplLoc[interval+1]<<"]";
//		std::cout<<" l["<<minID<<", "<<maxID<<"]"<<std::endl;		
			
		// sum weights for current interval
		for (int j=minID; j<=maxID; j++)
		{
//			std::cout<<"\tCalculating Weight "<<j<<" Interval=["<<minID<<", "<<maxID<<"]"<<std::endl;
			double addedVal=IntegrateLagrangianInterpolant(j, smplLoc, m, minID, maxID, interval);
			w[j] += addedVal;
//			std::cout<<"\t===== Weight["<<j<<"] = "<<w[j]<<" Increased by "<<addedVal<<std::endl;
		}
	}
	
	return;
}

double NumericalMethods::IntegrateLagrangianInterpolant(int id, std::vector<double> & smplLoc, 
														int m, int minID, int maxID, int interval)
{
	assert(2*m+2 == maxID-minID+1);
	// ∫ l_id(x) dx
	
	double denominator = 1;
	double xi = smplLoc[id];
		// find (constant) denominator of Lagrangian polynomial: ∏ (xi-xj)
	for (int j=minID; j<=maxID; j++)
	{
		if (xi!=smplLoc[j])
			denominator *= (xi-smplLoc[j]);
	}
	
	//go over all pairs in nominator of Lagrangian polynomial: 
	//		π (x-xj) = ∑(x or x_k) * (x or x_k+1) ... ~ (0 or 1) * (0 or 1)...
	int power = 2*m+1;	//how many possible powers of x are there?
	double * polynomial = new double[power+1];	//for each power store the weight
	
//	std::cout<<"\t\tInitializing polynomial to 0: power="<<power<<std::endl;
	for (int i=0; i<power+1; i++)		// initialize to 0 (sum)
		polynomial[i] = 0;

	for (int binaryLoop=0; binaryLoop<pow(2, power); binaryLoop++)
	{
		// calculate number of 0 bits
		int num0Bits = 0;
		double monomialWeight = 1;
		for (int p=0; p<power; p++)
		{
			if (((binaryLoop>>p) & 0x1) == 0)	// selected x
				num0Bits++;
			else			// selected x_k+n
			{
				int trueP = minID+p;
				if (trueP >= id)	//skip x_i
					trueP++;
				monomialWeight *= -smplLoc[trueP];	// update current weight
//				std::cout<<"\t\t\tWeight Updated By "<<-smplLoc[trueP]<<" At "<<trueP<<std::endl;
			}
		}
		polynomial[num0Bits] += monomialWeight/denominator;
//		std::cout<<"\t\tMonomial["<<num0Bits<<"] = "<<polynomial[num0Bits]<<" Increased by "<<monomialWeight<<"/"<<denominator<<"="<<monomialWeight/denominator<<std::endl;
	}
	
	// integrate polynomial
	double upperPoint = smplLoc[interval+1];
	double lowerPoint = smplLoc[interval+0];	
	double upperValue = 0;
	double lowerValue = 0;
//	std::cout<<"\t\tIntegrating over ["<<lowerPoint<<", "<<upperPoint<<"]"<<std::endl;
	for (int i=0; i<power+1; i++)
	{
		upperValue += polynomial[i] / (i+1) * pow(upperPoint, i+1);		// W_i / (i+1) * x_{id}^(i+1)
		lowerValue += polynomial[i] / (i+1) * pow(lowerPoint, i+1);		// W_i / (i+1) * x_{id+1}^(i+1)
	}
//	std::cout<<"\t\tIntegral="<<upperValue<<"-"<<lowerValue<<"="<<(upperValue - lowerValue)<<std::endl;
	
	delete [] polynomial;
	
	return (upperValue - lowerValue);	
}

double testIntegrationF1(double x)
{
	return x*x + 1;
}

double trueValueIntegrationF1(double xmin, double xmax)
{
	return (pow(xmax, 3.)/3.+xmax)-(pow(xmin, 3.)/3.+xmin);
}

void NumericalMethods::TestIntegration()
{
	std::vector<double> samples;
	samples.push_back(1);
	samples.push_back(1.25);	
	samples.push_back(1.33);	
	samples.push_back(1.5);		
	samples.push_back(1.66);	
	samples.push_back(1.75);			
	samples.push_back(2);	
	std::vector<double> w;
	IntegratedLagrangianWghts(samples, w, 2);
	double integral = 0;
	std::cout<<"w = [ "<<std::flush;
	for (int i=0; i<(int)w.size(); i++)
		std::cout<<w[i]<<" "<<std::flush;
	std::cout<<"] "<<std::endl;

	for (int i=0; i<(int)samples.size(); i++)
		integral+=testIntegrationF1(samples[i]) * w[i];

	std::cout<<"Approximation="<<integral<<" True="<<trueValueIntegrationF1(1, 2)<<std::endl;
}

#endif

