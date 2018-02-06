#include "LinAlgMatrixSparseReal.h"
#include "LinAlgVectorReal.h"
#include "LinAlgMatrixReal.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include "Sortable.h"


LinAlgMatrixSparseReal::LinAlgMatrixSparseReal(int m, int n, int formats, int maxNonzero)
{
	m_pSVDLibData = NULL;
	m_Formats = formats;
	this->m = m;
	this->n = n;
	m_bFreeData = true;
	int values = 1;	//init values
	int triplet = 1; //triplet
	if (IsFormatSupported(SUPPORT_LIBSVD))
	{
		assert(maxNonzero>0);	// need to know exact number of non-zero entries
		m_pSVDLibData = svdNewSMat(m, n, maxNonzero);
		m_iCurrSVDLibDataCol = -1;
		m_iCurrSVDLibDataRow = -1;
		m_iCurrentSVDLibDataN = 0;
	}
	if (IsFormatSupported(SUPPORT_CX_SPARSE))
	{
		if (maxNonzero==-1)
			maxNonzero = m * 10;
		m_pData = cs_di_spalloc (m, n, maxNonzero, values, triplet);
	}
	if (IsFormatSupported(SUPPORT_ALL_DYNAMIC_CAST))
	{
		m_ArrayOfColsMapRows = new std::map<int, double>[n]; 
	}
}

LinAlgMatrixSparseReal::~LinAlgMatrixSparseReal()
{
	if (m_bFreeData && IsFormatSupported(SUPPORT_CX_SPARSE))
		cs_di_spfree(m_pData);
	if (m_bFreeData && IsFormatSupported(SUPPORT_LIBSVD))
		svdFreeSMat(m_pSVDLibData);
	if (m_bFreeData && IsFormatSupported(SUPPORT_ALL_DYNAMIC_CAST))
		delete [] m_ArrayOfColsMapRows;
}

bool LinAlgMatrixSparseReal::IsFormatSupported(SupportedMatrixFormats format)
{
	return (format & m_Formats)!=0;
}

// Matrix operations, changing matrix
void LinAlgMatrixSparseReal::AddVal(int i, int j, double val)
{
	assert(i<m && j<n && i>=0 && j>=0);
	if (IsFormatSupported(SUPPORT_CX_SPARSE))
		cs_di_entry (m_pData, i, j, val) ;
	if (IsFormatSupported(SUPPORT_ALL_DYNAMIC_CAST))
	{
		std::map<int, double>::iterator iter = m_ArrayOfColsMapRows[j].find(i);
		if (iter == m_ArrayOfColsMapRows[j].end())
			(m_ArrayOfColsMapRows[j])[i] = val;
		else
			iter->second += val;
	}
	if (IsFormatSupported(SUPPORT_LIBSVD))
	{
		assert(m_pSVDLibData!=NULL);
		if (j < m_iCurrSVDLibDataCol )
		{
			std::cout<<"[ERROR] Sparse Matrix: Columns can only increase."<<std::endl;
			assert(false);
		}
		else if (j == m_iCurrSVDLibDataCol && i <= m_iCurrSVDLibDataRow)
		{
			std::cout<<"[ERROR] Sparse Matrix: Rows can only increase for the same column.";
			assert(false);
		}
		
		assert(m_iCurrentSVDLibDataN < m_pSVDLibData->vals);	// can fit non-zero entries
		if (j > m_iCurrSVDLibDataCol)	// column increased
		{
			for (int c=m_iCurrSVDLibDataCol+1; c<=j; c++)
				m_pSVDLibData->pointr[c] = m_iCurrentSVDLibDataN;
			m_iCurrSVDLibDataCol = j;
			m_iCurrSVDLibDataRow = i;
		}

		// add entry
		m_pSVDLibData->rowind[m_iCurrentSVDLibDataN] = i;
		m_pSVDLibData->value[m_iCurrentSVDLibDataN] = val;
		m_iCurrSVDLibDataRow = i;
		m_iCurrentSVDLibDataN++;
	}
}

double LinAlgMatrixSparseReal::GetVal(int i, int j, double val)
{
	assert(IsFormatSupported(SUPPORT_ALL));
	std::map<int, double>::iterator iter = m_ArrayOfColsMapRows[j].find(i);
	if (iter==m_ArrayOfColsMapRows[j].end())
		return 0;
	else
		return iter->second;
}

void LinAlgMatrixSparseReal::RowNormalize()
{
	assert(IsFormatSupported(SUPPORT_ALL));
	double * rowNorms = new double[m];
	for (int i=0; i<m; i++)
		rowNorms[i] = 0;
	
	for (int i=0; i<n; i++)
	{
		for (std::map<int, double>::iterator iter = m_ArrayOfColsMapRows[i].begin(); 
			 iter != m_ArrayOfColsMapRows[i].end(); iter++)
		{
			assert(iter->first < m);
			rowNorms[iter->first] += iter->second;
		}
	}		
	for (int i=0; i<n; i++)
	{
		for (std::map<int, double>::iterator iter = m_ArrayOfColsMapRows[i].begin(); 
			 iter != m_ArrayOfColsMapRows[i].end(); iter++)
			iter->second /= rowNorms[iter->first];
	}		
	delete [] rowNorms;
}


// System of Linear Equations
void LinAlgMatrixSparseReal::Solve(const LinAlgVectorReal & b, LinAlgVectorReal &x)
{
	assert(IsFormatSupported(SUPPORT_CX_SPARSE));
	assert(b.Rows()==m && b.Rows()==n 
		   && b.Rows()==x.Rows());
	//double * deleteMe = b.ToArray();
	
	double * deleteMe = new double[b.Rows()];
	for (int i=0; i<b.Rows(); i++)
		deleteMe[i] = b.GetValue(i);
	
	int order = 1;	// order: 0, 1, 2, 3 (?)
	double tol = 0;
	int ok;
	cs_di * deleteMeData = cs_di_compress(m_pData);
	assert(deleteMeData!=NULL);
	assert(cs_dupl(deleteMeData));
	ok = cs_lusol (order, deleteMeData, deleteMe, tol);
	
	assert(ok);

	for (int i=0; i<m; i++)
		x(i) = deleteMe[i];

	//delete [] deleteMe;
	cs_di_spfree(deleteMeData);
}

void LinAlgMatrixSparseReal::FillSVDLibData()
{
	assert(IsFormatSupported(SUPPORT_ALL_DYNAMIC_CAST));
	assert(m_pSVDLibData==NULL);	// for now only cast once
	int nz=0;
	for (int i=0; i<n; i++)
		nz += (int)m_ArrayOfColsMapRows[i].size();
	m_pSVDLibData = svdNewSMat(m, n, nz);
	m_iCurrentSVDLibDataN = 0;
	
	for (m_iCurrSVDLibDataCol=0; m_iCurrSVDLibDataCol<n; m_iCurrSVDLibDataCol++)
	{
		m_pSVDLibData->pointr[m_iCurrSVDLibDataCol] = m_iCurrentSVDLibDataN;		
		std::map<int, double> & currCol = m_ArrayOfColsMapRows[m_iCurrSVDLibDataCol];
		
		int row = -1;
		for (std::map<int, double>::iterator iter = currCol.begin(); 
			 iter != currCol.end(); iter++)
		{
			m_pSVDLibData->rowind[m_iCurrentSVDLibDataN] = iter->first;
			m_pSVDLibData->value[m_iCurrentSVDLibDataN] = iter->second;
			assert(row < iter->first);
			row = iter->first;	// make sure rows increase;
			m_iCurrentSVDLibDataN++;
		}
	}
}

void LinAlgMatrixSparseReal::WriteMatrixASCII(const char * rows, 
											  const char * cols,
											  const char * vals)
{
	assert(IsFormatSupported(SUPPORT_ALL_DYNAMIC_CAST));
	std::ofstream textStreamRows(rows);
	std::ofstream textStreamCols(cols);
	std::ofstream textStreamVals(vals);	
	assert(textStreamRows.is_open());
	assert(textStreamCols.is_open());
	assert(textStreamVals.is_open());

	for (int i=0; i<n; i++)
	{
		std::map<int, double> & currCol = m_ArrayOfColsMapRows[i];
		
		for (std::map<int, double>::iterator iter = currCol.begin(); 
			 iter != currCol.end(); iter++)
		{
			textStreamRows<<iter->first<<" ";
			textStreamCols<<i<<" ";
			textStreamVals<<iter->second<<" ";
		}
	}
}

void LinAlgMatrixSparseReal::WriteMatrixBin(const char * filename)
{
	assert(IsFormatSupported(SUPPORT_ALL_DYNAMIC_CAST));
	int nz=0;
	for (int i=0; i<n; i++)
	{
		for (std::map<int, double>::iterator iter = m_ArrayOfColsMapRows[i].begin(); 
			 iter != m_ArrayOfColsMapRows[i].end(); iter++)
			nz++;
	}
	int * rows = new int[nz];
	int * cols = new int[nz];
	double * vals = new double[nz];	
	nz = 0;
	for (int i=0; i<n; i++)
	{
		for (std::map<int, double>::iterator iter = m_ArrayOfColsMapRows[i].begin(); 
			 iter != m_ArrayOfColsMapRows[i].end(); iter++)
		{
			rows[nz] = iter->first;
			cols[nz] = i;
			vals[nz] = iter->second;			
			nz++;
		}
	}
	FILE * file = fopen(filename, "wb");
	assert(file!=NULL);
	fwrite(&nz, 1, sizeof(int), file);
	fwrite(rows, nz, sizeof(int), file);
	fwrite(cols, nz, sizeof(int), file);
	fwrite(vals, nz, sizeof(double), file);	
	fclose(file);
	delete [] rows;
	delete [] cols;
	delete [] vals;
}

bool LinAlgMatrixSparseReal::ReadMatrixBin(const char * filename)
{
	assert(IsFormatSupported(SUPPORT_ALL_DYNAMIC_CAST));
	FILE * file = fopen(filename, "rb");
	if (file==NULL)
		return false;

	int nz;
	fread(&nz, 1, sizeof(int), file);
	int * rows = new int[nz];
	int * cols = new int[nz];
	double * vals = new double[nz];	
	fread(rows, nz, sizeof(int), file);
	fread(cols, nz, sizeof(int), file);
	fread(vals, nz, sizeof(double), file);	
	fclose(file);

	for (int i=0; i<nz; i++)
		AddVal(rows[i], cols[i], vals[i]);
	
	delete [] rows;
	delete [] cols;
	delete [] vals;
	
	return true;
}


	// SPARSE SVD
void LinAlgMatrixSparseReal::EigendecompositionSymmetric(LinAlgVectorReal & eigenValues, 
														 LinAlgMatrixReal & eigenVectors,
														 int dimensions)
{
	EigendecompositionSymmetricSVDLib(eigenValues, eigenVectors, dimensions);
}


void LinAlgMatrixSparseReal::EigendecompositionSymmetricSVDLib(LinAlgVectorReal & eigenValues, 
															   LinAlgMatrixReal & eigenVectors,
															   int dimensions)
{
	if (!IsFormatSupported(SUPPORT_LIBSVD) && IsFormatSupported(SUPPORT_ALL_DYNAMIC_CAST))
		FillSVDLibData();
	else
		assert(IsFormatSupported(SUPPORT_LIBSVD));
	assert(m_pSVDLibData!=NULL);
	assert(m_iCurrentSVDLibDataN==m_pSVDLibData->vals);	// found all
	
	SVDRec svdRec = svdLAS2A(m_pSVDLibData, dimensions);
//	double end[2] = {-1.0e-30, 1.0e-30};
//	double kappa = 1e-6;
//	SVDRec svdRec = svdLAS2(m_pSVDLibData, dimensions, dimensions * 1000, end, kappa);
	
	int d = svdRec->d;	// dimensionality (rank)
	DMat Ut = svdRec->Ut;    // Transpose of left singular vectors. (d by m)  The vectors are the rows of Ut.
	double *S = svdRec->S;  // Array of singular values. (length d) 
	DMat Vt = svdRec->Vt;    // Transpose of right singular vectors. (d by n) The vectors are the rows of Vt. 
	
	if (eigenValues.Rows()!=d)
	{
		std::cout<<"[ERROR] Cannot fit eigenvalues: "<<eigenValues.Rows()<<" vs "<<d<<std::endl;
		assert(false);
	}
	
	std::vector<Sortable> sortedEigs;
	for (int i=0; i<d; i++)
		sortedEigs.push_back(Sortable(-vkAbs(S[i]), NULL, i));
	std::sort(sortedEigs.begin(), sortedEigs.end());
	
	for (int i=0; i<d; i++)
	{
		eigenValues(i) = S[sortedEigs[i].id];
		//eigenValues(i) = S[i];
		std::cout<<"eval["<<i<<"] = "<<S[sortedEigs[i].id]<<std::endl;
	}
	
	int r = Vt->rows;
	int c = Vt->cols;	
	assert(r==d);
	assert(eigenVectors.Cols()==r && eigenVectors.Rows()==c);
	for (int i=0; i<d; i++)
	{
		//std::cout<<"evec["<<i<<"]:";
		for (int j=0; j<c; j++)
		{
			eigenVectors(j,i) = Vt->value[sortedEigs[i].id][j];
		//	std::cout<<eigenVectors(j,i)<<" ";
		}
		//std::cout<<std::endl;
	}
			//eigenVectors(j,i) = Vt->value[i][j];	
}

void LinAlgMatrixSparseReal::EigendecompositionSymmetricMatlab(LinAlgVectorReal & eigenValues, 
															   LinAlgMatrixReal & eigenVectors,
															   int dimensions, const char * uniquePrefixChar)
{
	std::cout<<"RUNNING MATLAB"<<std::endl;
	// MATLAB VERSION
	std::ostringstream uniquePrefixStream;
	uniquePrefixStream<<uniquePrefixChar<<rand()<<time(NULL)<<(long int)this;
	std::string uniquePrefix = uniquePrefixStream.str();
//	WriteMatrixASCII((uniquePrefix+".rows.txt").c_str(), 
//					 (uniquePrefix+".cols.txt").c_str(),
//					 (uniquePrefix+".vals.txt").c_str());
	WriteMatrixBin((uniquePrefix+".matrix.bin").c_str());
	
	
	std::ostringstream cmd;	
	cmd<<"matlab -nodisplay -nosplash -r \"findAndWriteEigenvectors('";
	cmd<<uniquePrefix<<"', "<<dimensions<<", "<<m<<", "<<n<<");\"";
	//std::cout<<"running "<<cmd.str().c_str()<<std::endl;
	system(cmd.str().c_str());

	bool readEigenval = eigenValues.ReadVectorBin((uniquePrefix+".eigenvalues.bin").c_str());
	bool readEigenvec = eigenVectors.ReadMatrixBin((uniquePrefix+".eigenvectors.bin").c_str());	

	if (!readEigenval || !readEigenvec)
	{
		std::cout<<"[WARNING] Could not find 'matlab' trying: /usr/local/matlab/bin/matlab";
		std::cout<<std::endl;
		const char * matlabExe = "/usr/local/matlab/bin/matlab";
		std::ostringstream cmd;	
		cmd<<matlabExe<<" -nodisplay -nosplash -r \"";
		cmd<<"findAndWriteEigenvectors('";
		cmd<<uniquePrefix<<"', "<<dimensions<<", "<<m<<", "<<n<<");\"";
		system(cmd.str().c_str());
		
		readEigenval = eigenValues.ReadVectorBin((uniquePrefix+".eigenvalues.bin").c_str());
		readEigenvec = eigenVectors.ReadMatrixBin((uniquePrefix+".eigenvectors.bin").c_str());	
	}
	
	assert(readEigenval && readEigenvec);
	std::ostringstream cmd2;	
	cmd2<<"rm "<<uniquePrefix<<"*";
	system(cmd2.str().c_str());	
}

