#include "SurfacePerVertexValue.h"
#include "SampledSurface.h"
#include "VkFunctions.h"
#include "Sortable.h"
#include "SurfaceDistance.h"
#include "DistanceOnTheFly.h"
#include "AnalysisStats.h"
#include <algorithm>

std::map<int, double *> SurfacePerVertexValue::s_KnownValuesArrays;
std::map<int, double *> SurfacePerVertexValue::s_NormalizationArrays;
std::map<SampledSurface *, std::vector<SurfacePerVertexValue::PrecomputedWeights *> > 
								SurfacePerVertexValue::s_InterpolationWeights;

SurfacePerVertexValue::SurfacePerVertexValue(SampledSurface * surface, bool allocateAsNeeded)
{
	m_pSurface = surface;
	
	m_bAllVertexValuesAllocated = !allocateAsNeeded;
	m_PerVertexValues = NULL;
	m_PerVertexValuesExist = NULL;
	if (m_bAllVertexValuesAllocated)
		AllocatePerVertexStorage();
	
	m_iNumValuesNotSet = m_pSurface->GetMesh()->NVertices();

	m_fMinValue = FLT_MAX;
	m_fMaxValue = -FLT_MAX;	
}

SurfacePerVertexValue::~SurfacePerVertexValue()
{
	if (m_PerVertexValues!=NULL)
		delete [] m_PerVertexValues;
	
	if (m_PerVertexValuesExist!=NULL)
		delete [] m_PerVertexValuesExist;
}

void SurfacePerVertexValue::ClearValues()
{
	m_iNumValuesNotSet =m_pSurface->GetMesh()->NVertices(); 
	if (m_bAllVertexValuesAllocated)
		for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
			m_PerVertexValuesExist[i] = false;
	else
		m_VertexToValueMap.clear();
}

void SurfacePerVertexValue::AllocatePerVertexStorage()
{
//	std::cout<<"[WARNING] Allocating per vertex storage"<<std::endl;
//	std::cout<<"*"<<std::flush;
	m_bAllVertexValuesAllocated = true;
	int nvert = m_pSurface->GetMesh()->NVertices();	

	if (m_PerVertexValues==NULL)
		m_PerVertexValues = new double[nvert];
	
	if (m_PerVertexValuesExist==NULL)
		m_PerVertexValuesExist = new bool[nvert];	
	
	for (int i=0; i<nvert; i++)
	{
		m_PerVertexValues[i] = 0;
		m_PerVertexValuesExist[i] = false;
	}	
}

void SurfacePerVertexValue::AllocatePerVertexStorageFillFromKnown()
{
	if (m_bAllVertexValuesAllocated)
		return;
	
	m_bAllVertexValuesAllocated = true;
	AllocatePerVertexStorage();
	m_iNumValuesNotSet = m_pSurface->GetMesh()->NVertices();
	m_fMaxValue = -FLT_MAX;
	m_fMinValue = FLT_MAX;	
	for (std::map<int, double>::iterator iter = m_VertexToValueMap.begin(); 
		 iter != m_VertexToValueMap.end(); iter++)
		SetVertexValue(iter->first, iter->second);
}

void SurfacePerVertexValue::SetVertexValue(int vertexID, double value)
{
	m_fMinValue = vkMin(value, m_fMinValue);
	m_fMaxValue = vkMax(value, m_fMaxValue);	

	if (m_bAllVertexValuesAllocated)
	{
		m_PerVertexValues[vertexID] = value;
		if (!m_PerVertexValuesExist[vertexID])
			m_iNumValuesNotSet--;
		m_PerVertexValuesExist[vertexID] = true;
	}
	else
	{
		std::map<int, double>::iterator iter = m_VertexToValueMap.find(vertexID);
		
		if (iter == m_VertexToValueMap.end())
		{
			m_iNumValuesNotSet--;
			m_VertexToValueMap[vertexID] = value;
		}
		else
			iter->second = value;
	}
}

void SurfacePerVertexValue::SetVertexValue(const SurfaceSample & atVertexSample, double value)
{
	bool atVertex;
	int vertexID = atVertexSample.NearestVertex(&atVertex);
	assert(atVertex);
	SetVertexValue(vertexID, value);
}

double SurfacePerVertexValue::GetInterpolatedValue(const SurfaceSample & samp)
{
	bool atVertex;
	int vertexID = samp.NearestVertex(&atVertex);
	if (atVertex)
	{
		if (m_bAllVertexValuesAllocated)
		{
			assert(m_PerVertexValuesExist[vertexID]);
			return m_PerVertexValues[vertexID];
		}
		else
		{
			std::map<int, double>::iterator iter = m_VertexToValueMap.find(vertexID);
			assert(iter!=m_VertexToValueMap.end());
			return iter->second;
		}
	}
	else
	{
		if (m_bAllVertexValuesAllocated)
		{
			assert(m_PerVertexValuesExist[samp.VertID(0)]);
			assert(m_PerVertexValuesExist[samp.VertID(1)]);
			assert(m_PerVertexValuesExist[samp.VertID(2)]);
			return samp.Interpolate(m_PerVertexValues[samp.VertID(0)], 
									m_PerVertexValues[samp.VertID(1)], 
									m_PerVertexValues[samp.VertID(2)]);
		}
		else
		{
			std::map<int, double>::iterator iter1 = m_VertexToValueMap.find(vertexID);
			assert(iter1 != m_VertexToValueMap.end());
			std::map<int, double>::iterator iter2 = m_VertexToValueMap.find(vertexID);
			assert(iter2 != m_VertexToValueMap.end());			
			std::map<int, double>::iterator iter3 = m_VertexToValueMap.find(vertexID);
			assert(iter3 != m_VertexToValueMap.end());			
			
			return samp.Interpolate(iter1->second, iter2->second, iter3->second);
		}
		
	}
}

SurfaceSample SurfacePerVertexValue::NextAtVertexSampleRequiredForInterpolation(const SurfaceSample & samp)
{
	int vertexID = NextVertexRequiredForInterpolation(samp);
	if (vertexID == -1)
		return SurfaceSample();
	else 
		return SurfaceSample(vertexID, m_pSurface->GetMesh());
}

int SurfacePerVertexValue::NextVertexRequiredForInterpolation(const SurfaceSample & samp)
{
	bool atVertex;
	int vertexID = samp.NearestVertex(&atVertex);
	if (atVertex)
	{
		if (m_bAllVertexValuesAllocated)
			return m_PerVertexValuesExist[vertexID] ? -1 : vertexID;
		else
			return (m_VertexToValueMap.find(vertexID)!=m_VertexToValueMap.end()) ? -1 : vertexID;
	}
	else
		for (int i=0; i<3; i++)
			if ((m_bAllVertexValuesAllocated && !m_PerVertexValuesExist[samp.VertID(i)])
				|| (!m_bAllVertexValuesAllocated 
					&& m_VertexToValueMap.find(samp.VertID(i))==m_VertexToValueMap.end()))
				return samp.VertID(i);
	return -1;
}

void SurfacePerVertexValue::SmoothSimpleLaplacian(int numIterations)
{
	AllocatePerVertexStorageFillFromKnown();
	assert(AllValuesSet());
	R3Mesh * mesh = m_pSurface->GetMesh();
	double * vals = GetKnownValueArray(mesh->NVertices());
	for (int i=0; i<numIterations; i++)
	{
		for (int v=0; v<mesh->NVertices(); v++)	// mean shift
		{
			double norm=0;	
			vals[v] = 0;
			
			R3MeshVertex * vp = mesh->Vertex(v);
			for (int nb=0; nb<mesh->VertexValence(vp); nb++)
			{
				int otherV = mesh->VertexID(mesh->VertexAcrossEdge(mesh->EdgeOnVertex(vp, nb), vp));
				double weight = 1. / mesh->EdgeLength(mesh->EdgeOnVertex(vp, nb));
				vals[v] += m_PerVertexValues[otherV] * weight;
				norm+=weight;
			}
			vals[v] = 0.5 * vals[v] / norm + 0.5 * m_PerVertexValues[v];
		}
		
		for (int v=0; v<mesh->NVertices(); v++)	// re-assign values
			m_PerVertexValues[v] = vals[v];
	}	
}

void SurfacePerVertexValue::SmoothGaussianFromKnown(double sigma)
{
	AllocatePerVertexStorageFillFromKnown();	
	SurfaceSampleSet knownSet;
	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
		if (m_PerVertexValuesExist[i])
			knownSet.AddSample(SurfaceSample(i, m_pSurface->GetMesh()));
	SmoothGaussianFromKnown(&knownSet, sigma);
}

void SurfacePerVertexValue::SmoothGaussianFromKnown(const VKString & knownSetName, double sigma)
{
	std::cout<<"[WARNING] This version of the function SmoothGaussianFromKnown(VKString, double) ";
	std::cout<<"does not support splatting - might be slow. "<<std::endl;
	AllocatePerVertexStorageFillFromKnown();
	SurfaceSampleSet * knownSet = m_pSurface->GetSampleSet(knownSetName);
	assert(knownSet!=NULL);
	PrecomputedWeights * weights = GetInterpolationWeighting(knownSetName, m_pSurface, 
															 ReestimateSigma(sigma, knownSet->NumSamples()));

	double * knownValues = GetKnownValueArray(knownSet->NumSamples());
	for (int i=0; i<knownSet->NumSamples(); i++)
		knownValues[i] = GetInterpolatedValue(knownSet->GetSample(i));

	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
	{
		double * perSampleWeights=NULL;
		int * sampleIDs=NULL;
		int numSamples=0;
		weights->GetSampleWeighting(i, &perSampleWeights, &sampleIDs, &numSamples);
		double value = 0;
		for (int j=0; j<numSamples; j++)
			value += knownValues[sampleIDs[j]] * perSampleWeights[j];
		SetVertexValue(i, value);
	}
	
	// SmoothGaussianFromKnown(m_pSurface->GetSampleSet(knownSet), sigma);
}

void SurfacePerVertexValue::SmoothGaussianFromKnown(SurfaceSampleSet * knownSet, double sigma)
{
	AllocatePerVertexStorageFillFromKnown();	
	assert(knownSet!=NULL);
	
	// should do splatting if too many elements in a known set
	int numOperations = knownSet->NumSamples() * m_pSurface->GetMesh()->NVertices();
	if (numOperations > 62500)
	{
		//std::cout<<"[WARNING] Smoothing is too expensive ("<<knownSet->NumSamples()<<" x ";
		//std::cout<<m_pSurface->GetMesh()->NVertices()<<"). Switching to splatting."<<std::endl;
		double splattingRadius = EstimateRadiusForSplatting(knownSet->NumSamples());
		bool covered = SmoothGaussianFromKnownBySplatting(knownSet, sigma, splattingRadius);
		if (!covered)
		{
			//std::cout<<"[WARNING] Could not cover with first splatting - increasing the radius"<<std::endl;
			for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
				m_PerVertexValuesExist[i] = false;
			for (int i=0; i<knownSet->NumSamples(); i++)
				m_PerVertexValuesExist[knownSet->GetSample(i).NearestVertex()] = true;
			m_iNumValuesNotSet = m_pSurface->GetMesh()->NVertices() - knownSet->NumSamples();
			covered = SmoothGaussianFromKnownBySplatting(knownSet, sigma, splattingRadius*10);
			if (!covered)
			{
				std::cout<<"[WARNING] Could not cover using splatting - smoothing will be slow"<<std::endl;
				for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
					m_PerVertexValuesExist[i] = false;
				for (int i=0; i<knownSet->NumSamples(); i++)
					m_PerVertexValuesExist[knownSet->GetSample(i).NearestVertex()] = true;
				m_iNumValuesNotSet = m_pSurface->GetMesh()->NVertices() - knownSet->NumSamples();
				SmoothGaussianFromKnownExact(knownSet, sigma);
			}
		}
	}
	else
	{
		SmoothGaussianFromKnownExact(knownSet, sigma);
	}
}

double SurfacePerVertexValue::EstimateRadiusForSplatting(int numSamples)
{
	double areaPerSample = m_pSurface->Area() / numSamples;
	double distance = sqrt(areaPerSample/3.14159265);
	return distance * 3.;
}

bool SurfacePerVertexValue::SmoothGaussianFromKnownBySplatting(SurfaceSampleSet * knownSet, 
															   double sigma, double splatRadius)
{
	AllocatePerVertexStorageFillFromKnown();		
	double * knownValues = GetKnownValueArray(knownSet->NumSamples());
	double * normValues = GetNormalizationArray(m_pSurface->GetMesh()->NVertices());
	DistanceOnTheFly * tempDist = GetSurfaceDistance(splatRadius);

	double geodesicNorm = sqrt(m_pSurface->Area());	
	sigma = ReestimateSigma(sigma, knownSet->NumSamples());
	double sigma2 = sigma*sigma;	
		
	for (int i=0; i<knownSet->NumSamples(); i++)
		knownValues[i] = GetInterpolatedValue(knownSet->GetSample(i));
	
	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
	{
		m_PerVertexValues[i] = 0;
		m_PerVertexValuesExist[i] = false;
		normValues[i] = 0;	
	}
	m_iNumValuesNotSet = m_pSurface->GetMesh()->NVertices();
	
	for (int i=0; i<knownSet->NumSamples(); i++)
	{
		bool atVertex;
		std::map<int, double> & nhd = tempDist->GetNeighborhood(knownSet->GetSample(i).NearestVertex(&atVertex));
		assert(atVertex);
		for (std::map<int, double>::iterator iter = nhd.begin(); iter!=nhd.end(); iter++)
		{
			int neighborID = iter->first;
			double dist = iter->second / geodesicNorm;
			double weight = exp(-(dist*dist) / sigma2);
			normValues[neighborID] += weight;
			m_PerVertexValues[neighborID] += weight * knownValues[i];
		}
	}
	
	bool coveredAll = true;
	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
	{
		SetVertexValue(i, m_PerVertexValues[i] / normValues[i]);
		if (normValues[i]==0)
			coveredAll = false;
	}
	return coveredAll;
}

void SurfacePerVertexValue::SmoothGaussianFromKnownExact(SurfaceSampleSet * knownSet, double sigma)
{
	AllocatePerVertexStorageFillFromKnown();		
	DistanceOnTheFly * tempDist = GetSurfaceDistance();
	double geodesicNorm = sqrt(m_pSurface->Area());	
	sigma = ReestimateSigma(sigma, knownSet->NumSamples());
	std::cout<<"[WARNING] Gaussian Smoothing of signal ("<<m_pSurface->GetMesh()->NVertices()<<"): "<<std::flush;
	double sigma2 = sigma*sigma;
	double * knownValues = GetKnownValueArray(knownSet->NumSamples());
	for (int i=0; i<knownSet->NumSamples(); i++)
	{
		knownValues[i] = GetInterpolatedValue(knownSet->GetSample(i));
		bool atVertex;
		tempDist->PrecomputeRow(knownSet->GetSample(i).NearestVertex(&atVertex));
		assert(atVertex);
	}

	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
	{
		if (i%500==0)
			std::cout<<i<<" "<<std::flush;
		
		double vertexVal = 0;
		double vertexNorm=0;
		for (int j=0; j<knownSet->NumSamples(); j++)
		{
			double dist = tempDist->Distance(SurfaceSample(i, m_pSurface->GetMesh()), 
											 knownSet->GetSample(j)) / geodesicNorm;
			double weight = exp(-(dist*dist) / sigma2);
			vertexVal += knownValues[j] * weight;
			vertexNorm += weight;
		}
		SetVertexValue(i, vertexVal / vertexNorm);
	}
	std::cout<<std::endl;
	assert(AllValuesSet());
}

void SurfacePerVertexValue::SmoothPickNearestNeighborsForUnassigned()
{
	assert(false);
}

DistanceOnTheFly * SurfacePerVertexValue::GetSurfaceDistance(double maxDistance, 
															 const VKString & checkDistance)
{
	SurfaceDistance * currDistance = m_pSurface->GetDistanceMetric(checkDistance);
	if (currDistance==NULL)
		return m_pSurface->GetOnTheFlyDistanceMetric(maxDistance, "none", -1);

	return m_pSurface->GetOnTheFlyDistanceMetric(maxDistance, checkDistance, -1);
}

double SurfacePerVertexValue::ReestimateSigma(double sigma, int numSamples)
{
	assert(numSamples>=0);
	if (sigma < 0)
		return 0.1;
	else
		return sigma;
}

bool SurfacePerVertexValue::GetPerVertexValues(std::vector<double> & saveValues, 
											   double setMin, double setMax, bool maxToOne)
{
	if (!AllValuesSet())
	{
		std::cout<<"[WARNING] Not all vertices were defined: Undefined="<<m_iNumValuesNotSet<<" / ";
		std::cout<<m_pSurface->GetMesh()->NVertices()<<std::endl;
		if (m_iNumValuesNotSet==m_pSurface->GetMesh()->NVertices())
		{
			for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
				saveValues.push_back(-1);
			return false;
		}
		SmoothGaussianFromKnown(-1);
	}
	double maxVal=-FLT_MAX;
	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
	{
		double val = m_PerVertexValues[i];
		if (val < setMin)
			val = setMin;
		if (val > setMax)
			val = setMax;
		if (maxToOne)
			maxVal = vkMax(maxVal, val);
		saveValues.push_back(val);
	}
	
	if (maxToOne)
	{
		for (int i=0; i<(int)saveValues.size(); i++)
			saveValues[i] /= maxVal;
	}
	return true;
}

void SurfacePerVertexValue::LocalExtrema(SurfaceSampleSet * minSet, 
										 SurfaceSampleSet * maxSet, double minGeodesic)
{		
	assert(AllValuesSet());
	std::vector<Sortable> minsSorted;
	std::vector<Sortable> maxsSorted;
	R3Mesh * mesh = m_pSurface->GetMesh();
	
	for (int i=0; i<mesh->NVertices(); i++)
		mesh->SetVertexValue(mesh->Vertex(i), m_PerVertexValues[i]);
	
	VertexData * vertexTempArray = new VertexData[mesh->NVertices()];
	for (int v=0; v<mesh->NVertices(); v++)
	{
		vertexTempArray[v].vertex = mesh->Vertex(v);
		vertexTempArray[v].distance = 0;
		mesh->SetVertexData(mesh->Vertex(v), &(vertexTempArray[v]));
	}
	
	for (int v=0; v<mesh->NVertices(); v++)
	{
		R3MeshVertex * vp = mesh->Vertex(v);
		bool minValue = true;
		bool maxValue = true;
		// faster: check if local maxima
		IsExtrema(mesh, vertexTempArray, minValue, maxValue, vp, minGeodesic);
		
		if (minValue)	
		{
//			std::cout<<"Adding min vertex "<<v<<std::endl;
			minsSorted.push_back(Sortable(mesh->VertexValue(vp), NULL, v));
		}
		if (maxValue)
		{
//			std::cout<<"Adding max vertex "<<v<<std::endl;			
			maxsSorted.push_back(Sortable(mesh->VertexValue(vp), NULL, v));
		}
	}
	
	for (int v=0; v<mesh->NVertices(); v++)
		mesh->SetVertexData(mesh->Vertex(v), NULL);
	delete [] vertexTempArray;
	
	std::sort(minsSorted.begin(), minsSorted.end());
	std::sort(maxsSorted.begin(), maxsSorted.end());	

	for (int i=0; i<(int)minsSorted.size(); i++)
		minSet->AddSample(SurfaceSample(minsSorted[i].id, mesh));
	for (int i=(int)maxsSorted.size()-1; i>=0; i--)
		maxSet->AddSample(SurfaceSample(maxsSorted[i].id, mesh));
}

void SurfacePerVertexValue::CutToMaximum(double maxValue)
{
	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
		if (m_PerVertexValues[i] > maxValue)
			m_PerVertexValues[i] = maxValue;
}

void SurfacePerVertexValue::CutToMinimum(double minValue)
{
	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
		if (m_PerVertexValues[i] > minValue)
			m_PerVertexValues[i] = minValue;	
}

double SurfacePerVertexValue::MinValue()
{
	return m_fMinValue;
}

double SurfacePerVertexValue::MaxValue()
{
	return m_fMaxValue;
}

void SurfacePerVertexValue::WritePerVertexValuesInArff(const VKString & filename,
													   const VKString & relation, 
													   const VKString & attributeName)
{
	std::ofstream textStream(filename.c_str());
	assert(textStream.is_open());
	
	assert(AllValuesSet());
	
	textStream<<"@relation "<<relation.c_str()<<"\n\n";
	textStream<<"@attribute "<<attributeName.c_str()<<"\n";	
	textStream<<"@data"<<"\n";
	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
		textStream<<GetInterpolatedValue(SurfaceSample(i, m_pSurface->GetMesh()))<<"\n";
}

void SurfacePerVertexValue::WritePerVertexValues(std::ofstream & textStream)
{
	if (AllValuesSet())
	{
		textStream<<"AllValues\n";
		textStream<<m_pSurface->GetMesh()->NVertices()<<"\n";
		for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
			textStream<<m_PerVertexValues[i]<<" ";
		textStream<<"\n";
	}
	else
	{
		textStream<<"SomeValues\n";
		int valsSet = m_pSurface->GetMesh()->NVertices()-m_iNumValuesNotSet; 
		if ((int)m_VertexToValueMap.size() != valsSet)
		{
			std::cout<<"[ERROR] Number of values set is inconsistent: "<<valsSet<<" / ";
			std::cout<<m_pSurface->GetMesh()->NVertices()<<". Not set = "<<m_iNumValuesNotSet<<std::endl;
		}
		assert((int)m_VertexToValueMap.size() == valsSet);
		textStream<<valsSet<<"\n";
		for (std::map<int, double>::iterator iter = m_VertexToValueMap.begin();
			 iter != m_VertexToValueMap.end(); iter++)
			textStream<<iter->first<<" "<<iter->second<<" ";
		textStream<<"\n";
		
//		
//		int setVals=0;
//		for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
//		{
//			if (m_PerVertexValuesExist[i])
//			{
//				textStream<<i<<" "<<m_PerVertexValues[i]<<" ";
//				setVals++;
//			}
//		}
//		textStream<<"\n";
//		assert(setVals==(m_pSurface->GetMesh()->NVertices()-m_iNumValuesNotSet));		
	}
}

bool SurfacePerVertexValue::LoadPerVertexValues(std::ifstream & textStream)
{
	std::string tempStr;
	textStream>>tempStr;
	double val;	
	int numVals=0;	

	m_bAllVertexValuesAllocated = (VKString(tempStr.c_str())=="AllValues");
	
	if (m_bAllVertexValuesAllocated)
	{
		AllocatePerVertexStorage();
		textStream>>numVals;
		assert(numVals==m_pSurface->GetMesh()->NVertices());
		for (int i=0; i<numVals; i++)
		{
			textStream>>val;
			SetVertexValue(i, val);
		}
		return true;
	}
	else 
	{
		assert(VKString(tempStr.c_str())=="SomeValues");
		textStream>>numVals;
		int id;
		for (int i=0; i<numVals; i++)
		{
			textStream>>id;
			textStream>>val;
			SetVertexValue(id, val);
		}
		return true;
	}
	assert(false);
	return false;
}

bool SurfacePerVertexValue::AllValuesSet()
{
	assert(m_iNumValuesNotSet>=0);
	return (m_iNumValuesNotSet==0 && m_bAllVertexValuesAllocated);
}

double SurfacePerVertexValue::GetPatchArea(R3Mesh * mesh, int vertexID)
{
	return GetPatchArea(mesh, mesh->Vertex(vertexID));
}

double SurfacePerVertexValue::GetPatchArea(R3Mesh * mesh, R3MeshVertex * vertex)
{
	double area =0;
	for (int i=0; i<mesh->VertexValence(vertex); i++)
	{
		R3MeshFace * face = mesh->FaceOnVertex(vertex, mesh->EdgeOnVertex(vertex, i));
		if (face==NULL)
			continue;	
		area += mesh->FaceArea(face);
	}
	return area / 3.;
}


void SurfacePerVertexValue::IsExtrema(R3Mesh * mesh, VertexData * vertex_data, bool & isMin, bool &isMax, 
									  R3MeshVertex * source_vertex, double ringRadius)
{
	double myValue = mesh->VertexValue(source_vertex);
	isMin = true;
	isMax = true;
	
	// Re-set vertex data
	for (int j = 0; j < mesh->NVertices(); j++) 
	{
		vertex_data[j].distance = FLT_MAX;
		vertex_data[j].heappointer = NULL;
	}
	
	// Initialize priority queue
	VertexData *data = (VertexData *) mesh->VertexData(source_vertex);
	assert(data!=NULL);
	data->distance = 0;
	RNHeap<VertexData *> heap(data, &(data->distance), &(data->heappointer));
	heap.Push(data);		
	
	// Visit other nodes computing shortest distance
	while (!heap.IsEmpty()) 
	{
		VertexData *data = heap.Pop();
		R3MeshVertex *vertex = data->vertex;
		if (vertex!=source_vertex)
		{
			if (mesh->VertexValue(vertex) >= myValue)
				isMax = false;
			else if (mesh->VertexValue(vertex) == myValue)
				isMax = mesh->VertexID(vertex) < mesh->VertexID(source_vertex);

			if (mesh->VertexValue(vertex) <= myValue)
				isMin = false;
			else if (mesh->VertexValue(vertex) == myValue)
				isMin = mesh->VertexID(vertex) < mesh->VertexID(source_vertex);
						
			if (!isMin && !isMax)
				return;
		}
		
		for (int j = 0; j < mesh->VertexValence(vertex); j++) 
		{
			R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, j);
			R3MeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
			VertexData *neighbor_data = (VertexData *) mesh->VertexData(neighbor_vertex);
			RNScalar old_distance = neighbor_data->distance;
			
			RNScalar new_distance = mesh->EdgeLength(edge) + data->distance;
			assert(new_distance>=0);
			if (new_distance < old_distance 
				 && (new_distance < ringRadius 
					 || mesh->EdgeBetweenVertices(source_vertex, vertex)!=NULL)) 
			{
				neighbor_data->distance = new_distance;
				if (old_distance < FLT_MAX) heap.Update(neighbor_data);
				else heap.Push(neighbor_data);
			}
		}
	}	
}


double * SurfacePerVertexValue::GetKnownValueArray(int id)
{
	if (s_KnownValuesArrays.find(id)==s_KnownValuesArrays.end())
		s_KnownValuesArrays[id] = new double[id];
	return s_KnownValuesArrays[id];
}

double * SurfacePerVertexValue::GetNormalizationArray(int id)
{
	if (s_NormalizationArrays.find(id)==s_NormalizationArrays.end())
		s_NormalizationArrays[id] = new double[id];
	return s_NormalizationArrays[id];
}

/// static functions
void SurfacePerVertexValue::PrecomputedWeights::GetSampleWeighting(int vID, double ** savedWeights, 
																   int ** savedIDs, int * savedNumPnts)
{
	*savedWeights = m_WeightArray[vID];
	*savedIDs = m_SampleIDs[vID];
	*savedNumPnts = m_NumPointsInNhd[vID];
}

SurfacePerVertexValue::PrecomputedWeights * 
SurfacePerVertexValue::GetInterpolationWeighting(const VKString & sampleSetName,
												 SampledSurface * surf, double sigma)
{
	std::map<SampledSurface *, std::vector<PrecomputedWeights *> >::iterator iter;
	iter = s_InterpolationWeights.find(surf);
	
	if (iter!=s_InterpolationWeights.end())
	{
		for (int i=0; i<(int)iter->second.size(); i++)
			if (iter->second[i]->m_fSigma==sigma && iter->second[i]->m_sSampleSet==sampleSetName)
				return iter->second[i];
	}
	
	PrecomputedWeights * weights = new PrecomputedWeights();
	s_InterpolationWeights[surf].push_back(weights);
	
	weights->m_sSampleSet = sampleSetName;
	weights->m_fSigma = sigma;
	weights->m_pSurface = surf;
	
	weights->m_WeightArray = new double *[surf->GetMesh()->NVertices()];
	weights->m_SampleIDs = new int *[surf->GetMesh()->NVertices()];	
	weights->m_NumPointsInNhd = new int[surf->GetMesh()->NVertices()];

	SurfaceDistance * distanceMetric = surf->GetDistanceMetric("default");
	assert(distanceMetric!=NULL);
	SurfaceSampleSet * sampleSet = surf->GetSampleSet(sampleSetName);	
	assert(sampleSet!=NULL);
	double weightThreshold = 1. / (100. * (double)sampleSet->NumSamples());	
	double sigma2 = sigma * sigma;
	
	int minNumSamp=INT_MAX;
	int maxNumSamp=INT_MIN;
	int aveNumSamp=0;
	double geodesicNorm = sqrt(surf->Area());
		// for each vertex - find weights
	for (int i=0; i<surf->GetMesh()->NVertices(); i++)
	{
//		AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("PrepareSmoothing", 
//															i, surf->GetMesh()->NVertices());
		SurfaceSample vSamp(i, surf->GetMesh());
		// find unnormalized weights 
		double norm = 0;
		std::vector<double> unnormalizedWeights;
		for (int j=0; j<sampleSet->NumSamples(); j++)
		{
			double dist = distanceMetric->Distance(vSamp, sampleSet->GetSample(j)) / geodesicNorm;
			double currWeight = exp(-(dist*dist) / sigma2 );
			unnormalizedWeights.push_back(currWeight);
			norm += currWeight;
			//std::cout<<"\tWeight["<<j<<"] = "<<currWeight<<" = exp(-"<<dist<<"^2/"<<sigma<<"^2)"<<std::endl;
		}
		//std::cout<<"norm="<<norm<<std::endl;
		
		// normalize weights & prune based on a threshold
		double normAfterPruning=0;
		std::vector<double> normalizedWeights;	
		std::vector<int> acceptedWeightsID;
		for (int j=0; j<(int)unnormalizedWeights.size(); j++)
		{
			double currWeight = unnormalizedWeights[j] / norm;
			//std::cout<<"currweight ["<<currWeight<<"] "<<std::endl;
			if (currWeight >= weightThreshold)
			{
				acceptedWeightsID.push_back(j);
				normalizedWeights.push_back(currWeight);
				normAfterPruning += currWeight;
			}
		}
		
		// store normalized weights
		int N = (int)normalizedWeights.size();
		assert(N>0);		
		weights->m_NumPointsInNhd[i] = N;
		weights->m_WeightArray[i] = new double[N];
		weights->m_SampleIDs[i] = new int[N];
		minNumSamp=vkMin(minNumSamp, N);
		maxNumSamp=vkMax(maxNumSamp, N);		
		aveNumSamp+=N;
		
		for (int j=0; j<N; j++)
		{
			weights->m_SampleIDs[i][j] = acceptedWeightsID[j];
			weights->m_WeightArray[i][j] = (normalizedWeights[j] / normAfterPruning);
		}
	}
//	std::cout<<"Average Num Samples: "<<(double)aveNumSamp / (double)surf->GetMesh()->NVertices();
//	std::cout<<"\nRange Num Samples: ["<<minNumSamp<<", "<<maxNumSamp<<"]"<<std::endl;

	return weights;
}

void SurfacePerVertexValue::ClearPrecomputedWeights()
{
	std::map<SampledSurface *, std::vector<PrecomputedWeights *> >::iterator iter;
	for (iter = s_InterpolationWeights.begin(); iter != s_InterpolationWeights.end(); iter++)
	{
		for (int i=0; i<(int)iter->second.size(); i++)
		{
			for (int v=0; v<iter->first->GetMesh()->NVertices(); v++)
			{
				delete [] iter->second[i]->m_WeightArray[v];
				delete [] iter->second[i]->m_SampleIDs[v];
			}
			delete [] iter->second[i]->m_WeightArray;
			delete [] iter->second[i]->m_SampleIDs;
			delete [] iter->second[i]->m_NumPointsInNhd;
		}
	}
}


