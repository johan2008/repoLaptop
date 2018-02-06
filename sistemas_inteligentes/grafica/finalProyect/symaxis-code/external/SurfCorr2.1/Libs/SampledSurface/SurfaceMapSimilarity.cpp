#include "SurfaceMapSimilarity.h"

std::map<VKString, SurfaceMapSimilarity*> SurfaceMapSimilarity::m_mNameToSeed;

SurfaceMapSimilarity::SurfaceMapSimilarity(SurfaceMap * M1, SurfaceMap * M2,
										   bool smoothSimilarity, bool smoothDissimilarity,
										   double smoothingSigma,
										   const VKString & sampleSet1, const VKString & sampleSet2,
										   const VKString & distance1, const VKString & distance2,
										   double simThreshold, bool cacheSim, bool cacheDissim)
{
	assert((M1==NULL && M2==NULL) || (M1!=NULL && M2!=NULL));
	
	if (M1!=NULL && M2!=NULL)
		assert(M1->GetSurface(0)==M2->GetSurface(0) && M1->GetSurface(1)==M2->GetSurface(1));
	
	m_bCacheSimilarityValues = cacheSim;
	m_bCacheDissimilarityValues = cacheDissim;

	m_bSmoothSimilarity = smoothSimilarity;
	m_bSmoothDissimilarity = smoothDissimilarity;
	
	m_fSmoothSigma = smoothingSigma;
	if (smoothSimilarity)
		assert(m_bCacheSimilarityValues);
	if (smoothDissimilarity)
		assert(m_bCacheDissimilarityValues);
	
	m_fSimilarityThreshold = simThreshold;
	m_pSurfaceMap1 = M1;
	m_pSurfaceMap2 = M2;	
	
	m_CoarseSet1Name = sampleSet1;
	m_CoarseSet2Name = sampleSet2;
	m_Distance1Name = distance1;
	m_Distance2Name = distance2;
	
	m_fSimilarityValue = -1;
	m_fDissimlarityValue = -1;

	if (m_bCacheSimilarityValues && M1!=NULL)
	{
		for (int i=0; i<2; i++)
		{
			SurfacePerVertexValue * values = new SurfacePerVertexValue(M1->GetSurface(i), true);
			m_mSimilarityCache[m_pSurfaceMap1->GetSurface(i)] = values;
		}
	}

	if (m_bCacheDissimilarityValues && M1!=NULL)
	{
		for (int i=0; i<2; i++)
		{
			SurfacePerVertexValue * values = new SurfacePerVertexValue(M1->GetSurface(i), true);
			m_mDissimilarityCache[m_pSurfaceMap1->GetSurface(i)] = values;
		}
	}
}

SurfaceMapSimilarity::~SurfaceMapSimilarity()
{
	if (m_bCacheSimilarityValues && m_pSurfaceMap1!=NULL)
	{
		delete m_mSimilarityCache[m_pSurfaceMap1->GetSurface(0)];
		delete m_mSimilarityCache[m_pSurfaceMap1->GetSurface(1)];		
	}

	if (m_bCacheDissimilarityValues && m_pSurfaceMap1!=NULL)
	{
		delete m_mDissimilarityCache[m_pSurfaceMap1->GetSurface(0)];
		delete m_mDissimilarityCache[m_pSurfaceMap1->GetSurface(1)];		
	}
}

SurfaceMapSimilarity * SurfaceMapSimilarity::GetClonableSeed(const VKString & seedName)
{
	if (m_mNameToSeed.find(seedName) != m_mNameToSeed.end())
		return m_mNameToSeed[seedName];
	return NULL;
}

void SurfaceMapSimilarity::SetClonableSeed(const VKString & seedName, SurfaceMapSimilarity * seed)
{
	m_mNameToSeed[seedName] = seed;
}


double SurfaceMapSimilarity::Similarity()
{
	if (m_fSimilarityValue==-1.)
		m_fSimilarityValue = CalculateTotalSimilarity();
	return m_fSimilarityValue;
}

double SurfaceMapSimilarity::Dissimilarity()
{
	if (m_fDissimlarityValue==-1.)
		m_fDissimlarityValue = CalculateTotalDissimilarity();
	return m_fDissimlarityValue;
}

double SurfaceMapSimilarity::SimilarityAtSample(const SurfaceSample & sample)
{
	if (m_bSmoothSimilarity)
		PrecomputeSmoothSimilarities();
	
	bool atVertex;
	sample.NearestVertex(&atVertex);
	SampledSurface * surface = m_pSurfaceMap1->GetSurface(sample);
	R3Mesh * mesh = surface->GetMesh();
	
	if (!m_bCacheSimilarityValues)
	{
		if (atVertex)
			return CalculateSimilarityAtSample(sample);
		else
			return sample.Interpolate(CalculateSimilarityAtSample(SurfaceSample(sample.VertID(0), mesh)),
									  CalculateSimilarityAtSample(SurfaceSample(sample.VertID(1), mesh)),
									  CalculateSimilarityAtSample(SurfaceSample(sample.VertID(2), mesh)));
	}
	
	
	SurfaceSample getValues = m_mSimilarityCache[surface]->NextAtVertexSampleRequiredForInterpolation(sample);
	while (!getValues.Invalid())
	{
		m_mSimilarityCache[surface]->SetVertexValue(getValues,
													CalculateSimilarityAtSample(getValues));
		getValues = m_mSimilarityCache[surface]->NextAtVertexSampleRequiredForInterpolation(sample);
	}
	
	return m_mSimilarityCache[surface]->GetInterpolatedValue(sample);
}

double SurfaceMapSimilarity::DissimilarityAtSample(const SurfaceSample & sample)
{
	if (m_bSmoothDissimilarity)
		PrecomputeSmoothDissimilarities();
	
	bool atVertex;
	sample.NearestVertex(&atVertex);
	SampledSurface * surface = m_pSurfaceMap1->GetSurface(sample);
	R3Mesh * mesh = surface->GetMesh();
	
//	std::cout<<"Dissimilarity At Sample: "<<sample.NearestVertex()<<" Cache="<<m_bCacheDissimilarityValues<<std::endl;
	
	if (!m_bCacheDissimilarityValues)
	{
		if (atVertex)
			return CalculateDissimilarityAtSample(sample);
		else
			return sample.Interpolate(CalculateDissimilarityAtSample(SurfaceSample(sample.VertID(0), mesh)),
									  CalculateDissimilarityAtSample(SurfaceSample(sample.VertID(1), mesh)),
									  CalculateDissimilarityAtSample(SurfaceSample(sample.VertID(2), mesh)));
	}
	
	
	SurfaceSample getValues = m_mDissimilarityCache[surface]->NextAtVertexSampleRequiredForInterpolation(sample);
	while (!getValues.Invalid())
	{
		m_mDissimilarityCache[surface]->SetVertexValue(getValues,
													   CalculateDissimilarityAtSample(getValues));
		getValues = m_mDissimilarityCache[surface]->NextAtVertexSampleRequiredForInterpolation(sample);
	}
	
	return m_mDissimilarityCache[surface]->GetInterpolatedValue(sample);
}

double SurfaceMapSimilarity::GetSimilarityThreshold(SampledSurface * , SampledSurface * )
{
	return m_fSimilarityThreshold;
}

void SurfaceMapSimilarity::PrecomputeSmoothSimilarities()
{
	PrecomputeSmoothSimilarities(m_pSurfaceMap1->GetSurface(0));
	PrecomputeSmoothSimilarities(m_pSurfaceMap1->GetSurface(1));		
}

void SurfaceMapSimilarity::PrecomputeSmoothSimilarities(SampledSurface * surface)
{
	VKString sampleSetName;
	SurfaceSampleSet * sampleSet = NULL;
	if (surface==m_pSurfaceMap1->GetSurface(0))
		sampleSetName = m_CoarseSet1Name;
	else
		sampleSetName = m_CoarseSet2Name;
	
	assert(sampleSetName!="none" && sampleSetName!="AllVertices");
	
	sampleSet = surface->GetSampleSet(sampleSetName);
	assert(sampleSet!=NULL);
	
//	std::cout<<"[WARNING] Precomputing smooth confidence ("<<sampleSet->NumSamples()<<" / ";
//	std::cout<<surface->GetMesh()->NVertices()<<") : "<<std::flush;
	
	std::vector<double> samplesInSet;
	for (int i=0; i<sampleSet->NumSamples(); i++)
	{
//		if (i%100==0)
//			std::cout<<i<<" "<<std::flush;
		
		const SurfaceSample & currSamp = sampleSet->GetSample(i);
		m_mSimilarityCache[surface]->SetVertexValue(currSamp, CalculateSimilarityAtSample(currSamp));
	}
	
	m_mSimilarityCache[surface]->SmoothGaussianFromKnown(sampleSetName, m_fSmoothSigma);
}

void SurfaceMapSimilarity::PrecomputeSmoothDissimilarities()
{
	m_bSmoothDissimilarity = false;
	
	PrecomputeSmoothDissimilarities(m_pSurfaceMap1->GetSurface(0));
	PrecomputeSmoothDissimilarities(m_pSurfaceMap1->GetSurface(1));
}

void SurfaceMapSimilarity::PrecomputeSmoothDissimilarities(SampledSurface * surface)
{
	VKString sampleSetName;
	SurfaceSampleSet * sampleSet = NULL;
	if (surface==m_pSurfaceMap1->GetSurface(0))
		sampleSetName = m_CoarseSet1Name;
	else
		sampleSetName = m_CoarseSet2Name;
	
	assert(sampleSetName!="none" && sampleSetName!="AllVertices");
	
	sampleSet = surface->GetSampleSet(sampleSetName);
	assert(sampleSet!=NULL);
	
//	std::cout<<"[WARNING] Precomputing smooth similarity ("<<sampleSet->NumSamples()<<" / ";
//	std::cout<<surface->GetMesh()->NVertices()<<") : "<<std::flush;
	
	std::vector<double> samplesInSet;
	for (int i=0; i<sampleSet->NumSamples(); i++)
	{
//		if (i%100==0)
//			std::cout<<i<<" "<<std::flush;
		
		const SurfaceSample & currSamp = sampleSet->GetSample(i);
		m_mDissimilarityCache[surface]->SetVertexValue(currSamp, CalculateDissimilarityAtSample(currSamp));
	}
	
	m_mDissimilarityCache[surface]->SmoothGaussianFromKnown(sampleSetName, m_fSmoothSigma);
}

double SurfaceMapSimilarity::CalculateTotalSimilarity()
{
	bool disableSmoothing = m_bSmoothSimilarity;
	m_bSmoothSimilarity = false;
	int finalNormalization = 0;
	double totalSimilarity = 0;
	if (m_CoarseSet1Name!="none")
	{
		double locTotal=0;
		SurfaceSampleSet * sampleSet = m_pSurfaceMap1->GetSurface(0)->GetSampleSet(m_CoarseSet1Name);
		assert(sampleSet!=NULL);
		for (int i=0; i<sampleSet->NumSamples(); i++)
			locTotal += SimilarityAtSample(sampleSet->GetSample(i));
		totalSimilarity += locTotal / sampleSet->NumSamples();
		finalNormalization++;	
	}

	if (m_CoarseSet2Name!="none")
	{
		double locTotal=0;
		SurfaceSampleSet * sampleSet = m_pSurfaceMap1->GetSurface(1)->GetSampleSet(m_CoarseSet2Name);
		assert(sampleSet!=NULL);
		for (int i=0; i<sampleSet->NumSamples(); i++)
			locTotal += SimilarityAtSample(sampleSet->GetSample(i));
		totalSimilarity += locTotal / sampleSet->NumSamples();
		finalNormalization++;	
	}
	m_bSmoothSimilarity = disableSmoothing;
	return totalSimilarity/finalNormalization;
}

double SurfaceMapSimilarity::CalculateTotalDissimilarity()
{
	bool disableSmoothing = m_bSmoothDissimilarity;
	m_bSmoothDissimilarity = false;
	
	int finalNormalization = 0;
	double totalVal = 0;
	if (m_CoarseSet1Name!="none")
	{
		double locTotal=0;
		SurfaceSampleSet * sampleSet = m_pSurfaceMap1->GetSurface(0)->GetSampleSet(m_CoarseSet1Name);
		assert(sampleSet!=NULL);
		for (int i=0; i<sampleSet->NumSamples(); i++)
			locTotal += DissimilarityAtSample(sampleSet->GetSample(i));
		totalVal += locTotal / sampleSet->NumSamples();
		finalNormalization++;	
	}
	
	if (m_CoarseSet2Name!="none")
	{
		double locTotal=0;
		SurfaceSampleSet * sampleSet = m_pSurfaceMap1->GetSurface(1)->GetSampleSet(m_CoarseSet2Name);
		assert(sampleSet!=NULL);
		for (int i=0; i<sampleSet->NumSamples(); i++)
			locTotal += DissimilarityAtSample(sampleSet->GetSample(i));
		totalVal += locTotal / sampleSet->NumSamples();
		finalNormalization++;	
	}
	m_bSmoothDissimilarity = disableSmoothing;
	return totalVal/finalNormalization;
}

double SurfaceMapSimilarity::CalculateDissimilarityAtSample(const SurfaceSample & samp)
{
	return 1. - SimilarityAtSample(samp);
}

void SurfaceMapSimilarity::GetPerVertexSimilarity(SampledSurface * surf, 
												  std::vector<double> & vals, bool normalize)
{
	std::cout<<"Similarity = "<<Similarity()<<std::endl;
	assert(m_mSimilarityCache.find(surf)!=m_mSimilarityCache.end());
	m_mSimilarityCache[surf]->GetPerVertexValues(vals, normalize);
}

void SurfaceMapSimilarity::GetPerVertexDissimilarity(SampledSurface * surf, 
													 std::vector<double> & vals, bool normalize)
{
	std::cout<<"\nDissimilarity = "<<Dissimilarity()<<std::endl;
	assert(m_mDissimilarityCache.find(surf)!=m_mDissimilarityCache.end());
	m_mDissimilarityCache[surf]->GetPerVertexValues(vals, normalize);
}



