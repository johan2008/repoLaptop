#include "MapConfidenceDistanceToGen.h"

MapConfidenceDistanceToGen::MapConfidenceDistanceToGen(SurfaceMap * mapToAnalyze, 
													   const VKString & metric1Name, 
													   const VKString & metric2Name, 
													   const VKString & set1ToVerify, 
													   const VKString & set2ToVerify)	// arbitrary set
: SurfaceMapConfidence(mapToAnalyze, false, false, -1, set1ToVerify, set2ToVerify, true, false)
{
	m_sDistName1 = metric1Name;
	m_sDistName2 = metric2Name;

	if (mapToAnalyze!=NULL)
	{
		std::vector<int> * correspondences;
		SurfaceSampleSet * samplesFrom;
		SurfaceSampleSet * samplesTo;	

		assert(mapToAnalyze->GetSurfaceMapType()=="MapConformal");
		((MapConformal*)mapToAnalyze)->GetGeneratingCorrespondences(&correspondences, &samplesFrom, &samplesTo);
		SampledSurface * surfaceFrom = mapToAnalyze->GetSurface(0);
		SampledSurface * surfaceTo = mapToAnalyze->GetSurface(1);	
		
		for (int i=0; i<(int)(*correspondences).size(); i+=2)
		{
			m_mGenerators[surfaceFrom].push_back(samplesFrom->GetSample((*correspondences)[i+0]));
			m_mGenerators[surfaceTo].push_back(samplesTo->GetSample((*correspondences)[i+1]));
		}
		m_fEpsilon = 0.1;
	}
}

MapConfidenceDistanceToGen::~MapConfidenceDistanceToGen()
{
}

VKString MapConfidenceDistanceToGen::GetMapConfidenceType()
{
	return "MapConfidenceDistanceToGen";
}

SurfaceMapConfidence * MapConfidenceDistanceToGen::Clone(SurfaceMap * currMap)
{
	return new MapConfidenceDistanceToGen(currMap, m_sDistName1, m_sDistName2, m_SampleSet1, m_SampleSet2);
}

double MapConfidenceDistanceToGen::CalculateConfidenceAtSample(const SurfaceSample & sample)
{
	//std::cout<<"MapConfidenceDistanceToGen - Calc - 1 "<<std::endl;
	SampledSurface * surface = m_pSurfaceMap->GetSurface(sample);
	if (m_mMaxConfidence.find(surface)==m_mMaxConfidence.end())
	{
		m_mMaxConfidence[surface] = -FLT_MAX;
		
		VKString sampleSetName;
		VKString distanceMetricName;
		if (surface==m_pSurfaceMap->GetSurface(0))
		{
			sampleSetName = m_SampleSet1;
			distanceMetricName = m_sDistName1;
		}
		else
		{
			sampleSetName = m_SampleSet2;
			distanceMetricName = m_sDistName2;
		}
		
		SurfaceSampleSet * sampleSet = surface->GetSampleSet(sampleSetName);
		SurfaceDistance * distance = surface->GetDistanceMetric(distanceMetricName);
		assert(sampleSet!=NULL);
		assert(distance!=NULL);
		
		double norm = sqrt(surface->Area());
		for (int i=0; i<sampleSet->NumSamples(); i++)
		{
			const SurfaceSample & samp = sampleSet->GetSample(i);
			double value=0;
			std::vector<SurfaceSample> & sampleVec = m_mGenerators[m_pSurfaceMap->GetSurface(samp)];
			for (int j=0; j<(int)sampleVec.size(); j++)
				value += 1. / (m_fEpsilon + distance->Distance(samp, sampleVec[j]) / norm);
			
			m_CachedConfidenceOnSurf[surface]->SetVertexValue(samp, value);
			m_mMaxConfidence[surface] = vkMax(m_mMaxConfidence[surface], value);
		}
		
		for (int i=0; i<sampleSet->NumSamples(); i++)
		{
			const SurfaceSample & samp = sampleSet->GetSample(i);
			double currValue = m_CachedConfidenceOnSurf[surface]->GetInterpolatedValue(samp);
			m_CachedConfidenceOnSurf[surface]->SetVertexValue(samp, currValue / m_mMaxConfidence[surface]);
		}			
	}
	//std::cout<<"MapConfidenceDistanceToGen - Calc - DONE! "<<std::endl;
	assert(m_CachedConfidenceOnSurf[surface]->NextVertexRequiredForInterpolation(sample)==-1);
	return m_CachedConfidenceOnSurf[surface]->GetInterpolatedValue(sample);
}


