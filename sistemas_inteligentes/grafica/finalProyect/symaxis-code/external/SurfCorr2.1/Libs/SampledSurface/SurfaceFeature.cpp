#include "SurfaceFeature.h"
#include "SampledSurface.h"
#include "Sortable.h"
#include "AnalysisPipeline.h"
#include <algorithm>

////////////////////////// BASIC FUNCTIONALITIES ///////////////////////////////
SurfaceFeature::SurfaceFeature(SampledSurface * surface, bool perVertexFeature, VKString name)
{
	m_pSurface = surface;
	if (perVertexFeature)
		m_pPerVertexValues = new SurfacePerVertexValue(surface, false);
	else
		m_pPerVertexValues = NULL;
	
	m_sFeatureName = name;
}

void SurfaceFeature::WriteARFFPropertyFile(const VKString & filename)
{
	assert(m_pPerVertexValues!=NULL);
	
	m_pPerVertexValues->WritePerVertexValuesInArff(filename, "unknown", GetFeatureName());	
}

void SurfaceFeature::WriteFeature(const VKString & featureName)
{
	assert(m_pPerVertexValues!=NULL);
	std::ofstream dataStream(featureName.c_str(), std::ios::out);
	if(!dataStream.is_open())
	{
		std::cout<<"[ERROR] Could not write feature to file: "<<featureName.c_str()<<std::endl;
		assert(dataStream.is_open());
	}

	dataStream<<m_sFeatureName.c_str()<<"\n";
	m_pPerVertexValues->WritePerVertexValues(dataStream);
	dataStream.close();
}

bool SurfaceFeature::LoadFeature(const VKString & featureName)
{
	assert(m_pPerVertexValues!=NULL);	
	std::ifstream dataStream(featureName.c_str());	
	if (!dataStream.is_open())
		return false;
	VKString featureNameLoaded;
	featureNameLoaded.readToken(dataStream);
	if (featureNameLoaded != m_sFeatureName)
	{
		std::cout<<"[WARNING] Loading feature: File="<<featureName.c_str()<<std::endl;
		std::cout<<"\tLoaded names are different: "<<featureNameLoaded.c_str()<<" vs "<<m_sFeatureName.c_str()<<std::endl;
	}
	m_pPerVertexValues->LoadPerVertexValues(dataStream);

	dataStream.close();	
	return true;
}


double SurfaceFeature::Value(const SurfaceSample & sample)
{
	if (m_pPerVertexValues==NULL)
	{
		return CalculateValue(sample);
	}
	else
	{
		if (!m_pPerVertexValues->NextAtVertexSampleRequiredForInterpolation(sample).Invalid())
			PrecomputeVertexFeatures();	// not precomputed

		return m_pPerVertexValues->GetInterpolatedValue(sample);		
	}
}

double SurfaceFeature::GetMaximalValue()
{
	if (m_pPerVertexValues!=NULL)
		return m_pPerVertexValues->MaxValue();

	assert(false);	//otherwise child needs to implement it
	return -1;
}

double SurfaceFeature::GetMinimalValue()
{
	if (m_pPerVertexValues!=NULL)
		return m_pPerVertexValues->MinValue();
	assert(false);	//otherwise child needs to implement it
	return -1;
}

void SurfaceFeature::Extremas(SurfaceSampleSet * minSet, SurfaceSampleSet * maxSet, double maxInRing)
{
	m_pPerVertexValues->LocalExtrema(minSet, maxSet, maxInRing);
	VKString minName = GetFeatureName() + "_Min";
	VKString maxName = GetFeatureName() + "_Max";	
	minSet->SetSamplesType(minName);
	maxSet->SetSamplesType(maxName);	
}


void SurfaceFeature::PrecomputeVertexFeatures()
{
	if (AnalysisPipeline::Verb("Min"))
		std::cout<<"\tFeature: "<<m_sFeatureName.c_str()<<" Precomputing per vertex values ("<<m_pSurface->GetMesh()->NVertices()<<"):"<<std::flush;
	
	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
	{
		SurfaceSample samp(i, m_pSurface->GetMesh());
		m_pPerVertexValues->SetVertexValue(samp, CalculateValue(samp));
	}
	
	if (AnalysisPipeline::Verb("Min"))
		std::cout<<" - Done!"<<std::endl;
}


