#include "FeatureAGD.h"
#include "FeatureMGD.h"
#include "AnalysisStats.h"

#include "SampledSurface.h"
#include "SurfaceDistance.h"
#include "DistanceGeodesic.h"
#include "DistanceOnTheFly.h"

/////////////////// AGD - Slow calculation ///////////////
FeatureAGD::FeatureAGD(SampledSurface * surface, const VKString & distanceName, 
					bool vertexOnly, const VKString & featurename)
: SurfaceFeature(surface, vertexOnly, featurename)
{
	m_pDistance = surface->GetDistanceMetric(distanceName);
	m_fNorm = surface->Area();
	m_aVerticesDistances = NULL;
	m_aVertexData = NULL;
}

double FeatureAGD::CalculateValue(const SurfaceSample & sample) const
{
	if (m_pDistance!=NULL)
	{
		double value=0;
		double norm=0;
		for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
		{	
			double weight=SurfacePerVertexValue::GetPatchArea( m_pSurface->GetMesh(),  i) / m_fNorm;
			double dist=m_pDistance->Distance(sample, SurfaceSample(i, m_pSurface->GetMesh()));
			value+=dist*weight;
			norm+=weight;
		}
		assert(norm>0);
		return value / norm;
	}
	else
	{
		if (m_aVerticesDistances==NULL)
		{
			std::cout<<"\n[WARNING] FeatureAGD: no distances: USE FeatureFastAGD!!!!   "<<std::flush;
			assert(m_aVertexData==NULL);
				// hack since it should behave as constant, but we are caching some data
			((FeatureAGD*)this)->m_aVertexData = DistanceGeodesic::InitializeFunkPrecomputation(m_pSurface->GetMesh(), 
																								false);	
			((FeatureAGD*)this)->m_aVerticesDistances = new double[m_pSurface->GetMesh()->NVertices()];
		}
		
		double value=0;
		double norm=0;
		bool ok;
		DistanceGeodesic::PrecomputeFillFunkhouser((DistanceGeodesic::GeoVertexData*) m_aVertexData, 
												   m_pSurface->GetMesh(), sample.NearestVertex(&ok), 
												   m_aVerticesDistances);
		assert(ok);
		for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
		{	
			double weight = SurfacePerVertexValue::GetPatchArea(m_pSurface->GetMesh(), i) / m_fNorm;
			double dist = m_aVerticesDistances[i];
			value += dist*weight;
			norm += weight;
		}
		assert(norm>=0);

		return value / norm;
	}
}

VKString FeatureAGD::GetFeatureName()
{
	return "FeatureAGD";
}


///////////////////// AGD - Fast Calculation ///////////////
FeatureFastAGD::FeatureFastAGD(SampledSurface * surface, int numInitialSamples, const VKString & featurename)
: SurfaceFeature(surface, true, featurename)
{
	m_iNumInitialSamples = numInitialSamples;
}

FeatureFastAGD::~FeatureFastAGD()
{
}

VKString FeatureFastAGD::GetFeatureName()
{
	return "FeatureFastAGD";
}

void FeatureFastAGD::PrecomputeVertexFeatures()
{
	assert(m_pPerVertexValues!=NULL);
	//double norm = pow(m_pSurface->Area(), 1.5);	// why 1.5????
	double norm = pow(m_pSurface->Area(), .5);
	
	DistanceOnTheFly * distMetric = m_pPerVertexValues->GetSurfaceDistance();
		
	// spread N points
	std::vector<int> evenVertices;
	FeatureMGD::EvenlySpreadVertices(m_pSurface, m_iNumInitialSamples, evenVertices);
	std::cout<<"Calculating FastAGD"<<std::endl;	
	
	int MAX_ITER = evenVertices.size()*2;
	
	for (int i=0; i<(int)evenVertices.size(); i++)
	{
		if (i%10==0)		
			AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("CreatingSurface_FastAGD", i, MAX_ITER);
		distMetric->PrecomputeRow(evenVertices[i]);
	}
	
//	std::cout<<" ("<<evenVertices.size()<<"): "<<std::flush;
	// calculate AGD value at evenly spread vertices
	for (int i=0; i<(int)evenVertices.size(); i++)
	{
		if (i%10==0)
			AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("CreatingSurface_FastAGD", (int)evenVertices.size()+i, MAX_ITER);		
//		if (i%100==0)
//			std::cout<<"."<<std::flush;
		double value=0;
		SurfaceSample fromSamp(evenVertices[i], m_pSurface->GetMesh());
		for (int j=0; j<m_pSurface->GetMesh()->NVertices(); j++)	// distance every other point
		{
			double w = SurfacePerVertexValue::GetPatchArea(m_pSurface->GetMesh(), j);
			value += distMetric->Distance(fromSamp, SurfaceSample(j, m_pSurface->GetMesh())) * w / norm;
		}
		m_pPerVertexValues->SetVertexValue(evenVertices[i], value);
	}
//	std::cout<<" - Done!"<<std::endl;
	std::cout<<"  Smoothing"<<std::flush;			
	m_pPerVertexValues->SmoothGaussianFromKnown(-1);
	std::cout<<std::endl;	
}

double FeatureFastAGD::CalculateValue(const SurfaceSample & sample) const
{
	assert(false);
}




