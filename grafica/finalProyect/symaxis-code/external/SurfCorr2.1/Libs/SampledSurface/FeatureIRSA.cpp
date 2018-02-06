#include "FeatureIRSA.h"
#include "MapScoredCollection.h"
#include "SampledSurface.h"
#include "SurfaceMapConfidence.h"
#include "MapFlattening.h"
#include "DistanceOnTheFly.h"

FeatureIRSA::FeatureIRSA(SampledSurface * surface, const VKString & mapCollection,
						 const VKString & confidenceName, const VKString & distanceMetric,
						 const VKString & featureName)
: SurfaceFeature(surface, true, featureName)
{
	m_MapCollectionName = mapCollection;
	m_ConfidenceName = confidenceName;
	m_DistanceMetric = distanceMetric;
}

FeatureIRSA::~FeatureIRSA()
{
}

VKString FeatureIRSA::GetFeatureName()
{
	return "FeatureIRSA";
}

void FeatureIRSA::PrecomputeVertexFeatures()
{
	std::cout<<"Precomputing Intrinsic Reflective Symmetry Axis"<<std::endl;
	SurfaceMap * surfMap = m_pSurface->GetMapActingAsFrom(m_MapCollectionName);
	assert(surfMap!=NULL);
	assert(surfMap->GetSurfaceMapType()=="MapScoredCollection");
	MapScoredCollection * mapCollection = (MapScoredCollection*)surfMap;
	int N = mapCollection->GetNumMaps();
	int NVert = m_pSurface->GetMesh()->NVertices();
	double * values = new double[NVert];
	for (int i=0; i<NVert; i++)
		values[i] = 0;
	
	double distNorm = sqrt(m_pSurface->Area());
	//SurfaceDistance * dist = m_pSurface->GetDistanceMetric(m_DistanceMetric);
	DistanceOnTheFly * dist = m_pSurface->GetOnTheFlyDistanceMetric(-1, m_DistanceMetric, -1);
	assert(dist!=NULL);
//	DistanceOnTheFly * distShort = m_pSurface->GetOnTheFlyDistanceMetric(.25, m_DistanceMetric, -1);
//	assert(distShort!=NULL);
	
//	SurfacePerVertexValue amIStationary(m_pSurface, false);
//	SurfaceSampleSet 
	
	for (int i=0; i<N; i++)
	{
		if (i%10==0)
			AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("CalculatingPIRS", i, N);
		
		double score;
		SurfaceMap * currMap = mapCollection->GetMapByID(i, &score);
		assert(currMap->GetSurface(0)==currMap->GetSurface(1));	// symmetry
		
		// TODO; can have arbitrary check if defined by coarse set
		assert(currMap->GetSurfaceMapType()=="MapVia2DPlane");	
		
		MapVia2DPlane * confMap = (MapVia2DPlane*)currMap;
		
		std::vector<int> * corrs;
		SurfaceSampleSet * genSet1;
		SurfaceSampleSet * genSet2;
		confMap->GetGeneratingSet(&corrs, &genSet1, &genSet2);
		assert(genSet1==genSet2);
		
		assert(corrs!=NULL && genSet1!=NULL && genSet2!=NULL);
		assert(corrs->size()==6 && (*corrs)[0]==(*corrs)[3] 
			   && (*corrs)[1]==(*corrs)[2] && (*corrs)[4]==(*corrs)[5]);
		
		SurfaceSample s1 = genSet1->GetSample((*corrs)[0]);
		SurfaceSample s2 = genSet1->GetSample((*corrs)[1]);
		SurfaceSample s3 = genSet1->GetSample((*corrs)[4]);
		
		dist->PrecomputeRow(s1.NearestVertex());
		dist->PrecomputeRow(s2.NearestVertex());
		
		// OK params: 0.01, 0.1, 0.05 (conservative) - if stationary is the furthest
		// OK params: 0.01, 0.2, 0.2 (thicker line) - if stationary is the furthest + e.g. model 1.
		// 0.01, 0.3, 0.05, 0.5 - if stationary is in the middle of geodesic curve
		double sigmaDiffDist = 0.01;
		double sigmaConf = 0.1;
		double sigmaTotalConf = 0.1;
		//double sigmaDist = 0.1;

		SurfaceMapConfidence * confidence = currMap->GetMapConfidenceCalculator(m_ConfidenceName);
		double tC = confidence->Confidence();
		tC = exp((tC-1.) / sigmaTotalConf);
		
		//double pD = exp((dist->Distance(s1, s2) / distNorm) / sigmaDist);
		double pD = 1.;
		
		for (int v=0; v<NVert; v++)	// for every vertex cast a vote
		{
			SurfaceSample samp(v, m_pSurface->GetMesh());
			double d1 = dist->Distance(s1, samp) / distNorm;
			double d2 = dist->Distance(s2, samp) / distNorm;
			
			//SurfaceSample sampForw = currMap->ForwardMap(samp);
			//double dstat = dist->Distance(samp, sampForw);
			
			double c = confidence->ConfidenceAtSample(samp);
			//double c = confidence->ConfidenceAtSampleNoSmoothing(samp);
			
			c = exp((c - 1.) / sigmaConf);
			
			double d = exp(-vkAbs(d1-d2) / sigmaDiffDist);
			//double d = vkMin(d1 / (d2+0.00001), d2 / (d1+0.00001));
			//double d = exp(-dstat / .01);
			
			values[v] += d * c * tC * pD;
		}
	}
	for (int i=0; i<NVert; i++)
		m_pPerVertexValues->SetVertexValue(i, values[i]);
	delete [] values;
		
	AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("CalculatingPIRS", N, N);
	std::cout<<std::endl;

}

double FeatureIRSA::CalculateValue(const SurfaceSample & sample) const
{
	assert(false);
}



