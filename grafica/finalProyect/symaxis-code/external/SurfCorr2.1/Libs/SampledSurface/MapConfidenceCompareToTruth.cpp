#include "MapConfidenceCompareToTruth.h"
#include "DistanceOnTheFly.h"

MapConfidenceCompareToTruth::MapConfidenceCompareToTruth(SurfaceMap * mapToAnalyze, 
														 SurfaceMap * trueMap, 
														 const VKString & metric1Name, 
														 const VKString & metric2Name, 
														 const VKString & set1ToVerify,
														 const VKString & set2ToVerify, 
														 double unadjustedThreshold, 
														 bool bidir)	// arbitrary set
: SurfaceMapConfidence(mapToAnalyze, false, false, -1, set1ToVerify, set2ToVerify, true, true)
{
	if (mapToAnalyze!=NULL)
		assert(mapToAnalyze->GetSurface(0)==trueMap->GetSurface(0)
			   && mapToAnalyze->GetSurface(1)==trueMap->GetSurface(1));
	
	m_pTrueMap = trueMap;
	m_fUnadjustedThreshold = unadjustedThreshold;
	
	m_sDistName1 = metric1Name;
	m_sDistName2 = metric2Name;	
}

bool MapConfidenceCompareToTruth::DisplayProgressBar()
{
	return true;
}

VKString MapConfidenceCompareToTruth::GetMapConfidenceType()
{
	return "MapConfidenceCompareToTruth";
}

SurfaceMapConfidence * MapConfidenceCompareToTruth::Clone(SurfaceMap * currMap)
{
	return new MapConfidenceCompareToTruth(currMap, m_pTrueMap, m_sDistName1, m_sDistName2, 
										   m_SampleSet1, m_SampleSet2, m_fUnadjustedThreshold);	// arbitrary set
}

double MapConfidenceCompareToTruth::CalculateErrorAtSample(const SurfaceSample & sample)
{
	SampledSurface * surf = m_pSurfaceMap->GetSurface(sample);
	//double scale = m_pSurfaceMap->GetSurface(sample)->AdjustedRadius(m_fUnadjustedThreshold);	
	double scale = sqrt(surf->Area());
	if (surf==m_pSurfaceMap->GetSurface(1))
	{
		SurfaceDistance * dist = m_pSurfaceMap->GetSurface(1)->GetOnTheFlyDistanceMetric(-1., m_sDistName2, -1);
		SurfaceSample sampInv = m_pTrueMap->InverseMap(sample);
		assert(!sampInv.Invalid());
		SurfaceSample sampForw = m_pSurfaceMap->ForwardMap(sampInv);
		assert(!sampForw.Invalid());
		assert(sampForw.CheckMeshIsSame(surf->GetMesh()));
		assert(sample.CheckMeshIsSame(surf->GetMesh()));
		return dist->Distance(sample, sampForw) / scale;
	}
	else
	{
		assert(surf==m_pSurfaceMap->GetSurface(0));
		SurfaceDistance* dist = m_pSurfaceMap->GetSurface(0)->GetOnTheFlyDistanceMetric(-1., m_sDistName1, -1);
		SurfaceSample sampInv = m_pTrueMap->ForwardMap(sample);
		assert(!sampInv.Invalid());
		SurfaceSample sampForw = m_pSurfaceMap->InverseMap(sampInv);
		assert(!sampForw.Invalid());

		assert(sampForw.CheckMeshIsSame(surf->GetMesh()));
		assert(sample.CheckMeshIsSame(surf->GetMesh()));
		return dist->Distance(sample, sampForw) / scale;
	}
}

double MapConfidenceCompareToTruth::CalculateConfidenceAtSample(const SurfaceSample & sample)
{
//	SampledSurface * surf = m_pSurfaceMap->GetSurface(sample);
//	double threshold = surf->AdjustedRadius(m_fUnadjustedThreshold);
	
	if (ErrorAtSample(sample) > m_fUnadjustedThreshold)
		return 0.;
	return 1.;
}

double MapConfidenceCompareToTruth::ErrorThreshold()
{
	//return surf->AdjustedRadius(m_fUnadjustedThreshold);
	return m_fUnadjustedThreshold;
}

