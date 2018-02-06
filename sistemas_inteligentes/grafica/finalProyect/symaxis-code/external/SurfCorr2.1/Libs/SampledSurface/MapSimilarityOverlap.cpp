#include "MapSimilarityOverlap.h"

MapSimilarityOverlap::MapSimilarityOverlap(SurfaceMap * M1, SurfaceMap * M2,
										   const VKString & sampleSet1, const VKString & sampleSet2,
										   const VKString & distance1, const VKString & distance2)
: SurfaceMapSimilarity(M1, M2, false, false, -1, sampleSet1, sampleSet2, distance1, distance2,
					   -1, true, true)
{
	assert(false);
}

MapSimilarityOverlap::~MapSimilarityOverlap()
{
}

VKString MapSimilarityOverlap::GetSurfaceMapSimilarityType()
{
	return "MapSimilarityOverlap";
}

double MapSimilarityOverlap::CalculateSimilarityAtSample(const SurfaceSample & samp)
{
	assert(false);
}


SurfaceMapSimilarity * MapSimilarityOverlap::Clone(SurfaceMap * , SurfaceMap *)
{
	assert(false);
}


//#include "MapSimilarityOverlap.h"
//
////////////// MapSimilarityOverlap ////////////////////
//MapSimilarityOverlap::MapSimilarityOverlap(MapCoarse::FineToCoarseGeneration coarseMethod, 
//										   const VKString & sampleSet1, const VKString & sampleSet2,
//										   const VKString & distance1, const VKString & distance2)
//: SurfaceMapSimilarity(sampleSet1, sampleSet2, distance1, distance2)
//{
//	m_FineToCoarseMethod = coarseMethod;
//}
//
//MapSimilarityOverlap::~MapSimilarityOverlap()
//{
//}
//
//VKString MapSimilarityOverlap::GetMapSimilarityType()
//{
//	return "MapSimilarityOverlap";
//}
//
//double MapSimilarityOverlap::Similarity(SurfaceMap * M1, SurfaceMap * M2)
//{
//	assert(M1->GetSurface(0)==M2->GetSurface(0) && M1->GetSurface(1)==M2->GetSurface(1));
//	
//	if (M1->GetSurfaceMapType()=="MapCoarse")
//		return CalculateValueForCoarse((MapCoarse *) M1, (MapCoarse *) M2);
//	else
//	{
//		MapCoarse * coarseMap1 = M1->GetOrCreateCoarseMap(m_CoarseSet1Name, m_CoarseSet2Name,
//														  m_Distance1Name, m_Distance2Name, 
//														  (int)m_FineToCoarseMethod);
//		MapCoarse * coarseMap2 = M2->GetOrCreateCoarseMap(m_CoarseSet1Name, m_CoarseSet2Name,
//														  m_Distance1Name, m_Distance2Name, 
//														  (int)m_FineToCoarseMethod);
//		return CalculateValueForCoarse(coarseMap1, coarseMap2);
//	}
//}
//
//double MapSimilarityOverlap::CalculateValueForCoarse(MapCoarse * M1, MapCoarse * M2)
//{
//	if (M1->GetValidDomain()->NumSamples() > M2->GetValidDomain()->NumSamples())
//	{		MapCoarse * tmp = M1;		M1 = M2;		M2 = tmp;	}	// swap
//	
//	const SurfaceSampleSet * domainM1 = M1->GetValidDomain();
//	
//	double value = 0;
//	for (int i=0; i<domainM1->NumSamples(); i++)
//	{
//		const SurfaceSample & s1 = M1->ForwardMap(domainM1->GetSample(i));
//		const SurfaceSample & s2 = M2->ForwardMap(domainM1->GetSample(i));		
//		if (!s1.Invalid() && !s2.Invalid() && s1==s2)
//			value+=1;
//	}
//	if (domainM1->NumSamples()==0)
//		return 0;
//	return value / (double)domainM1->NumSamples();
//}
