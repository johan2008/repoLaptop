#include "MapConfidenceDiscreteArea.h"
MapConfidenceDiscreteArea::MapConfidenceDiscreteArea(SurfaceMap * corrMap,
													 double , 
													 const VKString & sampleSet1, 
													 const VKString & sampleSet2, 
													 const VKString & , 
													 const VKString & )
: SurfaceMapConfidence(corrMap, true, false, -1, sampleSet1, sampleSet2, true, false)
{
	assert(false);
}

MapConfidenceDiscreteArea::~MapConfidenceDiscreteArea()
{
}

VKString MapConfidenceDiscreteArea::GetMapValueType()
{
	return "MapConfidenceDiscreteArea";
}

SurfaceMapConfidence * MapConfidenceDiscreteArea::Clone(SurfaceMap * )
{
	assert(false);
}

double MapConfidenceDiscreteArea::CalculateConfidenceAtSample(const SurfaceSample & sample)
{
	assert(false);
}


//MapValueDiscreteArea::MapValueDiscreteArea(double neighborhoodRadius, 
//											   const VKString & sampleSet1, const VKString & sampleSet2, 
//											   MapCoarse::FineToCoarseGeneration nhdType, ValueType nhdValue,
//											   const VKString & dist1, const VKString & dist2)
//: SurfaceMapValue(sampleSet1, sampleSet2)
//{
//	m_fNeighborhoodRadius = neighborhoodRadius;
//	m_NeighborhoodType = nhdType;
//	m_NeighborhoodValue = nhdValue;
//	assert(m_NeighborhoodType==MapCoarse::F2C_FORWARD_ONLY 
//		   || m_NeighborhoodType==MapCoarse::F2C_BACKWARD_ONLY
//		   || m_NeighborhoodType==MapCoarse::F2C_FORWARD_BACKWARD);
//	
//	m_Distance1Name = dist1;
//	m_Distance2Name = dist2;
//}
//
//MapValueDiscreteArea::~MapValueDiscreteArea()
//{
//}
//
//VKString MapValueDiscreteArea::GetMapValueType()
//{
//	return "MapValueDiscreteArea";
//}
//
//double MapValueDiscreteArea::CalculateValue(SurfaceMap * map1)
//{	
//	if (map1->GetSurfaceMapType()=="MapCoarse")
//		return CalculateValueForCoarse((MapCoarse*)map1);
//	else 
//	{
//		MapCoarse coarseMap(map1->GetSurface(0)->GetSampleSet(m_SampleSet1), 
//							map1->GetSurface(1)->GetSampleSet(m_SampleSet2),
//							map1, m_NeighborhoodType, m_Distance1Name, m_Distance2Name);
//		return CalculateValueForCoarse(&coarseMap);
//	}
//}
//
//double MapValueDiscreteArea::CalculateValueForCoarse(MapCoarse * coarseMap)
//{
//	SurfaceSampleSet * samples1 = coarseMap->GetSurface(0)->GetSampleSet(m_SampleSet1);
//	SurfaceSampleSet * samples2 = coarseMap->GetSurface(1)->GetSampleSet(m_SampleSet2);
//	
//	SurfaceDistance * dist1 = coarseMap->GetSurface(0)->GetDistanceMetric(m_Distance1Name);
//	SurfaceDistance * dist2 = coarseMap->GetSurface(1)->GetDistanceMetric(m_Distance2Name);	
//	
//	double thresholdOn1 = coarseMap->GetSurface(0)->AdjustedRadius(m_fNeighborhoodRadius);
//	double thresholdOn2 = coarseMap->GetSurface(1)->AdjustedRadius(m_fNeighborhoodRadius);
//	
//	int norm=0;
//	double error=0;
//	
//	if (m_NeighborhoodType==MapCoarse::F2C_FORWARD_ONLY 
//		|| m_NeighborhoodType==MapCoarse::F2C_FORWARD_BACKWARD)
//	{
//		SurfaceSampleSet samplesOnS2;
//		for (int i=0; i<samples1->NumSamples(); i++)
//		{
//			const SurfaceSample & s1 = samples1->GetSample(i);
//			const SurfaceSample & s2 = coarseMap->ForwardMap(s1);
//			if (!s2.Invalid())
//				samplesOnS2.AddSample(s2);
//		}
//		
//		for (int i=0; i<samples2->NumSamples(); i++)
//		{
//			norm++;
//			const SurfaceSample & sEvenOn2 = samples2->GetSample(i);
//			// find local neighbors of sEvenOn2
//			int N1 = NumNeighbors(sEvenOn2, *samples2, dist2, thresholdOn2);
//			int N2 = NumNeighbors(sEvenOn2, samplesOnS2, dist2, thresholdOn2);
//			error += GetErr(N1, N2);
//		}
//	}
//	
//	if (m_NeighborhoodType==MapCoarse::F2C_BACKWARD_ONLY 
//		|| m_NeighborhoodType==MapCoarse::F2C_FORWARD_BACKWARD)
//	{
//		SurfaceSampleSet samplesOnS1;
//		for (int i=0; i<samples2->NumSamples(); i++)
//		{
//			const SurfaceSample & s2 = samples1->GetSample(i);
//			const SurfaceSample & s1 = coarseMap->InverseMap(s2);
//			if (!s1.Invalid())
//				samplesOnS1.AddSample(s1);
//		}
//		
//		for (int i=0; i<samples1->NumSamples(); i++)
//		{
//			norm++;
//			const SurfaceSample & sEvenOn1 = samples1->GetSample(i);
//			// find local neighbors of sEvenOn1
//			int N1 = NumNeighbors(sEvenOn1, *samples1, dist1, thresholdOn1);
//			int N2 = NumNeighbors(sEvenOn1, samplesOnS1, dist1, thresholdOn1);
//			error += GetErr(N1, N2);			
//		}
//	}			
//	return ((double)norm) / error;
//}
//
//int MapValueDiscreteArea::NumNeighbors(const SurfaceSample & s, const SurfaceSampleSet & sampleSet,
//										 SurfaceDistance * dist, double threshold, 
//										 std::vector<SurfaceSample> * saveNeighbors)
//{
//	int N=0;
//	for (int i=0; i<sampleSet.NumSamples(); i++)
//	{
//		if (dist->Distance(s, sampleSet.GetSample(i)) < threshold)
//		{
//			N++;
//			if (saveNeighbors!=NULL)
//				saveNeighbors->push_back(sampleSet.GetSample(i));
//		}
//	}
//	return N;
//}
//
//double MapValueDiscreteArea::GetErr(int N1, int N2)
//{
//	double val = 0;
//	if (m_NeighborhoodValue==NHD_VAL_DIFF)
//		val = ((double)N1-(double)N2);
//	else if (m_NeighborhoodValue==NHD_VAL_MIN_OVER_MAX)
//		val = 1-vkMin((double)N1, (double)N2)/vkMax((double)N1, (double)N2);
//	else
//		assert(false);
//	return val;
//}
//

