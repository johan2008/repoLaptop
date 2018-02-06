#include "MapConfidenceMutuallyClosest.h"

MapConfidenceMutuallyClosest::MapConfidenceMutuallyClosest(SurfaceMap * corrMap, 
														   const VKString & set1, const VKString & set2,
														   ValueType valueType, 
														   const VKString & dist1,
														   const VKString & dist2)
: SurfaceMapConfidence(corrMap, false, false, -1, set1, "none", true, true)
{
	m_pCorrespondingMCNMap = corrMap->GetOrCreateCoarseMap(set1, set2, dist1, dist2, 
														   MapCoarse::F2C_MUTUALLY_CLOSEST_NEIGHBORS);
	m_Distance1Name = dist1;
	m_Distance2Name = dist2;
	
	m_fEpsilon = 0.001;
	m_ValueType = valueType;	
}

MapConfidenceMutuallyClosest::~MapConfidenceMutuallyClosest()
{
}

VKString MapConfidenceMutuallyClosest::GetMapConfidenceType()
{
	return "MapConfidenceMutuallyClosest";
}

SurfaceMapConfidence * MapConfidenceMutuallyClosest::Clone(SurfaceMap * currMap)
{
	return new MapConfidenceMutuallyClosest(currMap, m_SampleSet1, m_SampleSet2,
											m_ValueType, m_Distance1Name, m_Distance2Name);
}

double MapConfidenceMutuallyClosest::Confidence()
{
	if (m_TotalConfidence==-1)
	{
		switch (m_ValueType)
		{
			case CMV_FRAC_MCN:
				m_TotalConfidence = CalculateTotalConfidence();
				break;
			case CMV_INV_DIST:
				m_TotalConfidence = 1. / (m_fEpsilon + Error());
				break;				
			default:
				assert(false);
		}
	}
	return m_TotalConfidence;
	
}

double MapConfidenceMutuallyClosest::CalculateConfidenceAtSample(const SurfaceSample & sample)
{
	if (m_pCorrespondingMCNMap->GetSurface(sample)==m_pCorrespondingMCNMap->GetSurface(0)
		&& !m_pCorrespondingMCNMap->ForwardMap(sample).Invalid())
		return 1.;
	else if (m_pCorrespondingMCNMap->GetSurface(sample)==m_pCorrespondingMCNMap->GetSurface(1)
				&& !m_pCorrespondingMCNMap->InverseMap(sample).Invalid())
		return 1.;
	
	return 0.;
}

double MapConfidenceMutuallyClosest::CalculateErrorAtSample(const SurfaceSample & sample)
{
	assert(false);
//	double err;
//	if (m_pCorrespondingMCNMap->GetSurface(sample)==m_pCorrespondingMCNMap->GetSurface(0)
//		&& !m_pCorrespondingMCNMap->ForwardMap(sample, &err).Invalid())
//		return err;
//	else if (m_pCorrespondingMCNMap->GetSurface(sample)==m_pCorrespondingMCNMap->GetSurface(1)
//			 && !m_pCorrespondingMCNMap->InverseMap(sample, &err).Invalid())
//		return err;
//	
//	return -2;
}

/////// Temporary some fast MCN stuff

void MapConfidenceMutuallyClosest::ConstructCoordsHolderForFastMCN(SurfaceMidEdgeConf * mc1, 
																   SurfaceMidEdgeConf * mc2,
																   SurfaceSampleSet * set1, 
																   SurfaceSampleSet * set2,
																   LinAlgComplex ** origCoords1, 
																   LinAlgComplex ** origCoords2,
																   LinAlgComplex ** transCoords1, 
																   LinAlgComplex ** transCoords2)
{
	*origCoords1 = new LinAlgComplex[set1->NumSamples()];
	*origCoords2 = new LinAlgComplex[set2->NumSamples()];
	*transCoords1 = new LinAlgComplex[set1->NumSamples()];
	*transCoords2 = new LinAlgComplex[set2->NumSamples()];
	
	for (int i=0; i<set1->NumSamples(); i++)
		(*origCoords1)[i] = mc1->GetConfCoordOriginal(set1->GetSample(i));
	
	for (int i=0; i<set2->NumSamples(); i++)	
		(*origCoords2)[i] = mc2->GetConfCoordOriginal(set2->GetSample(i));	
}

double MapConfidenceMutuallyClosest::GetFractionOfMutuallyClosest(SurfaceMidEdgeConf * mc1, 
																  SurfaceMidEdgeConf * mc2,
																  SurfaceSampleSet * set1, 
																  SurfaceSampleSet * set2,
																  LinAlgComplex * origCoords1,
																  LinAlgComplex * origCoords2,
																  LinAlgComplex * transCoords1, 
																  LinAlgComplex * transCoords2,												   
																  MapConformal * map1)
{
	MobiusTransformation m1 = mc1->GetCurrentTransform();
	MobiusTransformation m2 = mc2->GetCurrentTransform();
	if (mc1->GetCurrentTransform()!=map1->GetMobiusForSurface(0))
		mc1->Transform(map1->GetMobiusForSurface(0));
	if (mc2->GetCurrentTransform()!=map1->GetMobiusForSurface(1))
		mc2->Transform(map1->GetMobiusForSurface(1));
		
	for (int i=0; i<set1->NumSamples(); i++)
		transCoords1[i] = m1.Transform(origCoords1[i]);
	for (int i=0; i<set2->NumSamples(); i++)
		transCoords2[i] = m2.Transform(origCoords2[i]);

	std::vector<int> nearestNeihbor1;
	std::vector<int> nearestNeihbor2;
	std::vector<double> nearestNeihborDist2;
	for (int i=0; i<set1->NumSamples(); i++)
	{
		double minima = 0; 
		int minimaID = -1;
		LinAlgComplex & v1 = transCoords1[i];
		for (int j=0; j<set2->NumSamples(); j++)
		{
			LinAlgComplex & v2 = transCoords2[j];
			double dx12 = pow(v1.r-v2.r, 2) + pow(v1.i-v2.i, 2);
			if (dx12 < minima || minimaID==-1)
			{
				minimaID = j;
				minima = dx12;
			}
			
			if (i==0)
			{
				nearestNeihborDist2.push_back(dx12);
				nearestNeihbor2.push_back(i);
			}
			else
			{
				double currMin = nearestNeihborDist2[j];
				if (dx12 < currMin)
				{
					nearestNeihborDist2[j] = dx12;
					nearestNeihbor2[j] = i;
				}
			}
		}
		assert(minimaID!=-1);
		nearestNeihbor1.push_back(minimaID);
	}
	
	double mutuallyClosest =0;
	for (int i=0; i<(int)nearestNeihbor1.size(); i++)
	{
		int s1 = i;
		int s2 = nearestNeihbor1[s1];
		if (nearestNeihbor2[s2]==s1)
			mutuallyClosest+=1.;
	}
	return mutuallyClosest / (double)nearestNeihbor1.size();
}

