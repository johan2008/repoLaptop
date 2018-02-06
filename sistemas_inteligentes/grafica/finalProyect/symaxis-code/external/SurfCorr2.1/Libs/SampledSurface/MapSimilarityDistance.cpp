#include "MapSimilarityDistance.h"

MapSimilarityDistance::MapSimilarityDistance(SurfaceMap * M1, SurfaceMap * M2,
											 MapSimInvDistApproxType type, int simFlags, double maxDistance, 
											 double similarityThreshold, double sigma,
											 const VKString & sampleSet1, const VKString & sampleSet2,
											 const VKString & distance1, const VKString & distance2)
: SurfaceMapSimilarity(M1, M2, false, false, -1, sampleSet1, "none",
					   distance1, distance2, similarityThreshold, true, true)
{
	m_fSigma = sigma;
	m_fMaxDistance = maxDistance;
	m_ApproximationType = type;
	m_SimFlags = simFlags;	
	
	if ((m_ApproximationType==MAPSIM_BIDIR) 
		|| (m_ApproximationType==MAPSIM_FAST_DIST_SUBSETS_BIDIR))
		m_CoarseSet2Name = sampleSet2;
}

MapSimilarityDistance::~MapSimilarityDistance()
{
}

VKString MapSimilarityDistance::GetSurfaceMapSimilarityType()
{
	return "MapSimilarityDistance";
}

bool MapSimilarityDistance::IsFlagSet(SimilarityFlags flag)
{
	return (m_SimFlags & flag)!=0;
}

SurfaceMapSimilarity * MapSimilarityDistance::Clone(SurfaceMap * M1, SurfaceMap * M2)
{
	return new MapSimilarityDistance(M1, M2, m_ApproximationType, m_SimFlags, m_fMaxDistance, 
									 m_fSimilarityThreshold, m_fSigma,
									 m_CoarseSet1Name, m_CoarseSet2Name, 
									 m_Distance1Name, m_Distance2Name);
}

double MapSimilarityDistance::Similarity()
{
	if (m_fSimilarityValue==-1.)
	{
		if (IsFlagSet(SIMFLAG_NO_INCONSISTENT_GENS)
			&& m_pSurfaceMap1->GetSurfaceMapType()=="MapConformal" 
			&& m_pSurfaceMap2->GetSurfaceMapType()=="MapConformal"
			&& HaveInconsistentGenerators())
			m_fSimilarityValue = 0;	
		else
			m_fSimilarityValue = CalculateTotalSimilarity();
	}
	return m_fSimilarityValue;	
}

double MapSimilarityDistance::CalculateSimilarityAtSample(const SurfaceSample & samp)
{
//	static double sigma2 = pow(.5, 2.);
	return exp(-DissimilarityAtSample(samp) / (m_fSigma*m_fSigma));
//	return 1.;
}

double MapSimilarityDistance::FastSimilarity(double sigma2, MapConformal * M1, MapConformal * M2,
											 const SurfaceSample & samp, 
											 SurfaceDistance * dist1, double distNorm)
{
	if (M1==M2)
		return 1.;
	const SurfaceSample & sampForw = M1->ForwardMap(samp);
	if (sampForw.Invalid())
		return 0;
	const SurfaceSample & sampInv = M2->InverseMap(sampForw);
	if (sampInv.Invalid())
		return 0;
	double distance = dist1->Distance(samp, sampInv) / distNorm;
	return exp(-distance / sigma2);
}

double MapSimilarityDistance::CalculateDissimilarityAtSample(const SurfaceSample & samp)
{
	assert(!IsFlagSet(SIMFLAG_SCALE_BY_CONFIDENCE));
	
	SurfaceSample samp1, samp2;
	SurfaceDistance * dist = NULL;
	double distNorm =0;
	if (m_ApproximationType==MAPSIM_FORWARD_ONLY 
		|| m_ApproximationType==MAPSIM_BIDIR)
	{
		dist = m_pSurfaceMap1->GetOtherSurface(samp)->GetOnTheFlyDistanceMetric(m_fMaxDistance, m_Distance1Name, -1);
		distNorm = sqrt(m_pSurfaceMap1->GetOtherSurface(samp)->Area());
		if (m_pSurfaceMap1->GetSurface(0) == m_pSurfaceMap1->GetSurface(samp))
		{
			samp1 = m_pSurfaceMap1->ForwardMap(samp);
			if (!samp1.Invalid())
				samp2 = m_pSurfaceMap2->ForwardMap(samp);
		}
		else
		{
			samp1 = m_pSurfaceMap1->InverseMap(samp);
			if (!samp1.Invalid())			
				samp2 = m_pSurfaceMap2->InverseMap(samp);
		}
	}
	else if (m_ApproximationType==MAPSIM_FAST_DIST_SUBSETS 
			 || m_ApproximationType==MAPSIM_FAST_DIST_SUBSETS_BIDIR)
	{			
		dist = m_pSurfaceMap1->GetSurface(samp)->GetOnTheFlyDistanceMetric(m_fMaxDistance, m_Distance1Name, -1);
		distNorm = sqrt(m_pSurfaceMap1->GetSurface(samp)->Area());		
		if (m_pSurfaceMap1->GetSurface(0) == m_pSurfaceMap1->GetSurface(samp))
		{
			samp1 = samp;
			SurfaceSample otherSamp = m_pSurfaceMap1->ForwardMap(samp);
			if (!otherSamp.Invalid())
				samp2 = m_pSurfaceMap2->InverseMap(otherSamp);
		}
		else
		{
			samp1 = samp;
			SurfaceSample otherSamp = m_pSurfaceMap1->InverseMap(samp);
			if (!otherSamp.Invalid())
				samp2 = m_pSurfaceMap2->ForwardMap(otherSamp);
		}
	}
	
	assert(dist!=NULL);
	
	double retVal=-1.;
	if (samp1.Invalid() || samp2.Invalid())
	{
		retVal =  FLT_MAX;
	}
	else
	{
//		static bool checkDistanceOnce = true;
//		if (checkDistanceOnce)	// just a debug check
//		{
//			SurfaceDistance * existingDistance = NULL;
//			if (m_ApproximationType==MAPSIM_FORWARD_ONLY 
//				|| m_ApproximationType==MAPSIM_BIDIR)
//				existingDistance = m_pSurfaceMap1->GetOtherSurface(samp)->GetDistanceMetric(m_Distance1Name);
//			else if (m_ApproximationType==MAPSIM_FAST_DIST_SUBSETS 
//					 || m_ApproximationType==MAPSIM_FAST_DIST_SUBSETS_BIDIR)
//				existingDistance = m_pSurfaceMap1->GetSurface(samp)->GetDistanceMetric(m_Distance1Name);
//			else
//				assert(false);
//
//			if (!existingDistance->Valid(samp1, samp2))
//			{
//				checkDistanceOnce = false;
//				std::cout<<"[WARNING] MapSimilarityDistance - accessing non-precomputed distances"<<std::endl;
//			}
//		}		
		
//		std::cout<<"distVal = "<<dist->Distance(samp1, samp2)<<" / "<<distNorm<<std::endl;
		double distVal = dist->Distance(samp1, samp2)/distNorm;
		if (m_fMaxDistance>=0)
		{
			assert(false);	// verify scaling below, seems sketchy
			retVal = (distVal > m_fMaxDistance) ? 1. : (distVal/m_fMaxDistance);
		}
		else
		{
			retVal = distVal;
		}
	}		
	return retVal;
}

bool MapSimilarityDistance::HaveInconsistentGenerators()
{
	bool atVertex;
	std::map<int, int> vertexToVertex12;
	std::map<int, int> vertexToVertex21;
	
	SurfaceSampleSet * genSet1;
	SurfaceSampleSet * genSet2;
	std::vector<int> * correspondences;
	assert(m_pSurfaceMap1->GetSurfaceMapType()=="MapConformal");
	((MapConformal*)m_pSurfaceMap1)->GetGeneratingCorrespondences(&correspondences, &genSet1, &genSet2);
	for (int i=0; i<(int)correspondences->size(); i+=2)
	{
		int v1 = genSet1->GetSample((*correspondences)[i+0]).NearestVertex(&atVertex);
		assert(atVertex);
		int v2 = genSet2->GetSample((*correspondences)[i+1]).NearestVertex(&atVertex);
		assert(atVertex);
		
		vertexToVertex12[v1] = v2;
		vertexToVertex21[v2] = v1;
	}
	
	assert(m_pSurfaceMap2->GetSurfaceMapType()=="MapConformal");
	((MapConformal*)m_pSurfaceMap2)->GetGeneratingCorrespondences(&correspondences, &genSet1, &genSet2);
	for (int i=0; i<(int)correspondences->size(); i+=2)
	{
		int v1 = genSet1->GetSample((*correspondences)[i+0]).NearestVertex(&atVertex);
		int v2 = genSet2->GetSample((*correspondences)[i+1]).NearestVertex(&atVertex);
		assert(atVertex);
		
		if (vertexToVertex12.find(v1)!=vertexToVertex12.end() 
			&& vertexToVertex12[v1]!=v2)
			return true;
		
		if (vertexToVertex21.find(v2)!=vertexToVertex21.end() 
			&& vertexToVertex21[v2]!=v1)
			return true;
	}
	
	return false;	//no inconsistencies found
}

int MapSimilarityDistance::NumSharedCorrespondences(MapConformal * M1, MapConformal * M2)
{
	std::vector<int> * corrs1;
	std::vector<int> * corrs2;
	SurfaceSampleSet * set11;
	SurfaceSampleSet * set12;	
	SurfaceSampleSet * set21;
	SurfaceSampleSet * set22;	
	M1->GetGeneratingCorrespondences(&corrs1, &set11, &set12);
	M2->GetGeneratingCorrespondences(&corrs2, &set21, &set22);	
//	std::cout<<"Measuring consistency: "<<std::endl;
//	for (int i=0; i<6; i+=2)
//		std::cout<<"\tM1: ["<<(*corrs1)[i]<<" -> "<<(*corrs1)[i+1]<<"]"<<std::endl;
//	for (int i=0; i<6; i+=2)
//		std::cout<<"\tM2: ["<<(*corrs2)[i]<<" -> "<<(*corrs2)[i+1]<<"]"<<std::endl;
	
	std::map<int, int> set1map12;
	std::map<int, int> set1map21;	
	int numShared = 0;
	for (int i=0; i<(int)corrs1->size(); i+=2)
	{
		set1map12[(*corrs1)[i+0]] = (*corrs1)[i+1];
		set1map21[(*corrs1)[i+1]] = (*corrs1)[i+0];
	}
	
	for (int i=0; i<(int)corrs2->size(); i+=2)
	{
		std::map<int, int>::iterator iter12 = set1map12.find((*corrs2)[i+0]);
		std::map<int, int>::iterator iter21 = set1map21.find((*corrs2)[i+1]);
		if (iter12!=set1map12.end())
		{
			if (iter12->second==(*corrs2)[i+1])
				numShared++;
			else
				return -1;
		}
		
		if (iter21!=set1map21.end() 
			&& iter21->second!=(*corrs2)[i+0])
			return -1;		
	}
	return numShared;
}






