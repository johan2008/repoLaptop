#include "MapMultiConformal.h"
#include "MapConfidenceMultiConf.h"
#include "SurfaceMapSimilarity.h"
#include "SurfaceMapConfidence.h"
#include "Sortable.h"
#include "DistanceGeodesic.h"
#include "DistanceOnTheFly.h"
#include "AnalysisStats.h"
#include <algorithm>

int MapMultiConformal::s_iGlobalGeneratingCorrSetUniqueIndex = 0;

MapMultiConformal::MapMultiConformal(SampledSurface * surf1, SampledSurface * surf2)
: SurfaceMap(surf1, surf2)
{
//	m_bPreparedWeights = false;	
}

MapMultiConformal::MapMultiConformal(MapConformal * confMap, const VKString & multiConfConfidence)
: SurfaceMap(confMap->GetSurface(0), confMap->GetSurface(1))
{
	std::vector<MapConformal *> confMaps;	
	confMaps.push_back(confMap);
	Initialize(confMaps, multiConfConfidence, MULTICONF_GEODESIC_CENTROID);
}

MapMultiConformal::MapMultiConformal(std::vector<MapConformal*> & maps,
									 std::vector< std::map<int, double> > & weights,
									 const VKString & multiConfConfidence,
									 const VKString & interpMethodStr)
: SurfaceMap(maps[0]->GetSurface(0), maps[0]->GetSurface(1))
{
	Initialize(maps, multiConfConfidence, StrToInterMethod(interpMethodStr));
	SurfaceMapConfidence * conf = GetMapConfidenceCalculator(multiConfConfidence);
	assert(conf->GetMapConfidenceType()=="MapConfidenceMultiConf");
	((MapConfidenceMultiConf*)conf)->FillFinalWeights(weights);
}

MapMultiConformal::MapMultiConformal(MapScoredCollection * scoredCollection, 
									 const std::vector<int> & cluster, 
									 const VKString & multiConfConfidence,
									 const VKString & interpMethodStr)	// create multi-map
: SurfaceMap(scoredCollection->GetSurface(0), scoredCollection->GetSurface(1))
{
//	m_bPreparedWeights = false;	
	std::vector<MapConformal *> confMaps;
	for (int i=0; i<(int)cluster.size(); i++)
	{
		SurfaceMap * currMap = scoredCollection->GetMapByID(cluster[i]);
		assert(currMap->GetSurfaceMapType()=="MapConformal");
		confMaps.push_back((MapConformal*)currMap);
	}
	
	Initialize(confMaps, multiConfConfidence, StrToInterMethod(interpMethodStr));
}

MapMultiConformal::MapMultiConformal(MapCoarse * goodCorrs, 
									 const VKString & setName,
									 const VKString & multiConfConfidence,
									 const VKString & interpMethodStr,
									 int MAX_CONF_MAPS)
: SurfaceMap(goodCorrs->GetSurface(0), goodCorrs->GetSurface(1))
{
	bool ifTooManyTakeFarthest = true;
	std::vector<MapConformal *> confMaps;
	SurfaceSampleSet * set1 = GetSurface(0)->GetSampleSet(setName);
	SurfaceSampleSet * set2 = GetSurface(1)->GetSampleSet(setName);		
	
	std::vector<int> validCorrs = GetValidCorrs(goodCorrs, setName);
	
	assert(GetSurface(0)->GetSurfaceType()=="SurfaceMidEdgeConf");
	assert(GetSurface(1)->GetSurfaceType()=="SurfaceMidEdgeConf");	
	SurfaceMidEdgeConf * mConf1 = (SurfaceMidEdgeConf*)GetSurface(0);
	SurfaceMidEdgeConf * mConf2 = (SurfaceMidEdgeConf*)GetSurface(1);	
	int numCorrs = (int)validCorrs.size()/2;		
	
	int allConfMaps = numCorrs * (numCorrs-1) * (numCorrs-2) / 6;
	int allPairs = numCorrs * (numCorrs-1) / 2;
	//int MAX_CONF_MAPS = 10;
	
	double distNorm1 = sqrt(GetSurface(0)->Area());
	double distNorm2 = sqrt(GetSurface(1)->Area());		
//	SurfaceSampleSet * set1 = GetSurface(0)->GetSampleSet(setName);
//	SurfaceSampleSet * set2 = GetSurface(1)->GetSampleSet(setName);		
	DistanceOnTheFly * distanceMetric1 = GetSurface(0)->GetOnTheFlyDistanceMetric(-1, "default", -1);
	DistanceOnTheFly * distanceMetric2 = GetSurface(1)->GetOnTheFlyDistanceMetric(-1, "default", -1);
	bool atVertex;
	TimeProfiler distanceCalc;
	int MAX_N = set1->NumSamples() + set2->NumSamples();
	if (set1->NumSamples() > 256 || GetSurface(0)->GetMesh()->NVertices()
		|| set2->NumSamples() > 256 || GetSurface(1)->GetMesh()->NVertices())
		std::cout<<"[WARNING] Precomputing distances in MapMultiConformal"<<std::endl;
	for (int i=0; i<set1->NumSamples(); i++)
	{
		if (i%50==0)
			distanceCalc.WriteProgress("DistancePrecomp", i, MAX_N);
		distanceMetric1->PrecomputeRow(set1->GetSample(i).NearestVertex(&atVertex));
		assert(atVertex);	// if not - change to precompute all rows
	}
	for (int i=0; i<set2->NumSamples(); i++)
	{
		if (i%50==0)
			distanceCalc.WriteProgress("DistancePrecomp", set1->NumSamples()+i, MAX_N);
		distanceMetric2->PrecomputeRow(set2->GetSample(i).NearestVertex(&atVertex));
		assert(atVertex);	// if not - change to precompute all rows
	}
	std::cout<<std::endl;
	assert(set1!=NULL && set2!=NULL);
	assert(distanceMetric1!=NULL && distanceMetric2!=NULL);
	std::cout<<"Possible maps: ["<<allConfMaps<<" "<<allPairs<<" "<<numCorrs<<"] / "<<MAX_CONF_MAPS<<std::endl;
	if (allConfMaps <= MAX_CONF_MAPS)
	{
		// find set (?)
		for (int c1=0; c1<numCorrs; c1++)
		{
			for (int c2=c1+1; c2<numCorrs; c2++)
			{
				for (int c3=c2+1; c3<numCorrs; c3++)	// for every possible 3 corrs, create conf map
				{
					std::vector<int> corrs;
					corrs.push_back(validCorrs[c1*2+0]);			corrs.push_back(validCorrs[c1*2+1]);
					corrs.push_back(validCorrs[c2*2+0]);			corrs.push_back(validCorrs[c2*2+1]);
					corrs.push_back(validCorrs[c3*2+0]);			corrs.push_back(validCorrs[c3*2+1]);
					confMaps.push_back(new MapConformal(mConf1, mConf2, corrs, setName));
				}
			}
		}
	}
	else if (allPairs <= MAX_CONF_MAPS)
	{
		for (int c1=0; c1<numCorrs; c1++)
		{
			SurfaceSample c1_s1 = set1->GetSample(validCorrs[c1*2+0]);
			SurfaceSample c1_s2 = set2->GetSample(validCorrs[c1*2+1]);			
			for (int c2=c1+1; c2<numCorrs; c2++)
			{
				SurfaceSample c2_s1 = set1->GetSample(validCorrs[c2*2+0]);
				SurfaceSample c2_s2 = set2->GetSample(validCorrs[c2*2+1]);
				std::vector<Sortable> sortedDists;
				for (int c3=c2+1; c3<numCorrs; c3++)
				{
					SurfaceSample c3_s1 = set1->GetSample(validCorrs[c3*2+0]);
					SurfaceSample c3_s2 = set2->GetSample(validCorrs[c3*2+1]);					
					double d11 = distanceMetric1->Distance(c1_s1, c3_s1) / distNorm1;
					double d12 = distanceMetric2->Distance(c1_s2, c3_s2) / distNorm2;
					double d21 = distanceMetric1->Distance(c2_s1, c3_s1) / distNorm1;
					double d22 = distanceMetric2->Distance(c2_s2, c3_s2) / distNorm2;
					sortedDists.push_back(Sortable(d11 + d12 + d21 + d22, NULL, c3));
				}
				if (sortedDists.size()==0)
					continue;
				std::sort(sortedDists.begin(), sortedDists.end());
				int takeMap = ifTooManyTakeFarthest ? ((int)sortedDists.size()-1) : 0; 
				int c3 = sortedDists[takeMap].id;
				std::vector<int> corrs;
				corrs.push_back(validCorrs[c1*2+0]);			corrs.push_back(validCorrs[c1*2+1]);
				corrs.push_back(validCorrs[c2*2+0]);			corrs.push_back(validCorrs[c2*2+1]);
				corrs.push_back(validCorrs[c3*2+0]);			corrs.push_back(validCorrs[c3*2+1]);
				confMaps.push_back(new MapConformal(mConf1, mConf2, corrs, setName));
			}
		}
	}
	else if (!ifTooManyTakeFarthest)
	{
		if (numCorrs > MAX_CONF_MAPS)
			std::cout<<"[WARNING] Too many conformal maps"<<std::endl;
				
		for (int c1=0; c1<numCorrs; c1++)
		{
			SurfaceSample c1_s1 = set1->GetSample(validCorrs[c1*2+0]);
			SurfaceSample c1_s2 = set2->GetSample(validCorrs[c1*2+1]);
			std::vector<Sortable> sortedDists;
			for (int c2=0; c2<numCorrs; c2++)
			{
				if (c2==c1)
					continue;
				SurfaceSample c2_s1 = set1->GetSample(validCorrs[c2*2+0]);
				SurfaceSample c2_s2 = set2->GetSample(validCorrs[c2*2+1]);			
				double d1 = distanceMetric1->Distance(c1_s1, c2_s1) / distNorm1;
				double d2 = distanceMetric2->Distance(c1_s2, c2_s2) / distNorm2;
				sortedDists.push_back(Sortable(d1+d2, NULL, c2));
			}
			if (sortedDists.size()<2)
				continue;
			std::sort(sortedDists.begin(), sortedDists.end());
			
			int c2 = sortedDists[0].id;
			int c3 = sortedDists[1].id;
			std::vector<int> corrs;
			corrs.push_back(validCorrs[c1*2+0]);			corrs.push_back(validCorrs[c1*2+1]);
			corrs.push_back(validCorrs[c2*2+0]);			corrs.push_back(validCorrs[c2*2+1]);
			corrs.push_back(validCorrs[c3*2+0]);			corrs.push_back(validCorrs[c3*2+1]);
			confMaps.push_back(new MapConformal(mConf1, mConf2, corrs, setName));			
		}
	}
	else if (ifTooManyTakeFarthest)
	{
		if (numCorrs > MAX_CONF_MAPS)
			std::cout<<"[WARNING] Too many conformal maps"<<std::endl;
		
		for (int c1=0; c1<numCorrs; c1++)
		{
			SurfaceSample c1_s1 = set1->GetSample(validCorrs[c1*2+0]);
			SurfaceSample c1_s2 = set2->GetSample(validCorrs[c1*2+1]);
			std::vector<Sortable> sortedDists;
			for (int c2=0; c2<numCorrs; c2++)
			{
				if (c2==c1)
					continue;

				SurfaceSample c2_s1 = set1->GetSample(validCorrs[c2*2+0]);
				SurfaceSample c2_s2 = set2->GetSample(validCorrs[c2*2+1]);			
				double d1 = distanceMetric1->Distance(c1_s1, c2_s1) / distNorm1 + .0000001;
				double d2 = distanceMetric2->Distance(c1_s2, c2_s2) / distNorm2 + .0000001;
				
				for (int c3=0; c3<numCorrs; c3++)
				{
					if (c3==c2 || c3==c1)
						continue;
					
					SurfaceSample c3_s1 = set1->GetSample(validCorrs[c3*2+0]);
					SurfaceSample c3_s2 = set2->GetSample(validCorrs[c3*2+1]);			

					double d3 = distanceMetric1->Distance(c1_s1, c3_s1) / distNorm1 + .0000001;
					double d4 = distanceMetric2->Distance(c1_s2, c3_s2) / distNorm2 + .0000001;
					double d5 = distanceMetric1->Distance(c2_s1, c3_s1) / distNorm1 + .0000001;
					double d6 = distanceMetric2->Distance(c2_s2, c3_s2) / distNorm2 + .0000001;
					
					double val = 1./(d1*d1)+1./(d2*d2)+1./(d3*d3)+1./(d4*d4) + 1./(d5*d5) + 1./(d6*d6);
					//double val = 1./(d1)+1./(d2)+1./(d3)+1./(d4) + 1./(d5) + 1./(d6);
					
					sortedDists.push_back(Sortable(val, NULL, c2 * numCorrs + c3));
				}	
			}
			if (sortedDists.size()<2)
				continue;
			std::sort(sortedDists.begin(), sortedDists.end());
			
			int c3 = sortedDists[0].id % numCorrs;
			int c2 = sortedDists[0].id / numCorrs;
//			SurfaceSample c2_s1 = set1->GetSample(validCorrs[c2*2+0]);
//			SurfaceSample c2_s2 = set2->GetSample(validCorrs[c2*2+1]);			
//			SurfaceSample c3_s1 = set1->GetSample(validCorrs[c3*2+0]);
//			SurfaceSample c3_s2 = set2->GetSample(validCorrs[c3*2+1]);			
//			
//			std::cout<<"Distance["<<c1<<", "<<c2<<"] = "<<distanceMetric1->Distance(c1_s1, c2_s1);
//			std::cout<<" and "<<distanceMetric2->Distance(c1_s2, c2_s2)<<std::endl;
//			std::cout<<"Distance["<<c1<<", "<<c3<<"] = "<<distanceMetric1->Distance(c1_s1, c3_s1);
//			std::cout<<" and "<<distanceMetric2->Distance(c1_s2, c3_s2)<<std::endl;
//			std::cout<<"Distance["<<c3<<", "<<c2<<"] = "<<distanceMetric1->Distance(c3_s1, c2_s1);
//			std::cout<<" and "<<distanceMetric2->Distance(c3_s2, c2_s2)<<std::endl;

			std::vector<int> corrs;
			corrs.push_back(validCorrs[c1*2+0]);			corrs.push_back(validCorrs[c1*2+1]);
			corrs.push_back(validCorrs[c2*2+0]);			corrs.push_back(validCorrs[c2*2+1]);
			corrs.push_back(validCorrs[c3*2+0]);			corrs.push_back(validCorrs[c3*2+1]);
			confMaps.push_back(new MapConformal(mConf1, mConf2, corrs, setName));			
		}
	}
	else
		assert(false);
	std::cout<<"Initialized with "<<confMaps.size()<<" conformal maps"<<std::endl;
	Initialize(confMaps, multiConfConfidence, StrToInterMethod(interpMethodStr));
	MapConfidenceMultiConf * conf = ((MapConfidenceMultiConf *)GetMapConfidenceCalculator(multiConfConfidence));
	conf->ClearMultiConfWeights();
}

MapMultiConformal::MapMultiConformal(MapCoarse * goodCorrs, const VKString & setName,
									 const VKString & multiConfConfidence,
									 const VKString & interpMethodStr, 
									 std::vector<int> & triplets, std::vector<double> & weights)
: SurfaceMap(goodCorrs->GetSurface(0), goodCorrs->GetSurface(1))
{
	std::vector<int> validCorrs = GetValidCorrs(goodCorrs, setName);
	std::vector<MapConformal *> confMaps;
	m_MapWeight = weights;
	
	SurfaceMidEdgeConf * mConf1 = (SurfaceMidEdgeConf*)GetSurface(0);
	SurfaceMidEdgeConf * mConf2 = (SurfaceMidEdgeConf*)GetSurface(1);
	for (int i=0; i<(int)triplets.size(); i+=3)
	{
		std::vector<int> corrs;
		corrs.push_back(validCorrs[triplets[i+0]*2+0]);			corrs.push_back(validCorrs[triplets[i+0]*2+1]);
		corrs.push_back(validCorrs[triplets[i+1]*2+0]);			corrs.push_back(validCorrs[triplets[i+1]*2+1]);
		corrs.push_back(validCorrs[triplets[i+2]*2+0]);			corrs.push_back(validCorrs[triplets[i+2]*2+1]);
		
		confMaps.push_back(new MapConformal(mConf1, mConf2, corrs, setName));
	}
	
	Initialize(confMaps, multiConfConfidence, StrToInterMethod(interpMethodStr));
	MapConfidenceMultiConf * conf = ((MapConfidenceMultiConf *)GetMapConfidenceCalculator(multiConfConfidence));
	conf->ClearMultiConfWeights();
}

void MapMultiConformal::SetGeodesicCentroidOutlierThreshold(double threshold)
{
	m_fGeodesicCentroidOutlierThreshold = threshold;
}

std::vector<int> MapMultiConformal::GetValidCorrs(MapCoarse * goodCorrs, const VKString & setName)
{
	SurfaceSampleSet * set1 = GetSurface(0)->GetSampleSet(setName);
	SurfaceSampleSet * set2 = GetSurface(1)->GetSampleSet(setName);		
	const SurfaceSampleSet * domainCorrs = goodCorrs->GetValidDomain();
	
	std::vector<int> validCorrs;
	for (int i=0; i<domainCorrs->NumSamples(); i++)	// find valid samples
	{
		SurfaceSample original = domainCorrs->GetSample(i);
		SurfaceSample mapped = goodCorrs->ForwardMap(original);
		if (!mapped.Invalid())
		{
			int origID = -1;
			for (int j=0; j<set1->NumSamples(); j++)
			{
				if (original==set1->GetSample(j))
					origID = j;
			}
			
			int mapID = -1;
			for (int j=0; j<set2->NumSamples(); j++)
			{
				if (mapped==set2->GetSample(j))
					mapID = j;
			}
			assert(mapID!=-1 && origID!=-1);
			validCorrs.push_back(origID);	
			validCorrs.push_back(mapID);
		}
	}
	return validCorrs;
}

MapMultiConformal::~MapMultiConformal()
{
}

void MapMultiConformal::SetInterpolationDetails(InterpolationMethod interpolationMethod)
{
	m_InterpolationMethod = interpolationMethod;
}

void MapMultiConformal::SetInterpolationDetails(const VKString & interpMethodStr)
{
	SetInterpolationDetails(StrToInterMethod(interpMethodStr));
}

MapMultiConformal::InterpolationMethod MapMultiConformal::StrToInterMethod(const VKString & interpMethodStr)
{
	if (interpMethodStr=="GeodesicCentroid")
		return MapMultiConformal::MULTICONF_GEODESIC_CENTROID;
	else if (interpMethodStr=="AverageMobius")
		return MapMultiConformal::MULTICONF_AVERAGE_MOBIUS;
	else if (interpMethodStr=="TakeBest")
		return MapMultiConformal::MULTICONF_TAKE_BEST;
	else
		assert(false);
}

void MapMultiConformal::Initialize(std::vector<MapConformal *> & confMaps, 
								   const VKString & multiConfConfidence,
								   InterpolationMethod interpMethodStr)
{
	m_fGeodesicCentroidOutlierThreshold = 1.;
	
	m_bModifiedGenCorrSet = false;
	EnableCache();
	m_pCoarseForGenerators = NULL;
	SetInterpolationDetails(interpMethodStr);
	m_MultiConformalConfidenceName = multiConfConfidence;
	
	m_ConformalMaps = confMaps;	
}

SurfaceSample MapMultiConformal::ForwardMap(const SurfaceSample & s)
{
	if (IsCacheEnabled())
	{
		bool valid;
		SurfaceSample retVal = GetCachedForw(s, valid);
		if (valid)
			return retVal;
	}
	SurfaceSample retVal = MultiConfMap(s, true);

	if (IsCacheEnabled() && !IsCacheLockedForSave() && !retVal.Invalid())
		AddSampleToForwMap(s, retVal);
	
	return retVal;
}

SurfaceSample MapMultiConformal::InverseMap(const SurfaceSample & s)
{
	assert(!s.Invalid());	
	if (IsCacheEnabled())
	{
		bool valid;
		SurfaceSample retVal = GetCachedBack(s, valid);
		if (valid)
			return retVal;
	}
	
	SurfaceSample retVal = MultiConfMap(s, false);
	
//	assert(!retVal.Invalid());
	
	if (IsCacheEnabled() && !IsCacheLockedForSave() && !retVal.Invalid())
		AddSampleToBackMap(s, retVal);
	
	return retVal;	
}

void MapMultiConformal::FindOptimalWeights()
{
	SurfaceMapConfidence * confidence = GetMapConfidenceCalculator(m_MultiConformalConfidenceName);
	assert(confidence->GetMapConfidenceType()=="MapConfidenceMultiConf");
	MapConfidenceMultiConf * multiConfConfidence = (MapConfidenceMultiConf*)confidence;
	if (!multiConfConfidence->WeightsReady())
		multiConfConfidence->PrepareFinalWeights();
	
}

void MapMultiConformal::SmoothOptimalWeights()
{
	SurfaceMapConfidence * confidence = GetMapConfidenceCalculator(m_MultiConformalConfidenceName);
	assert(confidence->GetMapConfidenceType()=="MapConfidenceMultiConf");
	MapConfidenceMultiConf * multiConfConfidence = (MapConfidenceMultiConf*)confidence;
	assert(multiConfConfidence->WeightsReady());
	multiConfConfidence->SmoothFinalWeights();	
}

SurfaceSample MapMultiConformal::MultiConfMap(const SurfaceSample & s, bool forward)
{
	if (m_InterpolationMethod==MULTICONF_GEODESIC_CENTROID
		|| m_InterpolationMethod==MULTICONF_TAKE_BEST)
		return MultiConfMapGeodesicCentroid(s, forward);
	else if (m_InterpolationMethod==MULTICONF_AVERAGE_MOBIUS)
		return MultiConfMapAverageMobius(s, forward);
	else
		assert(false);
}

SurfaceSample MapMultiConformal::MultiConfMapGeodesicCentroid(const SurfaceSample & s, bool forward)
{
	std::vector<double> w;
	std::vector<SurfaceSample> samples;
	int bestID=0;
	for (int i=0; i<(int)m_ConformalMaps.size(); i++)
	{
		SurfaceSample currSamp;
		if (forward)
			currSamp =  m_ConformalMaps[i]->ForwardMap(s);
		else
			currSamp = m_ConformalMaps[i]->InverseMap(s);
		if (!currSamp.Invalid())
		{
			samples.push_back(currSamp);
			w.push_back(Weight(i, s));
			if (w[w.size()-1] > w[bestID])
				bestID = w.size()-1;
		}
	}
	
	if (samples.size()==0)
	{
		std::cout<<"[WARNING] MapMultiConformal does not have an appropriate conformal map.";
		std::cout<<" Vertex="<<s.NearestVertex()<<std::endl;
		return SurfaceSample();
	}
	
	if (m_InterpolationMethod==MULTICONF_TAKE_BEST)
	{
		assert(bestID>=0 && bestID < (int)samples.size());
		return samples[bestID];
	}
	
	SampledSurface::GeodesicCentroidAlgorithmType g=SampledSurface::GEO_CENTROID_EUCLIDEAN_APPROX;
	return GetOtherSurface(s)->WeightedGeodesicCentroid(samples, w, g, m_fGeodesicCentroidOutlierThreshold);
}

SurfaceSample MapMultiConformal::MultiConfMapAverageMobius(const SurfaceSample & s, bool forward)
{
	assert(m_ConformalMaps.size()>0);
	// Calculate normalized weights
	double norm = 0;
	std::vector<double> w;
	int bestID = 0;	
	for (int i=0; i<(int)m_ConformalMaps.size(); i++)	
	{
		w.push_back(Weight(i, s));
		if (w[i] > w[bestID])
			bestID = i;		
	}

	double thresh = .99 / m_ConformalMaps.size();
//	double thresh = .8 * w[bestID];
	int non_zero=0;
	for (int i=0; i<(int)m_ConformalMaps.size(); i++)
	{
//		if (w[i] < w[bestID] * nullBelow)
		if (w[i] < thresh)
			w[i] = 0;
		else
			non_zero++;
		norm += w[i];
	}
	
	assert (norm!=0);

	for (int i=0; i<(int)m_ConformalMaps.size(); i++)	
		w[i] /= norm;
		
	// Find Mobius transformations 
	//		(identity for 2nd surface, and just interpolate combined)
	std::vector<MobiusTransformation> xforms;
	for (int i=0; i<(int)m_ConformalMaps.size(); i++)
	{
		MapConformal * confMap = m_ConformalMaps[i];
		MobiusTransformation m_i = confMap->m_MapTransform2.Inverse() * confMap->m_MapTransform1;
		xforms.push_back(m_i);
	}
	
	// Interpolate Mobius Transformations
//	MobiusTransformation interpMobius;// = xforms[bestID];
//	if (m_InterpolationMethod==MULTICONF_AVERAGE_MOBIUS_EXP)
//	//	interpMobius = MobiusTransformationInterpolator::InterpExp(xforms, w);
//		interpMobius = MobiusTransformationInterpolator::InterpExpShiftI(xforms, w);
//	else if  (m_InterpolationMethod==MULTICONF_AVERAGE_MOBIUS_NAIVE)
//		interpMobius = MobiusTransformationInterpolator::InterpNaive(xforms, w);
//	else
//		assert(false);
	MobiusTransformation interpMobius;
	interpMobius = MobiusTransformationInterpolator::InterpExpShiftI(xforms, w);
	
	// map point by interpolated
	static MobiusTransformation identity;
	return MapConformal::MapSample(s, interpMobius, identity, 
								   m_ConformalMaps[0]->m_pConfM1, m_ConformalMaps[0]->m_pConfM2);
}

VKString MapMultiConformal::GetSurfaceMapType()
{
	return "MapMultiConformal";
}

MapCoarse * MapMultiConformal::GetCoarseMapForSamples()
{
	if (m_pCoarseForGenerators==NULL)
	{
		SurfaceSampleSet * sampSet1 = new SurfaceSampleSet();
		SurfaceSampleSet * sampSet2 = new SurfaceSampleSet();
		std::vector<int> corrs;
		for (int i=0; i<(int)m_ConformalMaps.size(); i++)
		{
			std::vector<int> * genCorrs;
			SurfaceSampleSet * genSamples1;
			SurfaceSampleSet * genSamples2;			
			m_ConformalMaps[i]->GetGeneratingCorrespondences(&genCorrs, &genSamples1, &genSamples2);
			assert(genSamples1!=NULL && genSamples2!=NULL && genCorrs!=NULL);
			for (int corrID=0; corrID<(int)(*genCorrs).size(); corrID+=2)
			{
				int fromSampID = (*genCorrs)[corrID+0];
				int toSampID = (*genCorrs)[corrID+1];
				assert(toSampID!=-1);
				const SurfaceSample & fromSamp = genSamples1->GetSample(fromSampID);
				const SurfaceSample & toSamp = genSamples2->GetSample(toSampID);
//				std::cout<<"\tMap: "<<fromSampID<<" -> "<<toSampID<<std::endl;

				if (!sampSet1->ContainsExact(fromSamp) && !sampSet2->ContainsExact(toSamp))
				{
					int newlyAddedFrom = sampSet1->AddSample(fromSamp);
					int newlyAddedTo = sampSet2->AddSample(toSamp);
					assert(newlyAddedFrom==(int)corrs.size());
					corrs.push_back(newlyAddedTo);
				}
				else if (!sampSet1->ContainsExact(fromSamp))
				{
					sampSet1->AddSample(fromSamp);
					corrs.push_back(-1);
				}
				else if (!sampSet2->ContainsExact(toSamp))
				{
					sampSet2->AddSample(toSamp);
				}			
			}
			m_pCoarseForGenerators = new MapCoarse(m_pM1, m_pM2, sampSet1, sampSet2);
			m_pCoarseForGenerators->SetFinalCorrMap(corrs);			
		}
	}
	return m_pCoarseForGenerators;
}

void MapMultiConformal::SaveMap(std::ofstream & textStream)
{
	textStream<<"Format MapMultiConformal\n";
	textStream<<"InterpMethod "<<m_InterpolationMethod<<"\n";

	textStream<<"MultiConfConfidenceName "<<m_MultiConformalConfidenceName.c_str()<<" ";
	
	textStream<<"ModifiedGenSet "<<m_bModifiedGenCorrSet<<"\n";
	if (m_bModifiedGenCorrSet)
	{
		textStream<<"ModifiedName "<<m_sModifiedGenCorrSetName.c_str()<<" ";
		textStream<<"MaxID "<<s_iGlobalGeneratingCorrSetUniqueIndex<<"\n";
		GetSurface(0)->GetSampleSet(m_sModifiedGenCorrSetName)->WriteSet(textStream);
		GetSurface(1)->GetSampleSet(m_sModifiedGenCorrSetName)->WriteSet(textStream);
	}
	
	textStream<<"ConfMaps "<<m_ConformalMaps.size()<<"\n";
	for (int i=0; i<(int)m_ConformalMaps.size(); i++)
		((SurfaceMap*)m_ConformalMaps[i])->SaveMap(textStream);
	
	// weights and confidence values
	textStream<<"SmoothedWeights "<<0<<"\n";
	SurfaceMapConfidence * confidence = GetMapConfidenceCalculator(m_MultiConformalConfidenceName);
	assert(confidence->GetMapConfidenceType()=="MapConfidenceMultiConf");
	MapConfidenceMultiConf * multiConfConfidence = (MapConfidenceMultiConf*)confidence;
	multiConfConfidence->SaveConfidence(textStream);	
	
	textStream<<"GeoCtrOutlierThreshold "<<m_fGeodesicCentroidOutlierThreshold<<"\n";

	textStream<<"CacheExists "<<IsCacheEnabled()<<"\n";
	if (IsCacheEnabled())
		WriteCache(textStream);
}

void MapMultiConformal::LoadMap(std::ifstream & textStream)
{
	int interpMethod;
	VKStringList confidences;
	VKString similarity;
	std::vector<MapConformal *> confMaps;
	
	std::string emptyStr;
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="Format");
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="MapMultiConformal");	
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="InterpMethod");	
	
	textStream>>interpMethod;
	
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="MultiConfConfidenceName");
	textStream>>emptyStr;	m_MultiConformalConfidenceName = VKString(emptyStr.c_str());
	
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="ModifiedGenSet");	
	textStream>>m_bModifiedGenCorrSet;
	if (m_bModifiedGenCorrSet)
	{
		textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="ModifiedName");	
		textStream>>emptyStr;	m_sModifiedGenCorrSetName = VKString(emptyStr.c_str());
		textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="MaxID");	
		int globalID = 0;		textStream>>globalID;
		s_iGlobalGeneratingCorrSetUniqueIndex = vkMax(s_iGlobalGeneratingCorrSetUniqueIndex, globalID);
		
		SurfaceSampleSet * sampleSet1 = new SurfaceSampleSet();
		sampleSet1->LoadSet(textStream, GetSurface(0)->GetMesh());
		SurfaceSampleSet * sampleSet2 = new SurfaceSampleSet();
		sampleSet2->LoadSet(textStream, GetSurface(1)->GetMesh());
		GetSurface(0)->AddSampleSet(m_sModifiedGenCorrSetName, sampleSet1);
		GetSurface(1)->AddSampleSet(m_sModifiedGenCorrSetName, sampleSet2);
	}
	
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="ConfMaps");
	int numMaps;	textStream>>numMaps;
	for (int i=0; i<numMaps; i++)
	{
		SurfaceMap * newMap = SurfaceMap::CreateMap(textStream, GetSurface(0), GetSurface(1));
		assert(newMap->GetSurfaceMapType()=="MapConformal");
		confMaps.push_back((MapConformal*)newMap);
	}
	
	// weights and confidence values
	textStream>>emptyStr;	assert(strcmp(emptyStr.c_str(), "SmoothedWeights")==0);
	bool nothing;
	textStream>>nothing;
	SurfaceMapConfidence * confidence = GetMapConfidenceCalculator(m_MultiConformalConfidenceName);
	assert(confidence->GetMapConfidenceType()=="MapConfidenceMultiConf");
	MapConfidenceMultiConf * multiConfConfidence = (MapConfidenceMultiConf*)confidence;
	multiConfConfidence->LoadConfidence(textStream);

	Initialize(confMaps, m_MultiConformalConfidenceName, (InterpolationMethod)interpMethod);
	
	textStream>>emptyStr;	assert(strcmp(emptyStr.c_str(), "GeoCtrOutlierThreshold")==0);	
	textStream>>m_fGeodesicCentroidOutlierThreshold;
	
	textStream>>emptyStr;	assert(strcmp(emptyStr.c_str(), "CacheExists")==0);	
	bool cacheSaved;	textStream>>cacheSaved;

	if (cacheSaved)
		ReadCache(textStream);
}

void MapMultiConformal::DrawMultiCorrespondenceForSelected(AnalysisWindow * window, ParamParser * params,
														   const VKString & renderingParams,
														   const VKString & surfaceName)
{
	DrawMultiCorrespondenceForSelected(window, params, renderingParams, 
									   surfaceName, GetSurface(0), GetSurface(1));
	DrawMultiCorrespondenceForSelected(window, params, renderingParams, 
									   surfaceName, GetSurface(1), GetSurface(0));
}

void MapMultiConformal::DrawMultiCorrespondenceForSelected(AnalysisWindow * window, ParamParser * params,
														   const VKString & renderingParams, 
														   const VKString & surfaceName,
														   SampledSurface * surfaceFrom, SampledSurface * surfaceTo)
{
	bool renderClosestOnly = true;
	bool renderBestTriplet = false;
	bool valid;
	VKString singleCorrColor = params->GetStrValue("RendererDefault", "SingleCorrColor", valid);

	int selectedVertex = surfaceFrom->GetSelectedVertex(params);
	if (selectedVertex!=-1)
	{
//		std::cout<<"singleCorrColor="<<singleCorrColor.c_str()<<std::endl;		
		std::vector<double> weights;
		SurfaceSample currSamp(selectedVertex, surfaceFrom->GetMesh());
		if (singleCorrColor=="none")
		{
			DrawCorrespondenceForSelected(params, renderingParams, surfaceName, surfaceFrom, surfaceTo);
		}
		else 
		{
			// if rendering only some conformal maps depending on distance to generators
			std::vector<double> minGeoDistances;
			int minGeoDistID = 0;
			std::vector<int> * corrs;
			SurfaceSampleSet * sampSet1;
			SurfaceSampleSet * sampSet2;
			int surfOffset = (surfaceFrom==m_ConformalMaps[0]->GetSurface(0)) ? 0 : 1;
			DistanceOnTheFly * distFrom = surfaceFrom->GetOnTheFlyDistanceMetric(-1, "default", -1);
			assert(selectedVertex!=-1);
			SurfaceSample sampFrom(selectedVertex, surfaceFrom->GetMesh());
			for (int i=0; i<(int)m_ConformalMaps.size(); i++)
			{
				// estimate distance to the nearest generator correspondence
				double currDist = -1;
				m_ConformalMaps[i]->GetGeneratingCorrespondences(&corrs, &sampSet1, &sampSet2);
				for (int j=0; j<3; j+=2)
				{
					SurfaceSample sampTo;
					if (surfOffset==0)
						sampTo = sampSet1->GetSample((*corrs)[j]);
					else
						sampTo = sampSet2->GetSample((*corrs)[j+1]);
					double myDist = distFrom->Distance(sampFrom, sampTo);
					if (j==0 || myDist < currDist)
						currDist = myDist;
				}
				//std::cout<<"DistToMap["<<i<<"] = "<<currDist<<std::endl;
				minGeoDistances.push_back(currDist);
				if (minGeoDistances[minGeoDistID] > minGeoDistances[minGeoDistances.size()-1])
					minGeoDistID = (int)(minGeoDistances.size()-1);
			}
			//std::cout<<"Best triplet ID = "<<minGeoDistID<<std::endl;
			//std::cout<<"Best Triplet Distance = "<<minGeoDistances[minGeoDistID]<<std::endl;
								
			
			if (singleCorrColor!="CombinedWeight")
			{
				for (int i=0; i<(int)m_ConformalMaps.size(); i++)
				{
					if (singleCorrColor=="Consistency")
						weights.push_back(GetSampleConsistency(i, currSamp));
					else if (singleCorrColor=="DistanceToGenerators")
						weights.push_back(GetSampleConfidence(i, 1, currSamp));
					else if (singleCorrColor=="MapConfidence")
						weights.push_back(GetSampleConfidence(i, 0, currSamp));
				}
			}
			else
			{				
				R3Point euclCentroid(0,0,0);				
				weights.clear();			
				double sumToOneNorm=0;
				
				for (int i=0; i<(int)m_ConformalMaps.size(); i++)
				{
					sumToOneNorm += Weight(i, currSamp);
					weights.push_back(Weight(i, currSamp));
				}
				
				int zerodOut=0;
				for (int i=0; i<(int)weights.size(); i++)
					if ((weights[i]/sumToOneNorm) < (1./(double)weights.size()))
					{
						weights[i] = 0;
						zerodOut++;
					}
				MapConfidenceMultiConf::NormalizeWeights(weights);
				for (int i=0; i<(int)weights.size(); i++)
				{
					if (weights[i]>0)
					{
						R3Point p = m_ConformalMaps[i]->ForwardMap(currSamp).GetPosition();
						euclCentroid = euclCentroid + weights[i] * p;
					}
				}
				double norm=0;
				std::vector<Sortable> sortedWeights;
				for (int i=0; i<(int)weights.size(); i++)
				{
					sortedWeights.push_back(Sortable(weights[i], NULL, i));
					norm = vkMax(weights[i], norm);
				}
				std::sort(sortedWeights.begin(), sortedWeights.end());
				int lastID = sortedWeights.size()-1;
				if (renderBestTriplet)
					m_ConformalMaps[sortedWeights[lastID].id]->DrawGenerators(window, params, renderingParams, surfaceName, 
																			  3., 0, 1., 0, true, true);
				SurfaceMapConfidence * confidence =  m_ConformalMaps[lastID]->GetMapConfidenceCalculator("MapConfidence_TriangleArea");
				assert(confidence!=NULL);
				std::cout<<"ConfidenceMax="<<confidence->ConfidenceAtSampleNoSmoothing(currSamp)<<std::endl;
				
				assert(norm>0);
				for (int i=0; i<(int)weights.size(); i++)
					weights[i] /= norm;
				
			}
			assert(weights.size()==m_ConformalMaps.size());
			double epsilon = 0;
			int numConfMaps=0;
			for (int i=0; i<(int)m_ConformalMaps.size(); i++)
			{
				glColor3d(1-weights[i], weights[i], 0);
				if (!renderClosestOnly || vkAbs(minGeoDistances[i] - minGeoDistances[minGeoDistID]) <= epsilon)
				{
					numConfMaps++;
					m_ConformalMaps[i]->DrawCorrespondenceForSelected(params, "RendererDefault", "none");
				}
			}
			std::cout<<"Rendered "<<numConfMaps<<" maps"<<std::endl;
		}
	}	
}

double MapMultiConformal::GetSampleConsistency(int confMapID, const SurfaceSample & samp)
{
	SurfaceMapConfidence * confidence = GetMapConfidenceCalculator(m_MultiConformalConfidenceName);
	assert(confidence->GetMapConfidenceType()=="MapConfidenceMultiConf");
	MapConfidenceMultiConf * multiConfConfidence = (MapConfidenceMultiConf*)confidence;
	
	return multiConfConfidence->GetConsistency(confMapID, samp) /  m_ConformalMaps.size();
}

double MapMultiConformal::GetSampleConfidence(int conformalID, int confidenceID, const SurfaceSample & samp)
{	
	SurfaceMapConfidence * confidence = GetMapConfidenceCalculator(m_MultiConformalConfidenceName);
	assert(confidence->GetMapConfidenceType()=="MapConfidenceMultiConf");
	MapConfidenceMultiConf * multiConfConfidence = (MapConfidenceMultiConf*)confidence;
	
	return multiConfConfidence->GetConfidence(conformalID, confidenceID, samp);
}

double MapMultiConformal::Weight(int mapID, const SurfaceSample & samp)
{
	FindOptimalWeights();
	SmoothOptimalWeights();
	// weights and confidence values
	SurfaceMapConfidence * confidence = GetMapConfidenceCalculator(m_MultiConformalConfidenceName);
	assert(confidence->GetMapConfidenceType()=="MapConfidenceMultiConf");
	MapConfidenceMultiConf * multiConfConfidence = (MapConfidenceMultiConf*)confidence;
	return multiConfConfidence->Weight(mapID, samp);
}

void MapMultiConformal::FillPerVertexWeights(int confMapID, const VKString & weightValues, SampledSurface * surface,
											 std::vector<double> & perVertexWeights)
{	
	for (int i=0; i<surface->GetMesh()->NVertices(); i++)
	{
		SurfaceSample atVertex = SurfaceSample(i, surface->GetMesh());
		if (weightValues=="MultiMapDistGen")
			perVertexWeights.push_back(GetSampleConfidence(confMapID%m_ConformalMaps.size(), 1, atVertex));
		else if (weightValues=="MultiMapConsistency")
			perVertexWeights.push_back(GetSampleConsistency(confMapID%m_ConformalMaps.size(), atVertex));
		else if (weightValues=="MultiMapConfidence")
			perVertexWeights.push_back(GetSampleConfidence(confMapID%m_ConformalMaps.size(), 0, atVertex));
		else if (weightValues=="MultiMapFinal")
			perVertexWeights.push_back(Weight(confMapID, atVertex));
		else
			assert(false);
	}
}

void MapMultiConformal::LocallyRefineCorrs(const VKString & confidenceName,
										   const VKString & searchInSet)
{
	std::vector<int> * corrs;
	SurfaceSampleSet * sampleSet1;
	SurfaceSampleSet * sampleSet2;
	
	if (m_pCoarseForGenerators!=NULL)
	{
		delete m_pCoarseForGenerators;
		m_pCoarseForGenerators = NULL;
	}
	ClearCache();
//	m_bPreparedWeights = false;
	m_bModifiedGenCorrSet = true;
	m_sModifiedGenCorrSetName = VKString("GenCorrSetRefinedForCluster_");
	m_sModifiedGenCorrSetName += VKString::number(s_iGlobalGeneratingCorrSetUniqueIndex++);
	
	// NOTE: after set has been modified - it should be added to surfaces
	//		under a different name, and saved under a different name
	
	// It should be sufficient to optimize just on one surface
	// Can iterate several times, but might be slow and irrelevant
	
	// assert that every conformal map has same generating set
	assert(m_ConformalMaps.size()>0);	
	VKString genSetName = m_ConformalMaps[0]->m_GenSetName;

	for (int i=0; i<(int)m_ConformalMaps.size(); i++)
		assert(genSetName == m_ConformalMaps[i]->m_GenSetName);

	// create new generating set (that can be modified just for this map)
	SurfaceSampleSet * oldGenSet1 = GetSurface(0)->GetSampleSet(genSetName);	
	SurfaceSampleSet * oldGenSet2 = GetSurface(1)->GetSampleSet(genSetName);
	assert(oldGenSet1!=NULL && oldGenSet2!=NULL);

	SurfaceSampleSet * generatingSet1 = new SurfaceSampleSet();
	SurfaceSampleSet * generatingSet2 = new SurfaceSampleSet();;
	
	for (int i=0; i<oldGenSet1->NumSamples(); i++)
		generatingSet1->AddSample(oldGenSet1->GetSample(i));
	
	for (int i=0; i<oldGenSet2->NumSamples(); i++)
		generatingSet2->AddSample(oldGenSet2->GetSample(i));
	
	GetSurface(0)->AddSampleSet(m_sModifiedGenCorrSetName, generatingSet1);
	GetSurface(1)->AddSampleSet(m_sModifiedGenCorrSetName, generatingSet2);	
		
	// find search set
	SurfaceSampleSet * searchSet2 = GetSurface(1)->GetSampleSet(searchInSet);
	assert(searchSet2!=NULL);
	
	// for each correspondence find relevant conformal map
	// create new conformal maps (do not delete the old ones, they are in the Collection)	
	std::vector<std::vector<int> > s1ToRelMaps;
	std::vector<int> s1ToS2;	
	for (int i=0; i<generatingSet1->NumSamples(); i++)
	{
		s1ToRelMaps.push_back(std::vector<int>());
		s1ToS2.push_back(-1);
	}
		
	std::cout<<std::endl;
	for (int i=0; i<(int)m_ConformalMaps.size(); i++)
	{
		m_ConformalMaps[i]->GetGeneratingCorrespondences(&corrs, &sampleSet1, &sampleSet2);
		assert(sampleSet1==oldGenSet1);
		assert(sampleSet2==oldGenSet2);	
		for (int j=0; j<6; j+=2)
		{
			s1ToRelMaps[(*corrs)[j]].push_back(i);
			assert((s1ToS2[(*corrs)[j]]==-1) || (s1ToS2[(*corrs)[j]]==(*corrs)[j+1]));
			s1ToS2[(*corrs)[j]] = (*corrs)[j+1];
		}
		
		m_ConformalMaps[i] = new MapConformal(m_ConformalMaps[i]->m_pConfM1, 
											  m_ConformalMaps[i]->m_pConfM2, 
											  *corrs, m_sModifiedGenCorrSetName);
	}
	
	// go over all correspondences
	ClearConfidenceCalculator(confidenceName);
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("InitialEnergyCalculation");
	std::cout<<"Optimizing Correspondencees ("<<generatingSet1->NumSamples()<<"): "<<std::flush;
	double bestValue = GetMapConfidenceCalculator(confidenceName)->Confidence();
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("InitialEnergyCalculation", true);	
//	std::cout<<"\nBestValue="<<bestValue<<std::endl;
	std::vector<int> evalNhd;
	for (int i=0; i<generatingSet1->NumSamples(); i++)	
	{
		if (s1ToRelMaps[i].size()>0)
		{
			// find correspondence
			int s2 = s1ToS2[i];
			
			assert(s2!=-1);
			
			while (true)
			{
				// find sufficiently close samples (to s2)
				evalNhd.clear();
				GetNeighborhoodInSet(searchSet2, generatingSet2->GetSample(s2), evalNhd);
				assert(searchSet2->GetSample(evalNhd[0]).NearestVertex()==generatingSet2->GetSample(s2).NearestVertex());
				
				// evalNhd at 0 there should be current correspondence (i.e. 0-distance)
				int bestID = evalNhd[0];	
				for (int j=1; j<(int)evalNhd.size(); j++)
				{
					std::cout<<j<<" "<<std::flush;
					int currID = evalNhd[j];
					double currVal = UpdateAndGetValueOfConfMaps(s1ToRelMaps[i], generatingSet2,
																 s2, searchSet2->GetSample(currID),
																 confidenceName);
					
					if (currVal > bestValue)
					{
						bestValue = currVal;
						bestID = currID;
					}
				}
				
//				AnalysisStats::m_GlobalStats.m_Timing.startedProcess("BestEnergyCalculation");
//				std::cout<<"\nChangedCorr BestValue="<<bestValue<<std::endl;				
				// replace samples: generatingSet2[s2] = searchSet2[bestID]
				bestValue = UpdateAndGetValueOfConfMaps(s1ToRelMaps[i], generatingSet2,
														s2, searchSet2->GetSample(bestID),
														confidenceName);
//				std::cout<<"Changed Vertex: "<<generatingSet1->GetSample(i).NearestVertex();
//				std::cout<<" -> "<<generatingSet2->GetSample(s2).NearestVertex()<<std::endl;
//				std::cout<<"BestID="<<bestID<<" BestVertex="<<searchSet2->GetSample(bestID).NearestVertex()<<std::endl;
				
//				AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("BestEnergyCalculation", true);
				if (bestID==evalNhd[0])
					break;
			}
		}
//		break;
	}
	std::cout<<std::endl;
}

double MapMultiConformal::UpdateAndGetValueOfConfMaps(std::vector<int> & updateConfMaps, 
													  SurfaceSampleSet * generatingSet2,
													  int s2, const SurfaceSample & newS2,
													  const VKString &confidenceName)
{
	ClearConfidenceCalculator(confidenceName);
	generatingSet2->ReplaceSample(s2, newS2);
	for (int i=0; i<(int)updateConfMaps.size(); i++)
		m_ConformalMaps[updateConfMaps[i]]->UpdateGeneratingCorrsChanged();
	return GetMapConfidenceCalculator(confidenceName)->Confidence();
}

void MapMultiConformal::GetNeighborhoodInSet(SurfaceSampleSet * searchInSet, 
											 const SurfaceSample & nearMe,
											 std::vector<int> & searchedIDs)
{
	// NOTE: serchedIDs should have the lowest-distance element 
	//			to have zero distance and be first (i.e. genSet <= searchSet)
	searchedIDs.clear();
	double distThresh = 0.25;	// TODO: read value from paramfile
	SampledSurface * surface = GetSurface(nearMe);
	DistanceOnTheFly & onTheFlyDist = *(surface->GetOnTheFlyDistanceMetric(distThresh, "default", -1));
	double maxDistance = surface->AdjustedRadius(distThresh);		
	
	bool addedFirst = false;
	for (int i=0; i<searchInSet->NumSamples(); i++)
	{
		double dist = onTheFlyDist.Distance(nearMe, searchInSet->GetSample(i));
		if (dist < maxDistance)
		{
			if (dist==0)
			{
				assert(!addedFirst);
				addedFirst = true;				
				searchedIDs.insert(searchedIDs.begin(), i);
			}
			else
				searchedIDs.push_back(i);
		}
	}
	assert(addedFirst);
}

void MapMultiConformal::PruneConfMaps(const VKString & confidenceName, int maxMaps)
{
	std::vector<int> * corrs;
	SurfaceSampleSet * sampleSet1;
	SurfaceSampleSet * sampleSet2;

	std::map<int, int> corrToNumMaps;
	
	for (int i=0; i<(int)m_ConformalMaps.size(); i++)
	{
		m_ConformalMaps[i]->GetGeneratingCorrespondences(&corrs, &sampleSet1, &sampleSet2);
		assert(corrs->size()==6);		
		if (corrToNumMaps.find((*corrs)[0])==corrToNumMaps.end())
			corrToNumMaps[(*corrs)[0]] = 0;
		if (corrToNumMaps.find((*corrs)[1])==corrToNumMaps.end())
			corrToNumMaps[(*corrs)[2]] = 0;
		if (corrToNumMaps.find((*corrs)[2])==corrToNumMaps.end())
			corrToNumMaps[(*corrs)[4]] = 0;
		
		corrToNumMaps[(*corrs)[0]]++;
		corrToNumMaps[(*corrs)[2]]++;
		corrToNumMaps[(*corrs)[4]]++;		
	}
	
	std::cout<<"Pruning Conformal Maps (Trg="<<maxMaps<<"): "<<std::flush;
	while ((int)m_ConformalMaps.size() > maxMaps)
	{
		double bestEnergy = -FLT_MAX;
		int removeConfID=-1;
		
		for (int i=0; i<(int)m_ConformalMaps.size(); i++)
		{
			m_ConformalMaps[i]->GetGeneratingCorrespondences(&corrs, &sampleSet1, &sampleSet2);
			bool canRemove = (corrToNumMaps[(*corrs)[0]]>1);
			canRemove = canRemove && (corrToNumMaps[(*corrs)[2]]>1);
			canRemove = canRemove && (corrToNumMaps[(*corrs)[4]]>1);
			if (!canRemove)
				continue;
			
			ClearConfidenceCalculator(confidenceName);
			MapConformal * confRemove = m_ConformalMaps[i];
			m_ConformalMaps.erase(m_ConformalMaps.begin()+i);
			double currEnergy = GetMapConfidenceCalculator(confidenceName)->Confidence();
			if (bestEnergy < currEnergy)
			{
				bestEnergy = currEnergy;
				removeConfID = i;
			}
			m_ConformalMaps.insert(m_ConformalMaps.begin()+i, confRemove);
		}
		std::cout<<"[N="<<m_ConformalMaps.size()<<" E="<<bestEnergy<<"] "<<std::flush;
		assert(removeConfID!=-1);
		
		m_ConformalMaps[removeConfID]->GetGeneratingCorrespondences(&corrs, &sampleSet1, &sampleSet2);
		corrToNumMaps[(*corrs)[0]]--;
		corrToNumMaps[(*corrs)[2]]--;
		corrToNumMaps[(*corrs)[4]]--;		
		
		m_ConformalMaps.erase(m_ConformalMaps.begin() + removeConfID);
	}
	std::cout<<std::endl;
}

void MapMultiConformal::Draw(AnalysisWindow * window, ParamParser * params,
							 const VKString & renderingParams,
							 const VKString & surfaceName)
{	
	
	DrawMultiCorrespondenceForSelected(window, params, renderingParams, surfaceName);
	bool valid;
	if (params->GetStrValues(renderingParams, "RenderCorrCoarse", valid).contains(GetSurfaceMapType()))
	{
		SurfaceSampleSet * genSet1 = NULL;
		SurfaceSampleSet * genSet2 = NULL;
		std::map<int, int> genCorrs;
		for (int i=0; i<(int) m_ConformalMaps.size(); i++)
		{
			SurfaceSampleSet * set1=NULL, * set2=NULL;
			std::vector<int> * corrs=NULL;
			m_ConformalMaps[i]->GetGeneratingCorrespondences(&corrs, &set1, &set2);
			assert(set1!=NULL && set2!=NULL && corrs!=NULL);
			if (i==0)
			{
				genSet1=set1;
				genSet2=set2;			
			}
			assert(genSet1==set1 && genSet2==set2);
			assert((int)corrs->size()==6);
			for (int c=0; c<6; c+=2)
				genCorrs[(*corrs)[c+0]] = (*corrs)[c+1];
		}
		
		bool mapRadical = ((int)genCorrs.size() <= 10);
		
		for (int i=0; i<(int) m_ConformalMaps.size(); i++)
			m_ConformalMaps[i]->DrawGenerators(window, params, renderingParams, surfaceName, 
											   mapRadical ? (3.) : 1, .8, 0, 0, true, false);

		glEnable(GL_LIGHTING);
		int colorID=0;
		for (std::map<int, int>::iterator iter = genCorrs.begin(); iter != genCorrs.end(); iter++)
		{
			double radius1 = GetSurface(0)->GetStandardRadius(params, renderingParams, surfaceName) * 3;
			double radius2 = GetSurface(1)->GetStandardRadius(params, renderingParams, surfaceName) * 3;
			
			static GLfloat material[4];	
			//window->MapIntToRadicalColor(colorID++, material);
			double r, g, b;
			if (mapRadical)
				window->MapIntToRadicalColor(colorID++, r, g, b);
			else
			{
				radius1 /= 3.;
				radius2 /= 3.;
				window->MapIntToColor(colorID++, r, g, b);
			}
			material[0] = r;		material[1] = g;		material[2] = b;		material[3] = 1;
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material);
					
			GetSurface(0)->LoadCamera();
			R3Point p1 = GetSurface(0)->GetDrawablePosition(genSet1->GetSample(iter->first), 
															params, renderingParams, surfaceName);				
			R3Sphere(p1, radius1).Draw();
			GetSurface(1)->LoadCamera();				
			R3Point p2 = GetSurface(1)->GetDrawablePosition(genSet2->GetSample(iter->second), 
															params, renderingParams, surfaceName);				
			R3Sphere(p2, radius2).Draw();
			
		}
		
		//bool valid;
		//bool renderCoarse = params->GetStrValues(renderingParams, "RenderCorrCoarse", valid).contains(GetSurfaceMapType());
		//if (renderCoarse && GetCurrentCoarseMap()!=NULL)
		//	GetCurrentCoarseMap()->Draw(window, params, renderingParams, surfaceName);
	}
}

int MapMultiConformal::GetNumConformalMaps()
{
	return (int)m_ConformalMaps.size();
}





