#include "MapGeoFeature.h"
#include "MapCoarse.h"
#include "MapConformal.h"
#include "MapConfidenceCompareToTruth.h"
#include "Sortable.h"
#include "SurfaceDistance.h"
#include <algorithm>

SurfaceSampleSet * MapCoarse::GetSetFeaturePoints(SampledSurface * M, 
												  const VKString & verticesToFeatPnts)
{
	std::vector<int> featurePointToVertex;
	FillFeatureToVertex(verticesToFeatPnts, featurePointToVertex);
	return new SurfaceSampleSet(M->GetMesh(), featurePointToVertex);
}

void MapCoarse::PreinitializeAllToNULL()
{
	m_pRenderAnotherCoarse = NULL;
	m_pSearchNodeArray1 = NULL;
	m_pSearchNodeArray2 = NULL;
	m_iSearchArraySize1 = -1;
	m_iSearchArraySize2 = -1;	
	
	m_fEpsilon = 0.00001;
	m_sTruthConfidence = "";
	m_pSet1 = NULL;
	m_pSet2 = NULL;
	m_VotedCorrMatrix = NULL;
	
	m_CorrColorR=0;
	m_CorrColorG=1;
	m_CorrColorB=0;		
	m_pTruthMapForDrawing = NULL;
}

MapCoarse::MapCoarse(SampledSurface * M1, SampledSurface * M2)
: SurfaceMap(M1, M2)
{
	PreinitializeAllToNULL();
}

MapCoarse::MapCoarse(SampledSurface * M1, SampledSurface * M2, 
					 SurfaceSampleSet * sampleSet1, SurfaceSampleSet * sampleSet2)
: SurfaceMap(M1, M2)
{
	PreinitializeAllToNULL();
	InitializeSets(sampleSet1, sampleSet2);
}

MapCoarse::MapCoarse(SampledSurface * M1, SampledSurface * M2, 
					 const VKString & sampleSet1, const VKString & sampleSet2)
: SurfaceMap(M1, M2)
{	
	PreinitializeAllToNULL();
	InitializeSets(M1->GetSampleSet(sampleSet1), M2->GetSampleSet(sampleSet2));
}

MapCoarse::MapCoarse(SurfaceSampleSet * sampleSet1, SurfaceSampleSet * sampleSet2,
					 SurfaceMap * surfaceMap, FineToCoarseGeneration method,
					 const VKString & distanceMetric1, const VKString & distanceMetric2)

: SurfaceMap(surfaceMap->GetSurface(0), surfaceMap->GetSurface(1))
{
	PreinitializeAllToNULL();	
	SurfaceSampleSet * sampleSetFrom = new SurfaceSampleSet();
	SurfaceSampleSet * sampleSetTo = new SurfaceSampleSet();
	std::vector<int> correspondenceMap;
	std::vector<double> errs;
	
	switch(method)
	{	
		case F2C_FORWARD_ONLY:
		case F2C_FORWARD_BACKWARD:
		case F2C_BACKWARD_ONLY:			// just do forward/backward map for some subset of points
			CreateSimpleCoarse(sampleSetFrom, sampleSetTo, correspondenceMap, 
							   sampleSet1, sampleSet2, surfaceMap, method);
			break;
		case F2C_MUTUALLY_CLOSEST_NEIGHBORS:	//choose mutually-closest neighbors, or error based on some metric
		case F2C_MCN_ONLY_ON_SURF2_FORWARD:
		case F2C_MCN_ONLY_ON_SURF1_BACKWARD:
			CreateProcessedCoarse(sampleSetFrom, sampleSetTo, correspondenceMap, errs, 
								  sampleSet1, sampleSet2, surfaceMap, method,
								  GetSurface(0)->GetDistanceMetric(distanceMetric1),
								  GetSurface(1)->GetDistanceMetric(distanceMetric2));
			break;			
		default:
			std::cout<<"[ERROR] Unknown fine-to-coarse method"<<std::endl;
			assert(false);
	}
	
	InitializeSets(sampleSetFrom, sampleSetTo);	
	SetFinalCorrMap(correspondenceMap);
}

void MapCoarse::AssociateCoarseToRender(MapCoarse * mapCoarse)
{
	m_pRenderAnotherCoarse = mapCoarse;
}

void MapCoarse::CreateSimpleCoarse(SurfaceSampleSet * sampleSetFrom, 
								   SurfaceSampleSet * sampleSetTo,
								   std::vector<int> &correspondenceMap, 
								   const SurfaceSampleSet * sampleSet1, 
								   const SurfaceSampleSet * sampleSet2,
								   SurfaceMap * surfaceMap, FineToCoarseGeneration method)
{
	TimeProfiler coarseMapProgress;
	if (method==F2C_FORWARD_ONLY || method==F2C_FORWARD_BACKWARD)
	{
		std::cout<<"Creating Simple Forward coarse map"<<std::endl;
		assert(sampleSet1!=NULL);
		for (int i=0; i<sampleSet1->NumSamples(); i++)
		{
			if (i % 100 == 0)
				coarseMapProgress.WriteProgress("CreateCoarseMapForward", i, sampleSet1->NumSamples());
			const SurfaceSample & s1 = sampleSet1->GetSample(i);
			const SurfaceSample & s2 = surfaceMap->ForwardMap(s1);	
			if (!s2.Invalid())
			{
				int id1, id2;
				if (!sampleSetFrom->ContainsExact(s1, &id1))
					id1 = sampleSetFrom->AddSample(s1);
				if (!sampleSetTo->ContainsExact(s2, &id2))
					id2 = sampleSetTo->AddSample(s2);				
				while ((int)correspondenceMap.size()<=id1)
					correspondenceMap.push_back(-1);
				correspondenceMap[id1]=id2;
			}
		}
		std::cout<<std::endl;	
	}
	
	if (method==F2C_BACKWARD_ONLY || method==F2C_FORWARD_BACKWARD)
	{
		std::cout<<"Creating Simple Forward coarse map"<<std::endl;
		assert(sampleSet2!=NULL);		
		for (int i=0; i<sampleSet2->NumSamples(); i++)
		{
			if (i % 100 == 0)
				coarseMapProgress.WriteProgress("CreateCoarseMapBackward", i, sampleSet2->NumSamples());			
			const SurfaceSample & s2 = sampleSet2->GetSample(i);
			const SurfaceSample & s1 = surfaceMap->InverseMap(s2);
			if (!s1.Invalid())
			{
				int id1, id2;
				if (!sampleSetFrom->ContainsExact(s1, &id1))
					id1 = sampleSetFrom->AddSample(s1);
				if (!sampleSetTo->ContainsExact(s2, &id2))
					id2 = sampleSetTo->AddSample(s2);				
				while ((int)correspondenceMap.size()<=id1)
					correspondenceMap.push_back(-1);
				correspondenceMap[id1]=id2;
			}
		}
		std::cout<<std::endl;
	}	
}

void MapCoarse::CreateProcessedCoarse(SurfaceSampleSet * sampleSetFrom, 
									  SurfaceSampleSet * sampleSetTo,
									  std::vector<int> &correspondenceMap, 
									  std::vector<double> &dists, 
									  const SurfaceSampleSet * sampleSet1, 
									  const SurfaceSampleSet * sampleSet2,
									  SurfaceMap * surfaceMap, FineToCoarseGeneration method,
									  SurfaceDistance * dist1, SurfaceDistance * dist2)
{
	assert(method==F2C_MUTUALLY_CLOSEST_NEIGHBORS
		   || method==F2C_MCN_ONLY_ON_SURF2_FORWARD
		   || method==F2C_MCN_ONLY_ON_SURF1_BACKWARD);
	assert(dist1!=NULL && dist2!=NULL);

	LinAlgMatrixReal distMatrix1To2(sampleSet1->NumSamples(), sampleSet2->NumSamples());
	LinAlgMatrixReal distMatrix2To1(sampleSet2->NumSamples(), sampleSet1->NumSamples());
	

	FillDistanceMatrix(sampleSet1, sampleSet2, distMatrix1To2, 
					   distMatrix2To1, dist1, dist2, surfaceMap, method);	
	
	std::vector<int> sampleMap;
	AssignMutuallyClosest(sampleSet1, sampleSet2, distMatrix1To2, distMatrix2To1, correspondenceMap, dists);
	
	for (int i=0; i<sampleSet1->NumSamples(); i++)
		sampleSetFrom->AddSample(sampleSet1->GetSample(i));
	for (int i=0; i<sampleSet2->NumSamples(); i++)
		sampleSetTo->AddSample(sampleSet2->GetSample(i));
}

void MapCoarse::FillDistanceMatrix(const SurfaceSampleSet * set1, const SurfaceSampleSet * set2, 
								   LinAlgMatrixReal & distMatrix1To2, LinAlgMatrixReal & distMatrix2To1,
								   SurfaceDistance * dist1, SurfaceDistance * dist2, SurfaceMap * surfaceMap,
								   FineToCoarseGeneration f2cmethod)

{
	std::vector<SurfaceSample> samplesOnS2;
	std::vector<SurfaceSample> samplesOnS1;

	if (f2cmethod==F2C_MUTUALLY_CLOSEST_NEIGHBORS || f2cmethod==F2C_MCN_ONLY_ON_SURF2_FORWARD)
	{
		for (int i=0; i<set1->NumSamples(); i++)
		{
			SurfaceSample sOnS2 = surfaceMap->ForwardMap(set1->GetSample(i));
			samplesOnS2.push_back(sOnS2);
		}
	}

	if (f2cmethod==F2C_MUTUALLY_CLOSEST_NEIGHBORS || f2cmethod==F2C_MCN_ONLY_ON_SURF1_BACKWARD)
	{
		for (int i=0; i<set2->NumSamples(); i++)
		{
			SurfaceSample invS2 = surfaceMap->InverseMap(set2->GetSample(i));
			samplesOnS1.push_back(invS2);
		}
	}

	double scale1 = sqrt(GetSurface(0)->Area());
	double scale2 = sqrt(GetSurface(1)->Area());

	for (int i=0; i<set1->NumSamples(); i++)	// fill geodesic distances
	{
		for (int j=0; j<set2->NumSamples(); j++)
		{
			if (f2cmethod==F2C_MUTUALLY_CLOSEST_NEIGHBORS 
				|| f2cmethod==F2C_MCN_ONLY_ON_SURF2_FORWARD)
			{
				assert(i < (int)samplesOnS2.size());
				assert(j < set2->NumSamples());
				if (samplesOnS2[i].Invalid() || set2->GetSample(j).Invalid())
					distMatrix1To2(i, j) = (double)FLT_MAX;
				else
					distMatrix1To2(i, j) = dist2->Distance(samplesOnS2[i], set2->GetSample(j)) / scale2;
				
			}
			else 
			{
				assert(f2cmethod==F2C_MCN_ONLY_ON_SURF1_BACKWARD);
				assert(i < set1->NumSamples());
				assert(j < (int)samplesOnS1.size());
				if (samplesOnS1[j].Invalid() || set1->GetSample(i).Invalid())
					distMatrix1To2(i, j) = (double)FLT_MAX;
				else					
					distMatrix1To2(i, j) = dist1->Distance(set1->GetSample(i), samplesOnS1[j]) / scale1;
			}

			if (f2cmethod==F2C_MUTUALLY_CLOSEST_NEIGHBORS 
				|| f2cmethod==F2C_MCN_ONLY_ON_SURF1_BACKWARD)
			{
				assert(j < (int)samplesOnS1.size());
				assert(i < set1->NumSamples());
				
				if (samplesOnS1[j].Invalid() || set1->GetSample(i).Invalid())
					distMatrix2To1(j, i) = (double)FLT_MAX;
				else
					distMatrix2To1(j, i) = dist1->Distance(samplesOnS1[j], set1->GetSample(i)) / scale1;
			}
			else 
			{
				assert(i < (int)samplesOnS2.size());
				assert(j < set2->NumSamples());
				assert(f2cmethod==F2C_MCN_ONLY_ON_SURF2_FORWARD);
				if (samplesOnS2[i].Invalid() || set2->GetSample(j).Invalid())
					distMatrix2To1(j, i) = (double)FLT_MAX;
				else
					distMatrix2To1(j, i) = dist2->Distance(set2->GetSample(j), samplesOnS2[i]) / scale2;
			}
		}
	}
}

void MapCoarse::ConformalEuclideanMCN(const SurfaceSampleSet * set1, const SurfaceSampleSet * set2,
									  MapConformal * map1,
									  std::vector<int> & sampleMap, std::vector<double> & corrDist)
{
	assert(GetSurface(0)->GetSurfaceType()=="SurfaceMidEdgeConf");
	assert(GetSurface(1)->GetSurfaceType()=="SurfaceMidEdgeConf");	
	
	SurfaceMidEdgeConf * mc1 = (SurfaceMidEdgeConf*)GetSurface(0);
	SurfaceMidEdgeConf * mc2 = (SurfaceMidEdgeConf*)GetSurface(1);
	
	MobiusTransformation m1 = mc1->GetCurrentTransform();
	MobiusTransformation m2 = mc2->GetCurrentTransform();
	if (mc1->GetCurrentTransform()!=map1->GetMobiusForSurface(0))
		mc1->Transform(map1->GetMobiusForSurface(0));
	if (mc2->GetCurrentTransform()!=map1->GetMobiusForSurface(1))
		mc2->Transform(map1->GetMobiusForSurface(1));
	
	sampleMap.clear();
	corrDist.clear();

	
//	// BELOW IS USING KD-STRUCTURE FOR MUTUALLY CLOSEST NEIGHBORS
//	// build kd-structure for 2D positions of points
//		// initialize data holders
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CastVoteMCNCalcPrepare");	
//	RNArray<FlatSearchNode*> samples1;	
//	RNArray<FlatSearchNode*> samples2;
//	if (m_pSearchNodeArray1==NULL)
//	{
//		m_pSearchNodeArray1 = new FlatSearchNode[set1->NumSamples()];
//		m_iSearchArraySize1 = set1->NumSamples();
//		for (int i=0; i<set1->NumSamples(); i++)
//			m_pSearchNodeArray1[i].m_iID = i;
//	}
//	if (m_pSearchNodeArray2==NULL)
//	{
//		m_pSearchNodeArray2 = new FlatSearchNode[set2->NumSamples()];
//		m_iSearchArraySize2 = set2->NumSamples();
//		for (int i=0; i<set2->NumSamples(); i++)
//			m_pSearchNodeArray2[i].m_iID = i;
//	}	
//	assert(m_iSearchArraySize1==set1->NumSamples());
//	assert(m_iSearchArraySize2==set2->NumSamples());	
//	
//		// for each set 1 -> put in rnarray -> kdtree
//	for (int i=0; i<set1->NumSamples(); i++)
//	{
//		LinAlgComplex v = mc1->GetConfCoord(set1->GetSample(i));
//		//LinAlgComplex v = m1.Transform(m_pOriginalComplexCoords1[i]);
//		m_pSearchNodeArray1[i].m_Pnt = R2Point(v.r, v.i);
//		samples1.InsertTail(&(m_pSearchNodeArray1[i]));
//	}
//	for (int i=0; i<set2->NumSamples(); i++)
//	{
//		LinAlgComplex v = mc2->GetConfCoord(set2->GetSample(i));
//		//LinAlgComplex v = m2.Transform(m_pOriginalComplexCoords2[i]);
//		m_pSearchNodeArray2[i].m_Pnt = R2Point(v.r, v.i);
//		samples2.InsertTail(&(m_pSearchNodeArray2[i]));
//	}
//	
//		// create kd trees
//	R2Kdtree<FlatSearchNode*> kdtree1(samples1, &global_FlatNodePosition, NULL);
//	R2Kdtree<FlatSearchNode*> kdtree2(samples2, &global_FlatNodePosition, NULL);	
//	
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("CastVoteMCNCalcPrepare");		
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CastVoteMCNCalcFind");			
//	// for each set 1  -> query w/ 2D location of s2, for s2 find closest
//		// record nearest neighbor in corrDist, sampleMap
//	for (int i=0; i<set1->NumSamples(); i++)
//	{
//		FlatSearchNode * closestInS2 = kdtree2.FindClosest(m_pSearchNodeArray1[i].m_Pnt);
//		// note: can cache closest in s1
//		FlatSearchNode * closestInS1 = kdtree1.FindClosest(closestInS2->m_Pnt);
//		if (closestInS1->m_iID==i)	// mutually closest
//		{
//			corrDist.push_back((m_pSearchNodeArray1[i].m_Pnt-closestInS2->m_Pnt).Length());
//			sampleMap.push_back(closestInS2->m_iID);
//		}
//		else
//		{
//			corrDist.push_back(0);
//			sampleMap.push_back(-1);
//		}
//	}
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("CastVoteMCNCalcFind");				
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CastVoteMCNCalc");	
	
	// ANOTHER EXHAUSTIVE SEARCH VERSION
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CastVoteMCNCalcPrepare");					
	for (int i=0; i<set1->NumSamples(); i++)
		m_pTransformedComplexCoords1[i] = m1.Transform(m_pOriginalComplexCoords1[i]);
	for (int i=0; i<set2->NumSamples(); i++)
		m_pTransformedComplexCoords2[i] = m2.Transform(m_pOriginalComplexCoords2[i]);
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("CastVoteMCNCalcPrepare");						
//	
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CastVoteMCNCalcFind");					
	std::vector<int> nearestNeihbor1;
	std::vector<int> nearestNeihbor2;
	std::vector<double> nearestNeihborDist2;
	for (int i=0; i<set1->NumSamples(); i++)
	{
		double minima = 0; 
		int minimaID = -1;
		LinAlgComplex & v1 = m_pTransformedComplexCoords1[i];
		for (int j=0; j<set2->NumSamples(); j++)
		{
			LinAlgComplex v2 = m_pTransformedComplexCoords2[j];
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
	
	for (int i=0; i<(int)nearestNeihbor1.size(); i++)
	{
		int s1 = i;
		int s2 = nearestNeihbor1[s1];
		if (nearestNeihbor2[s2]==s1)
		{
			corrDist.push_back(sqrt(nearestNeihborDist2[s2]));
			sampleMap.push_back(s2);
		}
		else
		{
			corrDist.push_back(0);
			sampleMap.push_back(-1);
		}
	}
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("CastVoteMCNCalcFind");						
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("CastVoteMCNCalc");		
}

void MapCoarse::AssignMutuallyClosest(const SurfaceSampleSet * set1, const SurfaceSampleSet * set2, 
									  LinAlgMatrixReal & distMatrix1To2, LinAlgMatrixReal & distMatrix2To1, 
									  std::vector<int> & sampleMap, std::vector<double> & corrError)
{
	std::vector<int> set1Closest;
	std::vector<double> dist1Closest;
	std::vector<int> set2Closest;
	std::vector<double> dist2Closest;
	// cast point from mesh 1 to mesh 2, find closest sample on mesh 2 (on mesh or on conformal space) 
	for (int i=0; i<set1->NumSamples(); i++)
	{	
		double minValue=-1;
		int minCorr=-1;
		for (int j=0; j<set2->NumSamples(); j++)
		{
			double val = distMatrix1To2(i, j);
			if (minCorr<0 || minValue>val)
			{
				minValue = val;
				minCorr = j;
			}
		}
		set1Closest.push_back(minCorr);
		dist1Closest.push_back(minValue);
	}
	
	// cast point from mesh 2 to mesh 1, find closest sample on mesh 1 (on mesh or on conformal space) 
	for (int i=0; i<set2->NumSamples(); i++)
	{	
		double minValue=-1;
		int minCorr=-1;
		for (int j=0; j<set1->NumSamples(); j++)
		{
			double val = distMatrix2To1(i, j);
			if (minCorr<0 || minValue>val)
			{
				minValue = val;
				minCorr = j;
			}
		}
		set2Closest.push_back(minCorr);
		dist2Closest.push_back(minValue);
	}
	
	for (int i=0; i<(int)set1Closest.size(); i++)
	{
//		if (set1->GetSample(i).NearestVertex()==10026)
//		{
//			std::cout<<"CN["<<set1->GetSample(i).NearestVertex()<<"] = ";
//			std::cout<<set2->GetSample(set1Closest[i]).NearestVertex();
//			std::cout<<" and CN[CN[i]]="<<set1->GetSample(set2Closest[set1Closest[i]]).NearestVertex()<<std::endl;
//		}
			
		if (set2Closest[set1Closest[i]]==i)	// mutually closest neighbor
		{
			sampleMap.push_back(set1Closest[i]);
			corrError.push_back(dist1Closest[i]+dist2Closest[set1Closest[i]]);	// error=d1+d2
		}
		else	// not mutually closest - do not assign
		{
			sampleMap.push_back(-1);
			corrError.push_back(-1);
		}
	}
}

MapCoarse::~MapCoarse()
{
	if (m_VotedCorrMatrix!=NULL)
		delete m_VotedCorrMatrix;
}

void MapCoarse::InitializeVoting()
{
	assert(GetSurface(0)->GetSurfaceType()=="SurfaceMidEdgeConf");
	assert(GetSurface(1)->GetSurfaceType()=="SurfaceMidEdgeConf");	
	
	SurfaceMidEdgeConf * mconf1= (SurfaceMidEdgeConf*)GetSurface(0);
	SurfaceMidEdgeConf * mconf2= (SurfaceMidEdgeConf*)GetSurface(1);	
	
	m_pOriginalComplexCoords1 = new LinAlgComplex[m_pSet1->NumSamples()];
	m_pOriginalComplexCoords2 = new LinAlgComplex[m_pSet2->NumSamples()];
	m_pTransformedComplexCoords1 = new LinAlgComplex[m_pSet1->NumSamples()];
	m_pTransformedComplexCoords2 = new LinAlgComplex[m_pSet2->NumSamples()];
	
	for (int i=0; i<m_pSet1->NumSamples(); i++)
		m_pOriginalComplexCoords1[i] = mconf1->GetConfCoordOriginal(m_pSet1->GetSample(i));

	for (int i=0; i<m_pSet2->NumSamples(); i++)	
		m_pOriginalComplexCoords2[i] = mconf2->GetConfCoordOriginal(m_pSet2->GetSample(i));		
	
	m_VotedCorrMatrix = new LinAlgMatrixReal(m_pSet1->NumSamples(), m_pSet2->NumSamples());
}

void MapCoarse::InitializeSets(SurfaceSampleSet * sampleSet1, SurfaceSampleSet * sampleSet2)
{
	m_sTruthConfidence="";
	m_OrderOfDomainRendering.clear();
	m_CorrColorR=0;
	m_CorrColorG=1;
	m_CorrColorB=0;	
	m_pTruthMapForDrawing = NULL;
	
	assert(sampleSet1!=NULL);
	assert(sampleSet2!=NULL);

	m_pSet1 = sampleSet1;
	m_pSet2 = sampleSet2;
	m_VotedCorrMatrix = NULL;
	
	m_mVertexIDToSetID1.clear();
	m_mVertexIDToSetID2.clear();
	for (int i=0; i<m_pSet1->NumSamples(); i++)
		m_mVertexIDToSetID1[m_pSet1->GetSample(i).NearestVertex()] = i;
	
	for (int i=0; i<m_pSet2->NumSamples(); i++)
		m_mVertexIDToSetID2[m_pSet2->GetSample(i).NearestVertex()] = i;
	
	ClearAllCorrespondences();	
}

void MapCoarse::ClearAllCorrespondences()
{
	ClearCache();
	if (m_VotedCorrMatrix!=NULL)
	{
		for (int i=0; i<m_pSet1->NumSamples(); i++)
		for (int j=0; j<m_pSet2->NumSamples(); j++)
			(*m_VotedCorrMatrix)(i,j) = 0;
	}
	
	m_mFinalCorrMapM1ToM2.clear();
	m_mFinalCorrValue.clear();
}

void MapCoarse::SetFinalCorrMap(const std::vector<int> & M1ToM2)
{
	assert((int)M1ToM2.size()==m_pSet1->NumSamples());
	m_mFinalCorrValue.clear();	   
	m_mFinalCorrMapM1ToM2 = M1ToM2;
	for (int i=0; i<(int)m_mFinalCorrMapM1ToM2.size(); i++)
		 m_mFinalCorrValue.push_back(1.);
}

double MapCoarse::CastVote(SurfaceMap * lowDimMap, FineToCoarseGeneration voteForMCN,
						   const VKString & confidenceName)
{
	assert(lowDimMap->GetSurfaceMapType()=="MapConformal");
	MapConformal * confMap = (MapConformal*)lowDimMap;
	std::vector<int> sampleMap;
	std::vector<double> corrDist;
	// find MCN in the set + confidence value
	if (confidenceName=="none" || confidenceName=="")
	{
		assert(voteForMCN==F2C_MCN_CONFORMAL_EUCLIDEAN);
		ConformalEuclideanMCN(m_pSet1, m_pSet2, confMap, sampleMap, corrDist);
	}
	else
	{
		assert(false);
//		conf = lowDimMap->GetMapConfidenceCalculator(confidenceName)->Confidence();
	}

	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CastVoteCasting");	
	// cast votes for MCNs
	double cumErr=0;
	double norm=0;
	for (int i=0; i<(int)corrDist.size(); i++)
	{
		if (sampleMap[i]!=-1)
		{
			cumErr+=corrDist[i];
			norm += 1.;
		}
	}
	
	double fracOverlap = (double)norm / (double)corrDist.size();
	
	if (fracOverlap > .4)
	{
		//std::cout<<"fracOverlap="<<fracOverlap<<std::endl;
		AnalysisStats::m_GlobalStats.Increase("Voting:AcceptedVotes");
		//double voteValue = 1. / (m_fEpsilon + (cumErr/norm));
		double voteValue = fracOverlap;
		for (int i=0; i<(int)sampleMap.size(); i++)
			if (sampleMap[i]!=-1)
				CastVote(voteValue, i, sampleMap[i]);
		return voteValue;
	}
	else
	{
		AnalysisStats::m_GlobalStats.Increase("Voting:RejectedVotes");				
		return 0;
	}
}

void MapCoarse::CastVote(double confValue, const SurfaceSample & s1, const SurfaceSample & s2)
{
	assert(m_mVertexIDToSetID1.find(s1.NearestVertex())!=m_mVertexIDToSetID1.end());
	assert(m_mVertexIDToSetID2.find(s2.NearestVertex())!=m_mVertexIDToSetID2.end());
	CastVote(confValue, m_mVertexIDToSetID1[s1.NearestVertex()], 
			 m_mVertexIDToSetID2[s2.NearestVertex()]);
}

void MapCoarse::CastVote(double confValue, int s1, int s2)
{
	
	(*m_VotedCorrMatrix)(s1, s2) += confValue;
//	assert(m_VotedCorrMatrix!=NULL);
//	switch(castType)
//	{
//		case VOTE_CAST_SUM:
//			(*m_VotedCorrMatrix)(s1, s2) += confValue;
//			break;
//		case VOTE_CAST_MIN_Val:
//			(*m_VotedCorrMatrix)(s1, s2) = vkMax(confValue, (*m_VotedCorrMatrix)(s1, s2));
//			break;
//		default:
//			assert(false);
//	}	
}

void MapCoarse::AssignGreedyCorr(double takeBestMaps, bool furtherPairsSearch)
{
//	std::cout<<"Assigning Greedy Corrs"<<std::endl;
	assert(m_VotedCorrMatrix!=NULL);
	for (int i=0; i<m_pSet1->NumSamples(); i++)
	{
		m_mFinalCorrMapM1ToM2.push_back(-1);
		m_mFinalCorrValue.push_back(-1);
	}

	LinAlgMatrixReal copy(*m_VotedCorrMatrix);

	int row, col;	
	double maxValue = copy.GetMax(&row, &col);
	if (maxValue==0)
	{
		std::cout<<"[WARNING] Cannot assign correspondences - no successful votes"<<std::endl;
		m_mFinalCorrMapM1ToM2[0] = 0;	// assign any corr - to avoid crashes
		m_mFinalCorrValue[0] = 1.;
		return;
	}

	std::vector<Sortable> sortedCorrs;
	while (maxValue!=0)
	{
		for (int i=0; i<copy.Rows(); i++)
			copy(i, col) = 0;

		for (int i=0; i<copy.Cols(); i++)
			copy(row, i) = 0;

		m_mFinalCorrMapM1ToM2[row] = col;
		m_mFinalCorrValue[row] = maxValue;
//		std::cout<<"AssignGreedy: "<<row<<" -> "<<col<<" val = "<<maxValue<<std::endl;
		maxValue = copy.GetMax(&row, &col);
		sortedCorrs.push_back(Sortable(maxValue, NULL, row));
	}
	
//	std::cout<<"Assigning Greedy Corrs - done matrix analysis"<<std::endl;	
	std::sort(sortedCorrs.begin(), sortedCorrs.end());

	// find good set (top N samples), find potential candidates for final pair refinement
	int numBestMaps = (int)(m_pSet1->NumSamples() * takeBestMaps);
	std::vector<int> goodSetFrom;
	std::vector<int> goodSetTo;
	std::vector<int> potentialCorrFrom;
	std::vector<int> potentialCorrTo;
	for (int i=0; i < (int)sortedCorrs.size(); i++)
	{
		int corrFromID = sortedCorrs[(int)sortedCorrs.size() - i - 1].id;
		int corrToID = m_mFinalCorrMapM1ToM2[corrFromID];
		double corrVal = sortedCorrs[(int)sortedCorrs.size() - i - 1].value;		
		
		if (corrVal>0 && corrToID!=-1)
		{
			if ((int)goodSetTo.size() < numBestMaps)	// adding good corrs
			{
				goodSetFrom.push_back(corrFromID);
				goodSetTo.push_back(corrToID);				
			}
			else	// adding potential further pairs
			{
				potentialCorrFrom.push_back(corrFromID);
				potentialCorrTo.push_back(corrToID);								
			}
		}
	}
	
//	std::cout<<"Assigning Greedy Corrs - done greedy best"<<std::endl;		
	// further pairs refinement (end of the section 8 in 'Mobius Voting' paper)
	assert((int)goodSetTo.size()<=numBestMaps);

	if (!furtherPairsSearch)
	{
		potentialCorrFrom.clear();
		potentialCorrTo.clear();
	}
	else
	{
		SurfaceDistance * distMetric1 = GetSurface(0)->GetDistanceMetric("default");
		SurfaceDistance * distMetric2 = GetSurface(1)->GetDistanceMetric("default");	
		double delta1 = GetSurface(0)->AdjustedRadius(.1);
		double delta2 = GetSurface(1)->AdjustedRadius(.1);
		double ** featureFrom = new double*[m_pSet1->NumSamples()];
		double ** featureTo = new double*[m_pSet2->NumSamples()];
		
//		std::cout<<"Assigning Greedy Corrs prepare feature vecs 1"<<std::endl;
		// prepare all feature vectors
		for (int i=0; i<m_pSet1->NumSamples(); i++)
		{
			std::vector<double> medianSearch1;
			featureFrom[i] = new double[goodSetFrom.size()];
			const SurfaceSample & samp = m_pSet1->GetSample(i);		
			for (int j=0; j<(int)goodSetFrom.size(); j++)
			{
//				std::cout<<"Getting sample from: "<<goodSetFrom[i]<<std::endl;
				const SurfaceSample & goodSampFrom = m_pSet1->GetSample(goodSetFrom[j]);
				featureFrom[i][j] = distMetric1->Distance(samp, goodSampFrom);
				medianSearch1.push_back(featureFrom[i][j]);
			}
			double median1 = medianSearch1[medianSearch1.size()/2];	
			for (int j=0; j<(int)goodSetFrom.size(); j++)
				featureFrom[i][j] /= median1;		
		}
//		std::cout<<"Assigning Greedy Corrs prepare feature vecs 2"<<std::endl;		
		for (int i=0; i<m_pSet2->NumSamples(); i++)
		{
			std::vector<double> medianSearch2;
			featureTo[i] = new double[goodSetTo.size()];
			const SurfaceSample & samp = m_pSet2->GetSample(i);		
			for (int j=0; j<(int)goodSetTo.size(); j++)
			{
				const SurfaceSample & goodSampTo = m_pSet2->GetSample(goodSetTo[j]);
				featureTo[i][j] = distMetric2->Distance(samp, goodSampTo);
				medianSearch2.push_back(featureTo[i][j]);
			}
			double median2 = medianSearch2[medianSearch2.size()/2];	
			for (int j=0; j<(int)goodSetTo.size(); j++)
				featureTo[i][j] /= median2;	
		}
		
//		std::cout<<"Assigning Greedy Corrs - check feature vecs"<<std::endl;
		int numGoodCorrs = (int)goodSetTo.size();
		// go over all potential corrs, check feature vectors
		for (int i=0; i<(int)potentialCorrTo.size(); i++)
		{
			// for each potential correspondence find closest feature value
			int potSampFromID = potentialCorrFrom[i];
			int potSampToID = potentialCorrTo[i];		
			const SurfaceSample & potSampFrom = m_pSet1->GetSample(potSampFromID);
			const SurfaceSample & potSampTo = m_pSet2->GetSample(potSampToID);
			
//			std::cout<<"1. Corr for "<<potSampFromID;
			bool b1 = FPS_NearestFeatureWithinThrehsold(m_pSet2, featureFrom[potSampFromID],
														potSampTo, featureTo, distMetric2, 
														delta2, numGoodCorrs);
			if (!b1)
			{
				potentialCorrFrom[i] = -1;
				potentialCorrTo[i] = -1;				
			}
			else
			{
//				std::cout<<"2. Corr for "<<potSampToID;				
				bool b2 = FPS_NearestFeatureWithinThrehsold(m_pSet1, featureTo[potSampToID],
															potSampFrom, featureFrom, distMetric1,
															delta1, numGoodCorrs);
				if (!b2)
				{
					potentialCorrFrom[i] = -1;
					potentialCorrTo[i] = -1;					
				}
			}
			
			// check that sample with closest feature is within the threshold
		}
		
		for (int i=0; i<m_pSet1->NumSamples(); i++)
			 delete [] featureFrom[i];
		for (int i=0; i<m_pSet2->NumSamples(); i++)
			 delete [] featureTo[i];
		delete [] featureFrom;
		delete [] featureTo;
	}
//	std::cout<<"Assigning Greedy Corrs - done all - ready to return"<<std::endl;			
	// add potential + good corrs
	for (int i=0; i<(int)m_mFinalCorrMapM1ToM2.size(); i++)
		m_mFinalCorrMapM1ToM2[i] = -1;

	assert(goodSetFrom.size()==goodSetTo.size());
	assert(potentialCorrFrom.size()==potentialCorrTo.size());	
	for (int i=0; i<(int)goodSetFrom.size(); i++)
		if (goodSetFrom[i]!=-1 && goodSetTo[i]!=-1)
		{
			m_mFinalCorrMapM1ToM2[goodSetFrom[i]] = goodSetTo[i];
			AnalysisStats::m_GlobalStats.Increase("MobiusVotingGoodCorrs");
		}
	for (int i=0; i<(int)potentialCorrFrom.size(); i++)
	{
		if (potentialCorrFrom[i]!=-1 && potentialCorrTo[i]!=-1)
		{
			m_mFinalCorrMapM1ToM2[potentialCorrFrom[i]] = potentialCorrTo[i];
			AnalysisStats::m_GlobalStats.Increase("MobiusVotingFurtherPairs");
		}
		else
			AnalysisStats::m_GlobalStats.Increase("MobiusVotingRejectedFurtherPairs");
	}

//	std::cout<<"Assigning Greedy Corrs - returning"<<std::endl;				
}

bool MapCoarse::FPS_NearestFeatureWithinThrehsold(SurfaceSampleSet * searchInSet,
												  double * refFeature,
												  const SurfaceSample & refSample,
												  double ** featuresInSearch,
												  SurfaceDistance * distMetric,
												  double threshold, int numGoodCorrs)
{
	double bestFeatureDistance = 0;
	int bestFeatureCorrID = -1;
	for (int j=0; j<searchInSet->NumSamples(); j++)
	{
		double * aj = refFeature;
		double * bl = featuresInSearch[j];
		double featDist = FPS_FeatureDistance(aj, bl, numGoodCorrs);
		if (featDist < bestFeatureDistance || bestFeatureCorrID==-1)
		{
			bestFeatureDistance = featDist;
			bestFeatureCorrID = j;
		}
	}
//	std::cout<<" closest feature = "<<bestFeatureCorrID<<std::flush;
	assert(bestFeatureCorrID!=-1);
	double dist = distMetric->Distance(refSample, searchInSet->GetSample(bestFeatureCorrID));
//	std::cout<<" Accepted="<<(dist < threshold)<<" dist="<<dist<<" tau="<<threshold<<std::endl;
	return (dist < threshold);
}

double MapCoarse::FPS_FeatureDistance(double * a, double * b, int vecSize)
{
	double diff2 = 0;
	double normA2 = 0;
	double normB2 = 0;	
	for (int i=0; i<vecSize; i++)
	{
		diff2 += pow(a[i] - b[i], 2.);
		normA2 += pow(a[i], 2.);
		normB2 += pow(b[i], 2.);
	}
	
	return diff2 / (sqrt(normA2) * sqrt(normB2));
}

void MapCoarse::FillBestCorr(SurfaceSampleSet ** saveSet1, 
							 SurfaceSampleSet ** saveSet2, int numBest) const 
{
	vector<Sortable> sortedCorrs;
	for (int i=0; i<(int)m_mFinalCorrValue.size(); i++)
		sortedCorrs.push_back(Sortable(m_mFinalCorrValue[i], NULL, i));
	sort(sortedCorrs.begin(), sortedCorrs.end());
	
	*saveSet1 = new SurfaceSampleSet();
	*saveSet2 = new SurfaceSampleSet();

	if (numBest<0)
		numBest = sortedCorrs.size();

	for (int i=(int)sortedCorrs.size()-1; i>=0; i--)
	{
		int s1 = sortedCorrs[i].id;
		int s2 = m_mFinalCorrMapM1ToM2[s1];
		if (s2==-1)
			continue;
		assert(s1>=0 && s1<m_pSet1->NumSamples());
		assert(s2>=0 && s2<m_pSet2->NumSamples());		
		(*saveSet1)->AddSample(m_pSet1->GetSample(s1));
		(*saveSet2)->AddSample(m_pSet2->GetSample(s2));
		if ((*saveSet1)->NumSamples()>=numBest)
			break;
	}
}

void MapCoarse::CleanupMap(double takeBestMaps)
{
	// remove double pointers
//	int rangeSize = GetValidRange()->NumSamples();
//	int * backPntr = new int[rangeSize];
//	for (int i=0; i<rangeSize; i++)
//		backPntr[i] = -1;
//	
//	for (int s1=0; s1<(int)m_mFinalCorrMapM1ToM2.size(); s1++)
//	{
//		int s2 = m_mFinalCorrMapM1ToM2[s1];
//		if (s2==-1)
//			continue;
//		assert(s2<rangeSize);
//		
//		double currValue = m_mFinalCorrValue[s1];
//		int oldS1 = backPntr[s2];
//		if (oldS1<0)	// no back pointer
//		{
//			backPntr[s2] = s1;
//		}
//		else if (m_mFinalCorrValue[oldS1] < currValue)	//new value is higher
//		{
//			backPntr[s2] = s1;
//			m_mFinalCorrMapM1ToM2[oldS1] = -1;
//			m_mFinalCorrValue[oldS1] = -1;
//		}
//		else		//old error is smaller
//		{	
//			m_mFinalCorrMapM1ToM2[s1] = -1;		
//			m_mFinalCorrValue[s1]=-1;
//		}		
//	}
//	delete [] backPntr;

	// remove worst-value maps
	std::vector<Sortable> sortedCorrespondences;
	
	for (int s1=0; s1<(int)m_mFinalCorrMapM1ToM2.size(); s1++)
	{
		if (m_mFinalCorrMapM1ToM2[s1]!=-1)
			sortedCorrespondences.push_back(Sortable(m_mFinalCorrValue[s1], NULL, s1));
	}
	std::sort(sortedCorrespondences.begin(), sortedCorrespondences.end());
	
	int numBest = (int)(m_mFinalCorrMapM1ToM2.size()*takeBestMaps);	
	
	for (int i=0; i<((int)m_mFinalCorrMapM1ToM2.size()-numBest); i++)
	{
//		std::cout<<"Cleanup "<<this<<" : ["<<sortedCorrespondences[i].id<<"] -> ";
//		std::cout<<m_mFinalCorrMapM1ToM2[sortedCorrespondences[i].id];
//		std::cout<<" Value = "<<m_mFinalCorrValue[sortedCorrespondences[i].id]<<std::endl;
		m_mFinalCorrMapM1ToM2[sortedCorrespondences[i].id]=-1;
		m_mFinalCorrValue[sortedCorrespondences[i].id]=-1;
	}
}
	
SurfaceSample MapCoarse::ForwardMap(const SurfaceSample & s)
{
	if(m_mVertexIDToSetID1.find(s.NearestVertex())==m_mVertexIDToSetID1.end())
		return SurfaceSample();

	assert((int)m_mFinalCorrMapM1ToM2.size()>0);
	
	int sID = m_mVertexIDToSetID1[s.NearestVertex()];

	if (m_mFinalCorrMapM1ToM2[sID]>=0)
		return m_pSet2->GetSample(m_mFinalCorrMapM1ToM2[sID]);
	else 
		return SurfaceSample();
}

SurfaceSample MapCoarse::InverseMap(const SurfaceSample & s)
{
	assert((int)m_mFinalCorrMapM1ToM2.size()>0);
	
	int sID = m_mVertexIDToSetID2[s.NearestVertex()];

	for (int i=0; i<(int)m_mFinalCorrMapM1ToM2.size(); i++)
		if (m_mFinalCorrMapM1ToM2[i]==sID)
			return m_pSet1->GetSample(i);

	return SurfaceSample();	
}

int MapCoarse::GetMappedDomainSize()
{
	int validSize=0; 
	for (int i=0; i<m_pSet1->NumSamples(); i++)
		if (!ForwardMap(m_pSet1->GetSample(i)).Invalid())
			validSize++;
	return validSize;
}

const SurfaceSampleSet * MapCoarse::GetValidDomain() const
{
	return m_pSet1;
}

const SurfaceSampleSet * MapCoarse::GetValidRange() const
{
	return m_pSet2;
}

bool MapCoarse::LoadMapInOneWayVertexIDsFormat(const VKString & filename, bool surf0ToSurf1)
{
	std::ifstream textStream(filename.c_str());
	if (!textStream.is_open())
		return false;
	
	assert(surf0ToSurf1);	// not implemented in the other direction
	
	SurfaceSampleSet * set1 = new SurfaceSampleSet(GetSurface(0)->GetMesh());
	SurfaceSampleSet * set2 = new SurfaceSampleSet(GetSurface(1)->GetMesh());	
	std::vector<int> finalMap;
	for (int i=0; i<GetSurface(0)->GetMesh()->NVertices(); i++)
	{
		int toID;			
		textStream>>toID;		
		finalMap.push_back(toID);
	}
	
	InitializeSets(set1, set2);
	SetFinalCorrMap(finalMap);
	
	return true;
}

bool MapCoarse::LoadMapInVertexByVertexFormat(const VKString & filename)
{
	std::ifstream textStream(filename.c_str());
	if (!textStream.is_open())
		return false;
	
	SurfaceSampleSet * set1 = new SurfaceSampleSet();
	SurfaceSampleSet * set2 = new SurfaceSampleSet();	
	std::vector<int> finalMap;
	
	while(!textStream.eof())
	{
		int v1, v2;
		textStream>>v1>>v2;
		set1->AddSample(SurfaceSample(v1, GetSurface(0)->GetMesh()));
		set2->AddSample(SurfaceSample(v2, GetSurface(1)->GetMesh()));		
		finalMap.push_back((int)finalMap.size());
	}
	InitializeSets(set1, set2);
	SetFinalCorrMap(finalMap);
	return true;
}

bool MapCoarse::LoadMapDensePreceiseFormat(const VKString & filename)
{
	SurfaceSampleSet * allVertices = GetSurface(0)->GetSampleSet("AllVertices");
	SurfaceSampleSet * mesh2Set = new SurfaceSampleSet();
	
	std::ifstream textStream(filename.c_str());
	if (!textStream.is_open())
		return false;
	
	std::vector<int> finalMap;
	int triID;
	double b1, b2, b3;
	for (int i=0; i<allVertices->NumSamples(); i++)
	{
		textStream>>triID>>b1>>b2>>b3;
		assert(textStream.is_open());
		assert(b1+b2+b3==1. && b1>=0 && b2>=0 && b3>=0);
		assert(triID < GetSurface(1)->GetMesh()->NFaces());
		mesh2Set->AddSample(SurfaceSample(triID, b1, b2, b3, GetSurface(1)->GetMesh()));
		finalMap.push_back((int)finalMap.size());
	}
	
	InitializeSets(allVertices, mesh2Set);
	SetFinalCorrMap(finalMap);
	return true;
}

void MapCoarse::SaveMapInVertexByVertexFormat(const VKString & filename)
{
	std::ofstream textStream(filename.c_str());
	assert(textStream.is_open());
	for (int i=0; i<(int)m_mFinalCorrMapM1ToM2.size(); i++)
	{
		if (m_mFinalCorrMapM1ToM2[i]>=0)
		{
			textStream<<m_pSet1->GetSample(i).NearestVertex()<<"\n";
			textStream<<m_pSet2->GetSample(m_mFinalCorrMapM1ToM2[i]).NearestVertex()<<"\n";			
		}
	}	
}

void MapCoarse::SaveMap(std::ofstream & textStream)
{
	WriteCorrByRow(textStream);
}

void MapCoarse::LoadMap(std::ifstream & textStream)
{
	LoadCorrByRow(textStream);
}

void MapCoarse::WriteCorrByRow(std::ofstream & textStream)
{
	std::vector<int> goodMapV1;
	std::vector<int> goodMapV2;
	std::vector<double> mapError;
	for (int i=0; i<(int)m_mFinalCorrMapM1ToM2.size(); i++)
	{
		if (m_mFinalCorrMapM1ToM2[i]>=0)
		{
			goodMapV1.push_back(m_pSet1->GetSample(i).NearestVertex());
			goodMapV2.push_back(m_pSet2->GetSample(m_mFinalCorrMapM1ToM2[i]).NearestVertex());
			mapError.push_back(m_mFinalCorrValue[i]);
		}
	}

	textStream<<"Format MapCoarse\n";
	textStream<<"CoarseType Pairwise\n";	
	textStream<<"NumMeshes 2\n";
	textStream<<"SampleTypes NearestVertexID_Error\n";
	textStream<<"NVertices1 "<<GetSurface(0)->GetMesh()->NVertices()<<"\n";
	textStream<<"NVertices2 "<<GetSurface(1)->GetMesh()->NVertices()<<"\n";	
	
	// num correspondences per entry
	textStream<<2<<" ";
	// mesh ids
	textStream<<0<<" "<<1<<" ";
	// num entries in this block 
	textStream<<goodMapV1.size()<<"\n";
	
	for (int i=0; i<(int)goodMapV1.size(); i++)
	{
		textStream<<goodMapV1[i]<<" "<<mapError[i]<<"\n";
		textStream<<goodMapV2[i]<<" "<<mapError[i]<<"\n";
	}
}

bool MapCoarse::LoadCorrByRow(std::ifstream & textStream)
{
	std::vector<int> vertexOnM1;
	std::vector<int> vertexOnM2;
	std::vector<double> errors;
	
//	std::ifstream textStream(filename.c_str(), ios::in);
//	if(!textStream.is_open())
//		return false;
	std::string emptyStr;
	int intVal;
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="Format");
	textStream>>emptyStr;	
	assert(VKString(emptyStr.c_str())=="MapCoarse");
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="CoarseType");
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="Pairwise");
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="NumMeshes");
	textStream>>intVal;		assert(intVal==2);
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="SampleTypes");
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="NearestVertexID_Error");
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="NVertices1");
	textStream>>intVal; assert(GetSurface(0)->GetMesh()->NVertices()==intVal);
	textStream>>emptyStr;	assert(VKString(emptyStr.c_str())=="NVertices2");
	textStream>>intVal; assert(GetSurface(1)->GetMesh()->NVertices()==intVal);	
	
	int numCorrPerEntry;
	textStream>>numCorrPerEntry;
	assert(numCorrPerEntry==2);		// bilateral symmetry
	
	int meshID1, meshID2;
	textStream>>meshID1>>meshID2;	// meshID
	assert(meshID1==0 && meshID2==1);
	int numEntriesInBlock;
	textStream>>numEntriesInBlock;
	
	std::vector<int> finalMap;
	for (int i=0; i<numEntriesInBlock; i++)
	{
		int v1ID, v2ID;
		double err;
		textStream>>v1ID>>err;
		textStream>>v2ID>>err;
		vertexOnM1.push_back(v1ID);
		vertexOnM2.push_back(v2ID);		
		errors.push_back(err);
		finalMap.push_back(i);
	}
	
	SurfaceSampleSet * set1 = new SurfaceSampleSet(GetSurface(0)->GetMesh(), vertexOnM1);
	SurfaceSampleSet * set2 = new SurfaceSampleSet(GetSurface(1)->GetMesh(), vertexOnM2);	
	
	InitializeSets(set1, set2);
	SetFinalCorrMap(finalMap);
	assert(errors.size()==m_mFinalCorrValue.size());
	for (int i=0; i<(int)errors.size(); i++)
		m_mFinalCorrValue[i] = errors[i];
	return true;
}

bool MapCoarse::FillFeatureToFeatureMap(const VKString & featurePntsMap, 
										std::map<int, int> & featureToFeatureMap)
{
	std::ifstream textStream(featurePntsMap.c_str(), ios::in);
	if(!textStream.is_open())
	{
		std::cout<<"[WARNING] Could not load "<<featurePntsMap.c_str()<<std::endl;
		return false;
	}
	VKString tempStr;
	while(!textStream.eof())
	{
		int feat1;
		int feat2;
		tempStr.readLine(textStream);
		
		VKStringList tempStrList = tempStr.split(" ");
		if (tempStrList.count()==0)
			continue;
		else if (tempStrList.count()==1)
		{
			feat1=tempStrList[0].toInt();
			feat2=feat1;
		}
		else if (tempStrList.count()==2)
		{
			feat1=tempStrList[0].toInt();
			feat2=tempStrList[1].toInt();			
		}
		else
			assert(false);
		featureToFeatureMap[feat1]=feat2;
		featureToFeatureMap[feat2]=feat1;
	}	
	return true;
}


bool MapCoarse::LoadCorrFeaturePntSymmetryMap(const VKString & verticesToFeatPnts, 
											  const VKString & featurePntsMap)
{
	std::map<int, int> featurePointToFeaturePoint;
	std::vector<int> featurePointToVertex;

	FillFeatureToVertex(verticesToFeatPnts, featurePointToVertex);
	
	if (!FillFeatureToFeatureMap(featurePntsMap, featurePointToFeaturePoint))
		return false;
	
	std::vector<int> finalMap;
	for (int i=0; i<(int)featurePointToVertex.size(); i++)
	{
		assert(featurePointToFeaturePoint.find(i)!=featurePointToFeaturePoint.end());
		finalMap.push_back(featurePointToFeaturePoint[i]);
	}
	
	SurfaceSampleSet * set1 = new SurfaceSampleSet(GetSurface(0)->GetMesh(), featurePointToVertex);
	SurfaceSampleSet * set2 = new SurfaceSampleSet(GetSurface(1)->GetMesh(), featurePointToVertex);
	InitializeSets(set1, set2);
	SetFinalCorrMap(finalMap);
	assert(m_mFinalCorrMapM1ToM2.size()>0);
	return true;
}

bool MapCoarse::LoadCorrSameFeaturePntIDs(const VKString & verticesToFeatPnts1, 
										  const VKString & verticesToFeatPnts2)
{
	std::vector<int> featurePointToVertex1;
	std::vector<int> featurePointToVertex2;	
	
	FillFeatureToVertex(verticesToFeatPnts1, featurePointToVertex1);
	FillFeatureToVertex(verticesToFeatPnts2, featurePointToVertex2);	
	
	SurfaceSampleSet * set1 = new SurfaceSampleSet(GetSurface(0)->GetMesh(), featurePointToVertex1);
	SurfaceSampleSet * set2 = new SurfaceSampleSet(GetSurface(1)->GetMesh(), featurePointToVertex2);
	InitializeSets(set1, set2);
	assert(set1->NumSamples()==set2->NumSamples());

	std::vector<int> finalMap;
	for (int i=0; i<(int)featurePointToVertex1.size(); i++)
		finalMap.push_back(i);	// 1 - 1 map
	SetFinalCorrMap(finalMap);
	return true;
}

bool MapCoarse::LoadCorrSameFeaturePntIDsFlipSymmetry(const VKString & verticesToFeatPnts1, 
													  const VKString & verticesToFeatPnts2,
													  const VKString & symInfoFile)
{
	std::vector<int> featurePointToVertex1;
	std::vector<int> featurePointToVertex2;	
	
	FillFeatureToVertex(verticesToFeatPnts1, featurePointToVertex1);
	FillFeatureToVertex(verticesToFeatPnts2, featurePointToVertex2);	
	
	SurfaceSampleSet * set1 = new SurfaceSampleSet(GetSurface(0)->GetMesh(), featurePointToVertex1);
	SurfaceSampleSet * set2 = new SurfaceSampleSet(GetSurface(1)->GetMesh(), featurePointToVertex2);
	InitializeSets(set1, set2);
	assert(set1->NumSamples()==set2->NumSamples());

	std::map<int, int> featureToFeatureMap;
	if (!FillFeatureToFeatureMap(symInfoFile, featureToFeatureMap))
		return false;
	
	std::vector<int> finalMap;
	
	for (int i=0; i<(int)featurePointToVertex1.size(); i++)
		finalMap.push_back(featureToFeatureMap[i]);	// 1 - 1 map
	SetFinalCorrMap(finalMap);
	return true;
	
}

bool MapCoarse::LoadCorrSameVertexIDs()
{
	SurfaceSampleSet * set1 = new SurfaceSampleSet(GetSurface(0)->GetMesh());
	SurfaceSampleSet * set2 = new SurfaceSampleSet(GetSurface(1)->GetMesh());
	
	if (set1->NumSamples()!=set2->NumSamples())
	{
		std::cout<<"[ERROR] Different number of ground corrs: ";
		std::cout<<set1->NumSamples()<<" vs "<<set2->NumSamples()<<std::endl;
		assert(false);
	}
	InitializeSets(set1, set2);
	std::vector<int> finalMap;
	for (int i=0; i<(int)set1->NumSamples(); i++)
		finalMap.push_back(i);	// 1 - 1 map	
	SetFinalCorrMap(finalMap);
	return true;
}


void MapCoarse::FillFeatureToVertex(const VKString & vertexIDsOfFeatPnts, 
									std::vector<int> &vertexIDs)
{
	std::ifstream textStream1(vertexIDsOfFeatPnts.c_str(), ios::in);
	
	if (!textStream1.is_open())
	{
		std::cout<<"[ERROR] Cannot open file: "<<vertexIDsOfFeatPnts.c_str()<<std::endl;
		assert(false);
	}
	
	while(!textStream1.eof())
	{
		int vertexID;
		double x, y, z;
		textStream1>>vertexID;
		
		textStream1>>x>>y>>z;
		if (!textStream1.fail())
			vertexIDs.push_back(vertexID);
	}
}

VKString MapCoarse::GetSurfaceMapType()
{
	return "MapCoarse";
}

VKString MapCoarse::GetTruthConfidenceName()
{
	return m_sTruthConfidence;
}

void MapCoarse::SetCompareToTruthRendering(const VKString & confidenceName,
										   MapCoarse * trueMap)
{
	m_pTruthMapForDrawing = trueMap;
	m_sTruthConfidence = confidenceName;
	
//	SurfaceMapConfidence * confidence = GetMapConfidenceCalculator(confidenceName);
//	double correspondenceRate, sumOfGeoDistances;
//
//	confidence->GetPerSampleErrors(trueMap->GetSurface(0), errors, false);
//	correspondenceRate = confidence->Confidence();
//	sumOfGeoDistances = confidence->Error();
	
//	assert(m_mFinalCorrMapM1ToM2.size()==errors.size() && m_mFinalCorrError.size()==errors.size());
//	m_fMaxErrorToTruth=0;
//	for (int i=0; i<(int)errors.size(); i++)
//		m_fMaxErrorToTruth = vkMax(errors[i], m_fMaxErrorToTruth);
//	
//	std::map<int, double> vertexToErrorMap;
//	for (int i=0; i<(int)errors.size(); i++)
//	{
//		m_mFinalCorrError[i] = errors[i];
//		bool atVertex;
//		vertexToErrorMap[m_pSet1->GetSample(i).NearestVertex(&atVertex)] = errors[i];
//		assert(atVertex);
//		//m_pTruthMapForDrawing->m_mFinalCorrError[i] = errors[i];
//	}
//	
//	const SurfaceSampleSet * truthDomain = GetValidDomain();
//	for (int i=0; i<truthDomain->NumSamples(); i++)
//	{
//		bool atVertex;
//		int vertexID = truthDomain->GetSample(i).NearestVertex(&atVertex);
//		assert(atVertex);
//		m_pTruthMapForDrawing->m_mFinalCorrError[i] = vertexToErrorMap[vertexID];
//	}
//	
//	
//	m_pTruthMapForDrawing->m_fMaxErrorToTruth = m_fMaxErrorToTruth;
}

void MapCoarse::SetCorrespondenceColor(double r, double g, double b)
{
	m_CorrColorR = r;
	m_CorrColorG = g;
	m_CorrColorB = b;	
}

void MapCoarse::InitializeOrderOfDomainRendering()
{
	for (int i=0; i<m_pSet1->NumSamples(); i++)
		 m_OrderOfDomainRendering.push_back(i);
}

void MapCoarse::Draw(AnalysisWindow * window,
					 ParamParser * params, 
					 const VKString & renderingParams, 
					 const VKString & surfaceName)
{
	if (m_pRenderAnotherCoarse!=NULL)
		m_pRenderAnotherCoarse->Draw(window, params, renderingParams, surfaceName);
	SurfaceMapConfidence * confidence=NULL;
	if (m_sTruthConfidence!="")
		confidence = GetMapConfidenceCalculator(m_sTruthConfidence);
		
	DrawCorrespondenceForSelected(params, renderingParams, surfaceName);
	
	static GLfloat material[4];	
	material[0] = m_CorrColorR;	material[1] = m_CorrColorG;	material[2] = m_CorrColorB;	material[3] = 1;	
	
	bool valid;
	bool drawEdges = params->GetStrValues(renderingParams, "MapFlags", valid).contains("Edges");
	bool drawColorBalls = params->GetStrValues(renderingParams, "MapFlags", valid).contains("Colors");
	bool benchOnlyError = params->GetStrValues(renderingParams, "MapFlags", valid).contains("BenchErrorOnly");
		
	int numColors = params->GetIntValue(renderingParams, "MaxNumMapColors", valid);
	if (!valid || numColors==-1)
		numColors = m_mFinalCorrMapM1ToM2.size()+1;
	int maxNumEdges = params->GetIntValue(renderingParams, "MaxNumMapEdges", valid);
	int maxNumElements = params->GetIntValue(renderingParams, "MaxCoarseElements", valid);
	if ((int)m_mFinalCorrMapM1ToM2.size() > maxNumElements)
	{
		drawEdges = false;
		drawColorBalls = false;
	}
	if ((int)m_mFinalCorrMapM1ToM2.size() > maxNumEdges)
		drawEdges = false;
	
	double radius1 = GetSurface(0)->GetStandardRadius(params, renderingParams, surfaceName);
	double radius2 = GetSurface(1)->GetStandardRadius(params, renderingParams, surfaceName);
	
	if (m_OrderOfDomainRendering.size()==0)
		InitializeOrderOfDomainRendering();
	
	int validCorrs = 0;
	for (int o=0; o<(int)m_mFinalCorrMapM1ToM2.size()+1; o++)
	{
		int i=-1;
		if (o<(int)m_mFinalCorrMapM1ToM2.size())
			i=(int)m_OrderOfDomainRendering[o];
		else
			break;
		
		if (m_mFinalCorrMapM1ToM2[i]==-1)
			continue;
					
		if (confidence!=NULL && benchOnlyError)
		{
			SurfaceSample sampOn1 = m_pSet1->GetSample(i);
			assert(m_pTruthMapForDrawing!=NULL);
			SurfaceSample sampOn2 = m_pTruthMapForDrawing->ForwardMap(sampOn1);
			if (confidence->ConfidenceAtSample(sampOn2)==1.)
				continue;
		}
		else if (validCorrs>=numColors)
			break;
		
		validCorrs++;
		
		R3Point p1 = GetSurface(0)->GetDrawablePosition(m_pSet1->GetSample(i), 
														params, renderingParams, surfaceName);
		R3Point p2 = GetSurface(1)->GetDrawablePosition(m_pSet2->GetSample(m_mFinalCorrMapM1ToM2[i]), 
														params, renderingParams, surfaceName);

		double truthErr = 1.;

		if (p1==p2)	// check somehow else that it's one mesh
		{
			glEnable(GL_LIGHTING);
			double scaleBy=1.5;
			if (confidence!=NULL)
			{
				material[0] = 1;
				material[1] = 1-truthErr;
				material[2] = 1-truthErr;
				scaleBy = 2;
			}
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material);
			R3Sphere(p1, radius1*scaleBy).Draw();
		}
		else
		{
			if (drawEdges)
			{
				glLineWidth(1.);
				glDisable(GL_LIGHTING);	
				glBegin(GL_LINES);				
				//TODO: can make color vary depending on error value
				//glColor3d(m_CorrColorR, m_CorrColorG, m_CorrColorB);
				glColor3d(.8, 0, 0);
				glVertex3d(p1.X(), p1.Y(), p1.Z());
				glVertex3d(p2.X(), p2.Y(), p2.Z());
				glEnd();
				
				if (confidence!=NULL)
				{
					glEnable(GL_LIGHTING);
					material[0] = 1;
					material[1] = 1-truthErr;
					material[2] = 1-truthErr;
					glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material);
					R3Sphere(p1, radius1).Draw();					
				}
			}
			
			if (drawColorBalls && confidence==NULL)
			{
				double radScale = 1. + (100.-(double)numColors)/70.;
				if (radScale<1)
					radScale = 1;
				if (radScale>2)
					radScale = 2;
				if (numColors<10)
					radScale = 3;
				else
					radScale = 2;
				glEnable(GL_LIGHTING);
				window->MapIntToColor(validCorrs, material);
				glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material);				
				GetSurface(0)->LoadCamera();
				R3Point p1 = GetSurface(0)->GetDrawablePosition(m_pSet1->GetSample(i), 
																params, renderingParams, surfaceName);				
				R3Sphere(p1, radius1*radScale).Draw();
				GetSurface(1)->LoadCamera();				
				R3Point p2 = GetSurface(1)->GetDrawablePosition(m_pSet2->GetSample(m_mFinalCorrMapM1ToM2[i]), 
																params, renderingParams, surfaceName);				
				R3Sphere(p2, radius2*radScale).Draw();
			}
		}
	}
	
	if (m_pTruthMapForDrawing!=NULL && drawEdges && confidence!=NULL)
	{
//		m_pTruthMapForDrawing->m_fTruthThreshold=m_fTruthThreshold;
//		m_pTruthMapForDrawing->SetCorrespondenceColor(0, 0, 1);
//		m_pTruthMapForDrawing->Draw(window, params, renderingParams, surfaceName);
//		m_pTruthMapForDrawing->m_fTruthThreshold=-1;		
	}
}


MapCoarse::FineToCoarseGeneration MapCoarse::StrToCoarseningMethod(const VKString & method)
{
	if (method=="F2C_MUTUALLY_CLOSEST_NEIGHBORS")
		return F2C_MUTUALLY_CLOSEST_NEIGHBORS;
	else if (method=="F2C_LINEAR_ASSIGNMENT")
		return F2C_LINEAR_ASSIGNMENT;			
	else if (method=="F2C_FORWARD_ONLY")
		return F2C_FORWARD_ONLY;
	else if (method=="F2C_BACKWARD_ONLY")
		return F2C_BACKWARD_ONLY;
	else if (method=="F2C_FORWARD_BACKWARD")
		return F2C_FORWARD_BACKWARD;
	else if (method=="F2C_MCN_CONFORMAL_EUCLIDEAN")
		return F2C_MCN_CONFORMAL_EUCLIDEAN;
	else
	{
		std::cout<<"[ERROR] AnalysisPipeline.cpp: Conformal to coarse method is not defined!"<<std::endl;
		assert(false);
	}
}

