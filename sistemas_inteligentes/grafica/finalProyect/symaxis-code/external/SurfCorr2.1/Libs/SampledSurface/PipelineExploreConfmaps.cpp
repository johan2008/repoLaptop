#include "PipelineExploreConfmaps.h"
#include "MapMultiConformal.h"
#include <algorithm>
#include "Sortable.h"
#include "LinAlgMatrixSparseReal.h"

R3Point global_SHREC10CorrVertex(SHREC10CorrVertex * vertex, void * )
{
	return vertex->pnt;
}

PipelineExploreConfmaps::PipelineExploreConfmaps(ParamParser & params) 
: AnalysisPipeline(params)
{	
	m_VertexToFuzzySetID1 = NULL;
	m_VertexToFuzzyVertexID2 = NULL;
	m_aFuzzyValues = NULL;

	
	m_sCoarseSetName = "none";
	m_pCollectionOfConformalMaps = NULL;
	m_CurrentSequence = NULL;
	m_pFinalCoarseMap = NULL;
	
	m_pIntegratedBlendingMatrix = NULL;
	m_pEigenvectorsIBM = NULL;
	m_pEigenvaluesIBM = NULL;
	m_pCachedConfidenceVals = NULL;
	
	VKString algorithmName = params.GetStrValue("Pipeline", "Algorithm", valid);
	assert(valid);
	
	SurfaceCorrespondenceConfMap(algorithmName);
}

PipelineExploreConfmaps::~PipelineExploreConfmaps()
{
}

void PipelineExploreConfmaps::SurfaceCorrespondenceConfMap(const VKString & algorithmName)
{
	WriteLog();	
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("SurfaceCorrespondenceConfMap");
	
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);	
	
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("AABasicCreateSurfaces");
	LoadAndProcessSurfaces(algorithmName);		// Load / create Features, Samples, Distances
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("AABasicCreateSurfaces");	
	
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("AABasicCreateMaps");		
	if (m_Params.GetStrValue(algorithmName, "Type", valid)=="OneStepSpectral")
	{
		LoadAndProcessOneStepSpectralMap(algorithmName);
	}
	else
	{
		LoadSetSamplers(algorithmName);				// Load / create generators for conformal maps
		LoadAndProcessMaps(algorithmName);			// Load / create maps
	}
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("AABasicCreateMaps");	
	LoadAndProcessBenchmark();					// Load / test on benchmark
	WriteLog();	
	
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("AABasicAdditionalOutput");
//	if (m_Params.GetStrValue("Pipeline", "LoadDensePreceise", valid)=="true" && valid)
//	{
//	}
	
	WritePostFinalProcessingResults();		// write various temp processing results
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("AABasicAdditionalOutput");		

	AnalysisStats::m_GlobalStats.PrintStats();	
}

void PipelineExploreConfmaps::LoadAndProcessMaps(const VKString & algorithmName)
{
	// Load and process maps
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);		
	if (loads.contains("Maps") || process.contains("Maps"))
	{
		bool needsCoarse, needsCoarse2Fine;
		LoadExternalMap(algorithmName, needsCoarse, needsCoarse2Fine);	// check if map is given
		if (needsCoarse || needsCoarse2Fine)
		{
			// load confidences and similarities
			VKStringList confidencesNames = m_Params.GetStrValues(algorithmName, "LoadConfidences", valid);	assert(valid);
			VKStringList similaritiesNames = m_Params.GetStrValues(algorithmName, "LoadSimilarities", valid);
			for (int i=0; i<confidencesNames.count(); i++)
				AddClonableConfidence(confidencesNames[i]);
			for (int i=0; i<similaritiesNames.count(); i++)
				AddClonableSimilarity(similaritiesNames[i]);	
			
			if (needsCoarse)	// still need to find coarse map - nothing was loaded
				FindCoarseCorrespondences(algorithmName);
			WriteLog();
			if (needsCoarse2Fine)	// still need to do coarse-to-fine - only coarse or nothing was loaded, 
				InterpolateCoarseCorrespondences(algorithmName);			
			WriteLog();
		}
	}
}

void PipelineExploreConfmaps::LoadExternalMap(const VKString & algorithmName,
											  bool & stillNeedCoarse, bool & stillNeedCoarseToFine)
{
	// LOAD EXTERNAL GROUND TRUTH FOR BENCHMARKING
	VKString loadExternal = m_Params.GetStrValue("Pipeline", "LoadExternalResult", valid);
	if (!valid)
		loadExternal = "none";
	stillNeedCoarse = true;
	stillNeedCoarseToFine = true;
		// loading from external file
	if (loadExternal!="none")
	{
		stillNeedCoarse = false;
		// load confidences and similarities
		AddClonableConfidence(m_Params.GetStrValue("Pipeline", "BenchmarkQuery", valid));
		AddClonableConfidence(m_Params.GetStrValue("Pipeline", "SecondBenchmark", valid));
		VKString externalMapFilename = m_Params.GetStrValue("Pipeline", "FinalExternalMap", valid);
		assert(valid);
		std::cout<<"Loading External Result: "<<loadExternal.c_str()<<std::endl;;
		// try to load time
		VKString timeFile = externalMapFilename;
		timeFile.replace(".map", ".time");
		std::ifstream tempTextStream(timeFile.c_str());
		if (tempTextStream.is_open())
		{
			tempTextStream.close();
			m_pExternalAlgorithmTime = new ParamParser(timeFile);
		}
		else
		{
			m_pExternalAlgorithmTime = new ParamParser();
			VKString zeroStr = VKString::number(0);
			m_pExternalAlgorithmTime->m_ParamModules["ExternalTiming"]["InitTime"].push_back(zeroStr);
			m_pExternalAlgorithmTime->m_ParamModules["ExternalTiming"]["FindMaps"].push_back(zeroStr);			
		}
		
		if (loadExternal=="FunkMap")
		{
			MapCoarse * theFinalMap = new MapCoarse(GetSurface(0), GetSurface(1));
			bool loaded = theFinalMap->LoadMapInOneWayVertexIDsFormat(externalMapFilename, true);
			assert(loaded);
			m_pTheFinalMap = theFinalMap;
			stillNeedCoarseToFine = false;
			
			m_pFinalCoarseMap = new MapCoarse(GetSurface(0), GetSurface(1));
			VKString coarseName = externalMapFilename;
			coarseName.replace(".map", ".cor");
			loaded = m_pFinalCoarseMap->LoadMapInVertexByVertexFormat(coarseName);
			if (!loaded)
				m_pFinalCoarseMap = NULL;
			else
				theFinalMap->AssociateCoarseToRender(m_pFinalCoarseMap);
		}
		else if (loadExternal=="FunkCoarse")
		{
			MapCoarse * finalCoarse = new MapCoarse(GetSurface(0), GetSurface(1));
			bool loaded = finalCoarse->LoadMapInVertexByVertexFormat(externalMapFilename);
			if (!loaded)
				std::cout<<"[ERROR] Could not load "<<externalMapFilename.c_str()<<std::endl;
			assert(loaded);
			m_pFinalCoarseMap = finalCoarse;
		}
		else
		{
			std::cout<<"[ERROR] Unknown External Map Type "<<loadExternal.c_str()<<std::endl;
			assert(false);
		}
	}
	
		// loading from truth
	if (m_Params.GetStrValue(algorithmName, "AggregationMethod", valid)=="LoadTruth" && valid)
	{
		LoadCoarseCorrsFromTruth();
		stillNeedCoarse = false;
	}
}

void PipelineExploreConfmaps::FindCoarseCorrespondences(const VKString & algorithmName)
{
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);	
	VKStringList processMaps = m_Params.GetStrValues("Pipeline", "MapProcess", valid);
	
	VKString coarseFinalMapExt = m_Params.GetStrValue(algorithmName, "CoarseFinalMapExt", valid);		
	VKString confCollExt = m_Params.GetStrValue(algorithmName, "ConfCollectExt", valid);
	VKString aggrCollExt = m_Params.GetStrValue(algorithmName, "MultiCollectExt", valid);
	VKString clusterExt = ".clustering";
		
	assert(loads.contains("Surface") && loads.contains("Samples") && loads.contains("Features"));
	
	if (Verb("Min"))
		std::cout<<"Creating Coarse Correspondences"<<std::endl;
	
	bool collectionAggregated=false;
	
	// LOADING
	if (!process.contains("Maps") || !processMaps.contains("ExploreConformals"))	// try loading conformal collection
		m_pCollectionOfConformalMaps = (MapScoredCollection*)SurfaceMap::CreateMap(LocalFilename(confCollExt), GetSurface(0), GetSurface(1));
	
	if (!process.contains("Maps") || !processMaps.contains("Aggregate"))
	{		// Try loading aggregated maps
		collectionAggregated = true;
		TryLoadCoarseMap(LocalFilename(coarseFinalMapExt), &m_pFinalCoarseMap);
		m_pCollectionOfMaps = (MapScoredCollection*)SurfaceMap::CreateMap(LocalFilename(aggrCollExt), GetSurface(0), GetSurface(1));
		if (m_pCollectionOfMaps==NULL)	// could not find collection of maps
		{								//		- set collection of maps to conformal maps
			collectionAggregated = false;
			m_pCollectionOfMaps = m_pCollectionOfConformalMaps;
		}
	}
	
	// PROCESSING (if necessary)
	double fracKeep = m_Params.GetDoubleValue(algorithmName, "PickBestConf", valid);
	VKString aggregateConfData = m_Params.GetStrValue(algorithmName, "AggregationMethod", valid);	

	bool exploredConformalMaps = false;
	if (m_pCollectionOfConformalMaps == NULL)
	{		// Explore Conformal maps
		AnalysisStats::m_GlobalStats.m_Timing.startedProcess("AABasicCreateMapsSearch");
		exploredConformalMaps = true;
		// initialize voting / correspondence extractor (scored collection)				
		m_pCollectionOfConformalMaps = new MapScoredCollection(GetSurface(0), GetSurface(1));
		ExploreConformalMaps();
		m_pCollectionOfConformalMaps->SortMaps();
		((SurfaceMap*)m_pCollectionOfConformalMaps)->SaveMap(LocalFilename(confCollExt));
		collectionAggregated = false;
		if (m_pCollectionOfMaps==NULL)
			m_pCollectionOfMaps = m_pCollectionOfConformalMaps;
		AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("AABasicCreateMapsSearch");		
	}
	
	// Aggregate conformal maps	
	if (aggregateConfData=="PickBest")		// pick best map
	{
		assert(m_pCollectionOfConformalMaps!=NULL);
		if (m_pFinalCoarseMap==NULL)
		{
			assert(m_pCollectionOfConformalMaps->GetNumMaps()>0);
			m_pFinalCoarseMap = m_pCollectionOfConformalMaps->GetBestCoarseMap();
			if (m_pFinalCoarseMap==NULL)
			{
				VKString sampleSetName = m_Params.GetStrValue(algorithmName, "SamplesFineCorr", valid);	assert(valid);
				SurfaceSampleSet * sampleSet1 = GetSurface(0)->GetSampleSet(sampleSetName);	assert(sampleSet1!=NULL);
				SurfaceSampleSet * sampleSet2 = GetSurface(1)->GetSampleSet(sampleSetName);	assert(sampleSet2!=NULL);				
				// get best map - make it into correspondence
				SurfaceMap * bestMap = m_pCollectionOfConformalMaps->GetBestMap();
				m_pFinalCoarseMap = new MapCoarse(sampleSet1, sampleSet2, bestMap, 
												  MapCoarse::F2C_MUTUALLY_CLOSEST_NEIGHBORS);
			}
			assert(m_pFinalCoarseMap!=NULL);
			if (exploredConformalMaps)
				m_pFinalCoarseMap->CleanupMap(fracKeep);
			((SurfaceMap*)m_pFinalCoarseMap)->SaveMap(LocalFilename(coarseFinalMapExt));
			m_pCollectionOfMaps = m_pCollectionOfConformalMaps;
		}
	}
	else if (aggregateConfData=="VoteInCorrMatrixGreedy")	// Mobius Voting
	{
		assert(m_pFinalCoarseMap!=NULL);
		((SurfaceMap*)m_pFinalCoarseMap)->SaveMap(LocalFilename(coarseFinalMapExt));				
		m_pCollectionOfMaps = m_pCollectionOfConformalMaps;		
	}
	else if (aggregateConfData!="LoadTruth")
	{
		std::cout<<"[ERROR] Unknown aggregation method: "<<aggregateConfData.c_str()<<std::endl;
		assert(false);
	}
}

void PipelineExploreConfmaps::InterpolateCoarseCorrespondences(const VKString & algorithmName)
{
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);	
	VKStringList processMaps = m_Params.GetStrValues("Pipeline", "MapProcess", valid);
	
	VKString finalMapExt = m_Params.GetStrValue(algorithmName, "FinalMapExt", valid);	
	
	// Try loading extrapolated maps	
	if (!process.contains("Maps") || !processMaps.contains("Extrapolate"))
		m_pTheFinalMap = SurfaceMap::CreateMap(LocalFilename(finalMapExt), GetSurface(0), GetSurface(1));
	
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("AABasicCreateMapsWeights");		
	VKString extrapolationType = m_Params.GetStrValue(algorithmName, "ExtrapolationMethod", valid);			
	// Extrapolate aggregated data
	if (m_pTheFinalMap==NULL)
	{
		std::cout<<"\tExtrapolating Aggregated Data: "<<extrapolationType.c_str()<<std::endl;
		assert(m_pFinalCoarseMap!=NULL);
		if (extrapolationType=="none")	// generally for the best conformal map
		{
			assert(m_pCollectionOfMaps!=NULL);				
			m_pTheFinalMap = m_pCollectionOfMaps->GetBestMap();
		}
		else if (extrapolationType=="GMDS")		// 
		{
			assert(m_pFinalCoarseMap!=NULL);
			
			double interpEpsilon = m_Params.GetDoubleValue(algorithmName, "GMDSInterpEpsilon", valid);
			if (!valid)
				interpEpsilon = 0;
			m_pTheFinalMap = new MapGeoFeature(GetSurface(0), GetSurface(1), interpEpsilon, 
											   m_pFinalCoarseMap, "default", "default", 
											   new SurfaceSampleSet(GetSurface(0)->GetMesh()),
											   new SurfaceSampleSet(GetSurface(1)->GetMesh()));	
		}
		else if (extrapolationType=="MultiConformalSmartTakeClosest")
		{
			InterpolateCoarseMultiConfSmartTakeClosest(algorithmName);
		}
		else if (extrapolationType=="MultiConformal")
		{
			VKString tripletsFilename = m_Params.GetStrValue(algorithmName, "Triplets", valid);
			std::vector<int> triplets;
			std::vector<double> weights;
			if (tripletsFilename!="none" && valid)
			{
				std::ifstream textStream(tripletsFilename.c_str());
				assert(textStream.is_open());
				while (!textStream.eof())
				{
					int t1, t2, t3;
					double w;
					
					textStream>>t1>>t2>>t3>>w;
										
					std::cout<<"Loaded Triplet=["<<t1<<", "<<t2<<", "<<t3<<"] w="<<w<<std::endl;
					if (w > 0)
					{
						triplets.push_back(t1);
						triplets.push_back(t2);
						triplets.push_back(t3);
						weights.push_back(w);
					}
				}
			}
			
			// Fill m_pMapCollection with multi-conformal maps
			int maxConformalMaps = m_Params.GetIntValue(algorithmName, "MaxConformalMaps", valid);
			assert(valid);
			VKString multiConfConfidence = m_Params.GetStrValue(algorithmName, "MultiConfWeights", valid);
			assert(valid);
			VKString multiConfInterp = m_Params.GetStrValue(algorithmName, "MultiConfInterp", valid);
			assert(valid);
			AddLoadedCoarse(algorithmName, extrapolationType);
			MapMultiConformal * newBIMMap = NULL;
			if (triplets.size()==0)
				newBIMMap = new MapMultiConformal(m_pFinalCoarseMap, m_sCoarseSetName, multiConfConfidence, 
													   multiConfInterp, maxConformalMaps);
			else
				newBIMMap = new MapMultiConformal(m_pFinalCoarseMap, m_sCoarseSetName, multiConfConfidence, 
													   multiConfInterp, triplets, weights);

			double geoOutlier = m_Params.GetDoubleValue(algorithmName, "GeoCtrOutlierThreshold", valid);
			if (valid)
				newBIMMap->SetGeodesicCentroidOutlierThreshold(geoOutlier);
			m_pTheFinalMap = newBIMMap;
		}			
		else
		{
			std::cout<<"[ERROR] Unknown extrapolation method: "<<extrapolationType.c_str()<<std::endl;
			assert(false);
		}
		
		assert(m_pTheFinalMap!=NULL);
		m_pTheFinalMap->SaveMap(LocalFilename(finalMapExt));
	}
	else
	{		// if set was created by extrapolation - might need to reaload the coarse set
		AddLoadedCoarse(algorithmName, extrapolationType);
	}
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("AABasicCreateMapsWeights");
	
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("SurfaceCorrespondenceConfMap");	
}

int PipelineExploreConfmaps::AddLoadedCoarse(const VKString & algorithmName, 
											 const VKString & extrapolationType)
{
	if (extrapolationType!="MultiConformal" && extrapolationType!="MultiConformalSmartTakeClosest")
		return -1;
	
	if (m_sCoarseSetName!="none")
		return -1;
	
	double ADD_DIST_THRESH = .05;
		
	m_sCoarseSetName = "LoadedCoarseSet";
	const SurfaceSampleSet * domain = m_pFinalCoarseMap->GetValidDomain();
	const SurfaceSampleSet * range = m_pFinalCoarseMap->GetValidRange();	
	SurfaceSampleSet * domainCopy = new SurfaceSampleSet();
	SurfaceSampleSet * rangeCopy = new SurfaceSampleSet();
	for (int i=0; i<domain->NumSamples(); i++)
		domainCopy->AddSample(domain->GetSample(i));
	for (int i=0; i<range->NumSamples(); i++)
		rangeCopy->AddSample(range->GetSample(i));
	
	if (extrapolationType != "MultiConformal")
	{
		SurfaceDistance * distance1 = NULL;
		SurfaceDistance * distance2 = NULL;
		double distNorm1 = sqrt(m_pFinalCoarseMap->GetSurface(0)->Area());
		double distNorm2 = sqrt(m_pFinalCoarseMap->GetSurface(1)->Area());
		if (ADD_DIST_THRESH>0)
		{
			distance1 = m_pFinalCoarseMap->GetSurface(0)->GetOnTheFlyDistanceMetric(ADD_DIST_THRESH, "default", -1);
			distance2 = m_pFinalCoarseMap->GetSurface(1)->GetOnTheFlyDistanceMetric(ADD_DIST_THRESH, "default", -1);
			assert(distance1!=NULL && distance2!=NULL);
		}
		
		assert(extrapolationType.lastIndexOf("Smart")!=-1);	// an extrapolation with additional stuff
		VKString featureSamplesInit = m_Params.GetStrValue(algorithmName, "SamplesFirstStepMCN", valid);
		SurfaceSampleSet * featureSampleSet1 = GetSurface(0)->GetSampleSet(featureSamplesInit);
		SurfaceSampleSet * featureSampleSet2 = GetSurface(1)->GetSampleSet(featureSamplesInit);
		if (featureSamplesInit!="none" && valid)
		{
			assert(featureSampleSet1!=NULL && featureSampleSet2!=NULL);
			
			domainCopy->Union(*featureSampleSet1, ADD_DIST_THRESH, distance1, distNorm1);
			rangeCopy->Union(*featureSampleSet2, ADD_DIST_THRESH, distance2, distNorm2);		
		}
		
		VKString featureSamplesName = m_Params.GetStrValue(algorithmName, "SamplesMobVoteCast", valid);
		featureSampleSet1 = GetSurface(0)->GetSampleSet(featureSamplesName);
		featureSampleSet2 = GetSurface(1)->GetSampleSet(featureSamplesName);
		assert(featureSampleSet1!=NULL && featureSampleSet2!=NULL);
				
		domainCopy->Union(*featureSampleSet1, ADD_DIST_THRESH, distance1, distNorm1);
		rangeCopy->Union(*featureSampleSet2, ADD_DIST_THRESH, distance2, distNorm2);					
	}
	GetSurface(0)->AddSampleSet(m_sCoarseSetName, domainCopy);
	GetSurface(1)->AddSampleSet(m_sCoarseSetName, rangeCopy);
//	std::cout<<"[DELETEME] Domain / Range size: "<<domainCopy->NumSamples()<<" / "<<rangeCopy->NumSamples()<<std::endl;
	return domain->NumSamples();
}

void PipelineExploreConfmaps::InterpolateCoarseMultiConfSmartTakeClosest(const VKString & algorithmName)
{		
	int MAX_COUNTER=15;
	bool SET_NEW_AS_KNOWN = true;
	double MAP_ONLY_IF_CONFIDENCE_THRESHOLD = 0;
	double VALID_MAP_MULTICONFORMAL_THRESHOLD = .2;
	
	VKString initialMCNSetName = m_Params.GetStrValue(algorithmName, 
													  "SamplesFirstStepMCN", valid);
	SurfaceSampleSet * initialMCNSet1 = NULL;
	SurfaceSampleSet * initialMCNSet2 = NULL;
	if (valid && initialMCNSetName!="none")
	{
		initialMCNSet1 = GetSurface(0)->GetSampleSet(initialMCNSetName);
		initialMCNSet2 = GetSurface(1)->GetSampleSet(initialMCNSetName);
		assert(initialMCNSet1!=NULL && initialMCNSet2!=NULL);
	}
	VKString relaxedMCNProblemSet = m_Params.GetStrValue(algorithmName, 
														 "SamplesMobVoteCast", valid);
	SurfaceSampleSet * iterativeMCNSet1 = GetSurface(0)->GetSampleSet(relaxedMCNProblemSet);
	SurfaceSampleSet * iterativeMCNSet2 = GetSurface(1)->GetSampleSet(relaxedMCNProblemSet);
	assert(iterativeMCNSet1!=NULL && iterativeMCNSet2!=NULL);
	
	VKString extrapolationType = m_Params.GetStrValue(algorithmName, 
													  "ExtrapolationMethod", valid);	
	assert(extrapolationType=="MultiConformalSmartTakeClosest");
	int knownCorrs = AddLoadedCoarse(algorithmName, extrapolationType);
	assert(knownCorrs>0);
	// Fill m_pMapCollection with multi-conformal maps
	VKString multiConfConfidence = m_Params.GetStrValue(algorithmName, 
														"MultiConfWeights", valid);
	assert(valid);
	VKString multiConfInterp = m_Params.GetStrValue(algorithmName, "MultiConfInterp", valid);
	assert(valid);
	int maxConformalMaps = m_Params.GetIntValue(algorithmName, "MaxConformalMaps", valid);
	assert(valid);

	MapMultiConformal * currMap = new MapMultiConformal(m_pFinalCoarseMap, m_sCoarseSetName,
														multiConfConfidence, multiConfInterp, 
														maxConformalMaps);

	double geoOutlier = m_Params.GetDoubleValue(algorithmName, "GeoCtrOutlierThreshold", valid);
	if (valid)
		currMap->SetGeodesicCentroidOutlierThreshold(geoOutlier);
	
	MapConfidenceMultiConf * currMapConfidence = (MapConfidenceMultiConf*)((SurfaceMap*)currMap)->GetMapConfidenceCalculator(multiConfConfidence);
	currMapConfidence->SetIfNoGoodConformalMapThreshold(VALID_MAP_MULTICONFORMAL_THRESHOLD);
	
	SurfaceSampleSet * coarseSet1 = GetSurface(0)->GetSampleSet(m_sCoarseSetName);
	SurfaceSampleSet * coarseSet2 = GetSurface(1)->GetSampleSet(m_sCoarseSetName);

	std::vector<int> unknownSamples1IDs;	
	std::vector<int> unknownSamples2IDs;	
	SurfaceSampleSet * unknownSamples1 = new SurfaceSampleSet();
	SurfaceSampleSet * unknownSamples2 = new SurfaceSampleSet();
	DistanceOnTheFly * dist1 = GetSurface(0)->GetOnTheFlyDistanceMetric(-1, "default", -1);
	DistanceOnTheFly * dist2 = GetSurface(1)->GetOnTheFlyDistanceMetric(-1, "default", -1);
	for (int i=knownCorrs; i<coarseSet1->NumSamples(); i++)
		dist1->PrecomputeRow(coarseSet1->GetSample(i).NearestVertex());
	for (int i=knownCorrs; i<coarseSet2->NumSamples(); i++)
		dist2->PrecomputeRow(coarseSet2->GetSample(i).NearestVertex());
	
	std::map<int, int> originalSampleMap;
	std::map<int, bool> originalSampleOnSurface2Mapped;
	
	std::map<int, bool> sampleOnSurface2Mapped;
	std::map<int, int> knownSamplesMap;
	for (int i=0; i<knownCorrs; i++)
	{
		originalSampleMap[i] = i;
		originalSampleOnSurface2Mapped[i] = true;
	}
	knownSamplesMap = originalSampleMap;
	sampleOnSurface2Mapped = originalSampleOnSurface2Mapped;
	
	std::map<int, int> previousUnknownMap;
	
	bool changedCorrs = true;
	int counter = 0;

	delete m_pFinalCoarseMap;
	m_pFinalCoarseMap = new MapCoarse(GetSurface(0), GetSurface(1), m_sCoarseSetName, m_sCoarseSetName);

	VKString confidenceName = m_Params.GetStrValue(algorithmName, "RefinementMapEval", valid);	// area distor
	//VKString confidenceName = m_Params.GetStrValue(algorithmName, "FineMapEvaluator", valid); // 'multi-conf err'
	while (changedCorrs && counter++ < MAX_COUNTER)
	{
		changedCorrs = false;
		std::cout<<"ITERATION = "<<(counter-1)<<" NUM_KNOWN="<<knownSamplesMap.size()<<std::endl;
		
		SurfaceMapConfidence * confidence = currMap->GetMapConfidenceCalculator(confidenceName);
		
		// set unknown samples
		unknownSamples1->ClearSamples();
		unknownSamples2->ClearSamples();
		for (int i=0; i<coarseSet1->NumSamples(); i++)
		{
			std::map<int, int>::iterator iter1 = knownSamplesMap.find(i);			
			if (iter1==knownSamplesMap.end())
			{
				bool confident = true;
				bool apropriateSet = true;
				if (MAP_ONLY_IF_CONFIDENCE_THRESHOLD>0)
				{
					double confVal = 1.;
					assert(confidence!=NULL);
					confVal = confidence->ConfidenceAtSampleNoSmoothing(coarseSet1->GetSample(i));
					confident = (confVal>= MAP_ONLY_IF_CONFIDENCE_THRESHOLD);
				}
				
				if (initialMCNSet1!=NULL && initialMCNSet2!=NULL)
				{
					if (counter==1)
						apropriateSet = initialMCNSet1->ContainsExact(coarseSet1->GetSample(i));
					else
						apropriateSet = iterativeMCNSet1->ContainsExact(coarseSet1->GetSample(i));
				}
				
				if (confident && apropriateSet)
					unknownSamples1->AddSample(coarseSet1->GetSample(i));
			}
		}
		for (int i=0; i<coarseSet2->NumSamples(); i++)
		{
			bool appropriateSet = true;
			if (initialMCNSet1!=NULL && initialMCNSet2!=NULL)
			{
				if (counter==1)
					appropriateSet = initialMCNSet2->ContainsExact(coarseSet2->GetSample(i));
				else
					appropriateSet = iterativeMCNSet2->ContainsExact(coarseSet2->GetSample(i));
			}
						
			std::map<int, bool>::iterator iter2 = sampleOnSurface2Mapped.find(i);
			if (iter2==sampleOnSurface2Mapped.end() && appropriateSet)
				unknownSamples2->AddSample(coarseSet2->GetSample(i));
		}
		std::cout<<"Unknown Samples: "<<unknownSamples1->NumSamples()<<" and "<<unknownSamples2->NumSamples()<<std::endl;
		if (unknownSamples1->NumSamples()==0 || unknownSamples2->NumSamples()==0)
			break;
		
		// find mutually closest
		MapCoarse * newCorrs = NULL;
		if (initialMCNSet1!=NULL && initialMCNSet2!=NULL && counter==1)
			newCorrs = new MapCoarse(unknownSamples1, unknownSamples2, currMap, 
									 MapCoarse::F2C_MUTUALLY_CLOSEST_NEIGHBORS);
		else
			newCorrs = new MapCoarse(unknownSamples1, unknownSamples2, currMap, 
									 MapCoarse::F2C_MCN_ONLY_ON_SURF2_FORWARD);
		assert(newCorrs!=NULL);
		// set all mutually closest to 'known' set
		for (int i=0; i < coarseSet1->NumSamples(); i++)
		{
			SurfaceSample sample1 = coarseSet1->GetSample(i);
			SurfaceSample sample2 = newCorrs->ForwardMap(sample1);
			if (!sample2.Invalid())
			{
				int id2 = -1;
				bool hasExact = coarseSet2->ContainsExact(sample2, &id2);
				assert(hasExact && id2>=0);
				knownSamplesMap[i] = id2;
				
				std::map<int, int>::iterator iter1 = previousUnknownMap.find(i);

				if (iter1==previousUnknownMap.end() || iter1->second!=id2)
					changedCorrs = true;					
				
				std::cout<<"KnownMap["<<i<<"] -> "<<id2<<std::endl;
			}
		}
		previousUnknownMap.clear();
		for (std::map<int, int>::iterator iter1 = knownSamplesMap.begin(); iter1!=knownSamplesMap.end(); iter1++)
		{
			previousUnknownMap[iter1->first] = iter1->second;
			sampleOnSurface2Mapped[iter1->second] = true;
		}
		
		// create coarse map based on MCN
		std::vector<int> currCorrs;
		m_pFinalCoarseMap->ClearAllCorrespondences();
		for (int i=0; i<coarseSet1->NumSamples(); i++)
		{
			std::map<int, int>::iterator iter1 = knownSamplesMap.find(i);
			
			if (iter1 != knownSamplesMap.end())
				currCorrs.push_back(iter1->second);
			else
				currCorrs.push_back(-1);
		}
		m_pFinalCoarseMap->SetFinalCorrMap(currCorrs);
		
		// create multi-conformal map
		delete currMap;
		currMap = new MapMultiConformal(m_pFinalCoarseMap, m_sCoarseSetName, multiConfConfidence, multiConfInterp, maxConformalMaps);
		MapConfidenceMultiConf * currMapConfidence = (MapConfidenceMultiConf*)((SurfaceMap*)currMap)->GetMapConfidenceCalculator(multiConfConfidence);
		currMapConfidence->SetIfNoGoodConformalMapThreshold(VALID_MAP_MULTICONFORMAL_THRESHOLD);
		
		double geoOutlier = m_Params.GetDoubleValue(algorithmName, "GeoCtrOutlierThreshold", valid);
		if (valid)
			currMap->SetGeodesicCentroidOutlierThreshold(geoOutlier);		
		
		// set only input to be known (if necessary)
		if (!SET_NEW_AS_KNOWN)
		{
			knownSamplesMap.clear();
			sampleOnSurface2Mapped.clear();
			knownSamplesMap = originalSampleMap;
			sampleOnSurface2Mapped = originalSampleOnSurface2Mapped;
		}
		delete newCorrs;
	}
	std::cout<<"TERMINATED = "<<counter<<" NUM_KNOWN="<<knownSamplesMap.size()<<std::endl;
	assert(currMap!=NULL);
	m_pTheFinalMap = currMap;
}

void PipelineExploreConfmaps::WritePostFinalProcessingResults()
{
	std::cout<<"WRITING POST-PROCESSING STUFF"<<std::endl;
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);
	assert(valid);
	VKString surfDescr = m_Params.GetStrValue(algorithmName, "Surface", valid);
	assert(valid);
	VKStringList meshNames = m_Params.GetStrValues(surfDescr, "MeshName", valid);
	if (meshNames.count()==1)
		meshNames.push_back(meshNames[0]);
	assert(meshNames.count()==2);	// for now correspondence between two only		
	
	VKString algOutputName = m_Params.GetStrValue(algorithmName, "OutputName", valid);	
	if (!valid)
		algOutputName = algorithmName;
	
	VKStringList postFin = m_Params.GetStrValues("Pipeline", "PostFinal", valid);

	if (postFin.contains("CacheVtoV") || postFin.contains("DenseMap")
		|| postFin.contains("DenseSHREC10Map") || postFin.contains("DenseMapPreceise"))
	{
		if (!m_pTheFinalMap->IsCachedVertexToVertex())
		{
			std::cout<<"[WARNING] Caching vertex to vertex map"<<std::endl;
			m_pTheFinalMap->CacheVertexToVertexMap();
		}
	}
		
	if (postFin.contains("DenseMap"))
	{
		VKString denseMapName = LocalFilename(meshNames[0], algOutputName+"_DenseMap_"+meshNames[0]+"_to_"+meshNames[1], ".dense.gaps.map");
		//SurfaceSampleSet * allVertices = new SurfaceSampleSet(GetSurface(0)->GetMesh());
		//MapCoarse * denseMap = new MapCoarse(allVertices, NULL, m_pTheFinalMap, 
		//									 MapCoarse::F2C_FORWARD_ONLY);
		//denseMap->SaveMapInVertexByVertexFormat(denseMapName);
		std::ofstream outStream(denseMapName.c_str());
		assert(outStream.is_open());
		for (int i=0; i<GetSurface(0)->GetMesh()->NVertices(); i++)
		{
			SurfaceSample samp(i, GetSurface(0)->GetMesh());
			SurfaceSample sampOut = m_pTheFinalMap->ForwardMap(samp);
			outStream<<sampOut.NearestVertex()<<"\n";
		}
	}
	
	if (postFin.contains("DenseMapPrecise"))
	{
		std::cout<<"Writing Dense Map"<<std::endl;
		VKString denseMapName = LocalFilename(meshNames[0], algOutputName+"_DenseMapPreceise_"+meshNames[0]+"_to_"+meshNames[1], ".dense.preceise.map");
		std::ofstream outStream(denseMapName.c_str());
		assert(outStream.is_open());
		for (int i=0; i<GetSurface(0)->GetMesh()->NVertices(); i++)
		{
			SurfaceSample samp(i, GetSurface(0)->GetMesh());
			SurfaceSample sampOut = m_pTheFinalMap->ForwardMap(samp);
			outStream<<sampOut.TriID()<<" "<<sampOut.B(0)<<" "<<sampOut.B(1)<<" "<<sampOut.B(2)<<"\n";
		}
	}
	
	if (postFin.contains("DenseMapPreciseAll"))
	{
		if (m_MapIDToEigenvalue.size()>0)
		{
			std::cout<<"Writing Dense Map"<<std::endl;
			VKString fileConfidences = algOutputName+"_DenseMapPreceiseAll_";
			fileConfidences += VKString("_")+meshNames[0]+"_to_"+meshNames[1]+ ".dense.preceise.confidences.txt";
			std::ofstream outStreamConfidences(LocalFilename(meshNames[0], fileConfidences).c_str());
			assert(outStreamConfidences.is_open());
			assert(m_pCollectionOfMaps!=NULL);
			int nMaps = m_pCollectionOfMaps->GetNumMaps();
			outStreamConfidences<<"NumMaps "<<nMaps<<"\n";
			
			for (int i=0; i<nMaps; i++)
			{
				double score;
				SurfaceMap * currMap = m_pCollectionOfMaps->GetMapByID(i, &score);
				outStreamConfidences<<score<<" ";
				double eval = -1;
				if (m_pEigenvaluesIBM!=NULL)
				{
					assert(i < (int)m_MapIDToEigenvalue.size());
					int eigID = m_MapIDToEigenvalue[i].id;
					assert(eigID < m_pEigenvaluesIBM->Rows());
					eval = (*m_pEigenvaluesIBM)(eigID);
				}
				outStreamConfidences<<eval<<" ";
				
				VKString denseMapName = algOutputName+"_DenseMapPreceiseAll_" + VKString::number(i);
				denseMapName += VKString("_")+meshNames[0]+"_to_"+meshNames[1]+ ".dense.preceise.map";
				
				std::ofstream outStreamDense(LocalFilename(meshNames[0], denseMapName).c_str());
				assert(outStreamDense.is_open());
				for (int i=0; i<GetSurface(0)->GetMesh()->NVertices(); i++)
				{
					SurfaceSample samp(i, GetSurface(0)->GetMesh());
					SurfaceSample sampOut = currMap->ForwardMap(samp);
					outStreamDense<<sampOut.TriID()<<" "<<sampOut.B(0)<<" "<<sampOut.B(1)<<" "<<sampOut.B(2)<<"\n";
				}
			}
		}
		else
		{
			std::cout<<"[WARNING] Could not write dense preceise maps: eigenvalues are missing";
			std::cout<<std::endl;
		}
	}
	
	if (postFin.contains("CoarseMapWithConfidencesAll"))
	{
		assert(false);	// writes indefinitely
		VKString fileConfidences = algOutputName+"_MapWithConfidencesAll_";
		fileConfidences += VKString("_")+meshNames[0]+"_to_"+meshNames[1]+ ".coarse.maps.confidences.txt";
		std::ofstream outStreamConfidences(LocalFilename(meshNames[0], fileConfidences).c_str());
		assert(outStreamConfidences.is_open());
		assert(m_pCollectionOfMaps!=NULL);
		int nMaps = m_pCollectionOfMaps->GetNumMaps();
		outStreamConfidences<<"NumMaps "<<nMaps<<"\n";

		VKString sampleSetName = m_Params.GetStrValue("Pipeline", "EvalOnPntSet", valid);
		SurfaceSampleSet * sampleSet1 = NULL;
		std::ifstream textStreamInSet(sampleSetName.c_str());
		if (!valid || sampleSetName == "none")
		{
			VKString sampleSetName = m_Params.GetStrValue(algorithmName, "SamplesFineCorr", valid);	
			assert(valid);
			sampleSet1 = GetSurface(0)->GetSampleSet(sampleSetName);	
			assert(sampleSet1!=NULL);		
		}
		else
		{
			assert(textStreamInSet.is_open());
		}
		
		VKString confName = m_Params.GetStrValue(algorithmName, "RefinementMapEval", valid);
		assert(valid);
		
		for (int i=0; i<nMaps; i++)
		{
			double score;
			SurfaceMap * currMap = m_pCollectionOfMaps->GetMapByID(i, &score);
			outStreamConfidences<<score<<" ";
			double eval = -1;
			if (m_pEigenvaluesIBM!=NULL)
			{
				assert(i < (int)m_MapIDToEigenvalue.size());
				int eigID = m_MapIDToEigenvalue[i].id;
				assert(eigID < m_pEigenvaluesIBM->Rows());
				eval = (*m_pEigenvaluesIBM)(eigID);
			}
			outStreamConfidences<<eval<<" ";
			
			VKString denseMapName = algOutputName+"_MapWithConfidencesAll_" + VKString::number(i);
			denseMapName += VKString("_")+meshNames[0]+"_to_"+meshNames[1]+ ".coarse.pnts.map";
			std::ofstream textStreamOutCorrsConfidences(denseMapName.c_str());
			
			SurfaceMapConfidence * confCalc = currMap->GetMapConfidenceCalculator(confName);
			
			int pID = 0;
			while (true)
			{
				SurfaceSample samp;
				if (sampleSet1==NULL)
				{
					samp = sampleSet1->GetSample(pID);
					pID++;
					if (pID >= sampleSet1->NumSamples())
						break;;
				}
				else 
				{
					int vID;
					textStreamInSet>>vID;
					if (textStreamInSet.eof())
						break;
					samp = SurfaceSample(vID, GetSurface(0)->GetMesh());
				}
				
				
				int vID = samp.NearestVertex();
				
				SurfaceSample sampOut = currMap->ForwardMap(samp);
				double val = confCalc->CalculateConfidenceAtSample(samp);
				
				textStreamOutCorrsConfidences<<"vID "<<vID<<"\n";
				textStreamOutCorrsConfidences<<"Corr "<<sampOut.TriID()<<" ";
				textStreamOutCorrsConfidences<<sampOut.B(0)<<" "<<sampOut.B(1)<<" "<<sampOut.B(2)<<"\n";
				textStreamOutCorrsConfidences<<"Val "<<val<<"\n\n";
			}
		}
	}

	
	if (postFin.contains("SymmetryFunction"))
	{
		assert(false);		/// apparently this is not implemented right
		WriteSymmetryFunction(LocalFilename(meshNames[0], meshNames[0], ".symmetry.arff"), 
							  VKString("../../")+meshNames[0]+".off");
	}

	
	
	bool validBenchmark = true;
	VKString benchmarkName = m_Params.GetStrValue("Pipeline", "BenchmarkQuery", valid);	
	validBenchmark = validBenchmark && valid;
	VKStringList benchPath = m_Params.GetStrValues(benchmarkName, "Path", valid);	
	validBenchmark = validBenchmark && valid && benchPath.count()==2;	
	VKStringList benchFile = m_Params.GetStrValues(benchmarkName, "FileVertexIDs", valid);
	validBenchmark = validBenchmark && valid && benchFile.count()==2;	
	
	VKString benchMeshFromNull;
	VKString benchMeshFrom;
	if (validBenchmark)
	{
		benchMeshFromNull = benchPath[0] + benchFile[0];
		benchMeshFrom = benchPath[1] + benchFile[1];
	}
	
	if (postFin.contains("DenseSHREC10Map"))
	{
		assert(validBenchmark);
		VKString filebody = VKString("SHREC10Corr_Dense_")+meshNames[1];
		VKString locFilename = LocalFilename(meshNames[0], filebody, ".corr");
		WriteDenseSHREC10(benchMeshFromNull, benchMeshFrom, locFilename);
	}
	
	for (int i=0; i<postFin.count(); i++)
	{
		if (postFin[i].startsWith("SparseSHREC10Map"))
		{
			assert(validBenchmark);				
			VKString tempStrCopy = postFin[i];
			VKString numPntsStr = tempStrCopy.replace("SparseSHREC10Map_", "");
			bool ok;
			int numPnts = numPntsStr.toInt(&ok);
			assert(ok);
			VKString filebody = VKString("SHREC10Corr_Sparse_")+VKString::number(numPnts)+"_"+meshNames[1];
			VKString locFilename = LocalFilename(meshNames[0], filebody, ".corr");
			
			VKString confidenceName = m_Params.GetStrValue(algorithmName, "FineMapEvaluator", valid);
			assert(valid);
			VKString sampleSetName = m_Params.GetStrValue(algorithmName, "SamplesForFineErr", valid);
			assert(valid);
			
			WriteSparseSHREC10(benchMeshFromNull, benchMeshFrom, locFilename, numPnts, 
							   sampleSetName, confidenceName);
		}
	}
	
	VKStringList writeConfidences = m_Params.GetStrValues("Pipeline", "WriteFinalConfidence", valid);
	for (int i=0; i<writeConfidences.count(); i++)
	{
		assert(m_pTheFinalMap!=NULL);
		VKString confName = writeConfidences[i];
		SurfaceMapConfidence * confidence = m_pTheFinalMap->GetMapConfidenceCalculator(confName);
		VKString locFilename = LocalFilename(meshNames[0], algOutputName+"_"+confName, ".conf.vals");
		confidence->Confidence();
		//confidence->Error();	
		confidence->ConfidenceAtSample(SurfaceSample(0, GetSurface(0)->GetMesh()));
		std::ofstream textStream(locFilename.c_str());
		assert(textStream.is_open());
		confidence->SaveConfidenceValues(textStream);
		textStream.close();
	}
	
	if (postFin.contains("FuzzyMaps"))
	{
		WriteFuzzymaps();
	}
}

void PipelineExploreConfmaps::WriteFuzzymaps()
{
	int nMaps = m_pCollectionOfMaps->GetNumMaps();
	if (nMaps==0)
		return;
	
//	std::cout<<"nMaps="<<nMaps<<std::endl;
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);
	assert(valid);
	VKString surfDescr = m_Params.GetStrValue(algorithmName, "Surface", valid);
	assert(valid);
	VKStringList meshNames = m_Params.GetStrValues(surfDescr, "MeshName", valid);
	assert(valid);
	if (meshNames.count()==1)
		meshNames.push_back(meshNames[0]);
	
//	std::cout<<"fuzzy - 1"<<std::endl;
	VKString confName = m_Params.GetStrValue(algorithmName, "RefinementMapEval", valid);
	assert(valid);
	assert(meshNames.count()==2);	// for now correspondence between two only		
	VKString distMetricName=m_Params.GetStrValue(surfDescr, "DefaultMetric", valid);
	assert(valid);

	VKString filebody = VKString("FuzzyMap_") + meshNames[0] + "_to_"+meshNames[1];
	VKString algOutputName = m_Params.GetStrValue(algorithmName, "OutputName", valid);	
	if (!valid)
		algOutputName = algorithmName;

//	std::cout<<"fuzzy - 2"<<std::endl;	
	double maxDist = m_Params.GetDoubleValue("Pipeline", "FuzzyDistance", valid);
	assert(valid);
	double sigma2 = pow(maxDist, 2.);
	
	VKString sampleSetName = m_Params.GetStrValue(algorithmName, "SamplesFineCorr", valid);	
	SurfaceSampleSet * sampleSet1 = GetSurface(0)->GetSampleSet(sampleSetName);
	assert(sampleSet1!=NULL);
	SurfaceSampleSet * sampleSet2 = GetSurface(1)->GetSampleSet(sampleSetName);
	assert(sampleSet2!=NULL);	
	
//	std::cout<<"fuzzy - 3"<<std::endl;		
	// write Fuzzy.txt (meta file)
	VKString fuzzyFile = LocalFilename(meshNames[0], filebody, ".txt");
	std::ofstream textStream(fuzzyFile.c_str());
	assert(textStream.is_open());
	textStream<<"NumFuzzyMaps "<<nMaps<<"\n";
	SurfaceDistance * dist2 = GetSurface(1)->GetOnTheFlyDistanceMetric(maxDist,
																	  distMetricName, 256);
	assert(dist2!=NULL);
	VKString mapFile = LocalFilename(meshNames[0], filebody, ".concat.bin");
	FILE * fuzzyCorrFile = fopen(mapFile.c_str(), "wb");
	assert(fuzzyCorrFile!=NULL);
	fwrite(&nMaps, sizeof(int), 1, fuzzyCorrFile);
	
//	std::cout<<"fuzzy - 4"<<std::endl;			
	for (int i=0; i<nMaps; i++)
	{
//		std::cout<<"fuzzy - ALLMAPS. "<<i<<" / "<<nMaps<<std::endl;					
		double score=0;
		SurfaceMap * currMap = m_pCollectionOfMaps->GetMapByID(i, &score);
//		VKString mapID = VKString(".map")+VKString::number(i)+".fuzzy";
//		VKString mapFile = LocalFilename(meshNames[0], filebody, mapID);
//		textStream<<mapFile.c_str()<<"\n"<<score<<"\n";
		textStream<<score<<"\n";
		SurfaceMapConfidence * confidence = currMap->GetMapConfidenceCalculator(confName);
		assert(confidence!=NULL);
//		std::cout<<"fuzzy - ALLMAPS - 1"<<std::endl;
		
		int nModels=1;
		fwrite(&nModels, sizeof(int), 1, fuzzyCorrFile);
		int model1ID=0;
		int model2ID=1;
		int nPnts1=sampleSet1->NumSamples();
		int nPnts2=sampleSet2->NumSamples();
		fwrite(&model1ID, sizeof(int), 1, fuzzyCorrFile);
		fwrite(&nPnts1, sizeof(int), 1, fuzzyCorrFile);
		
//		std::cout<<"fuzzy - ALLMAPS - 2"<<std::endl;
		bool atVtx;		
		for (int p1ID=0; p1ID < nPnts1; p1ID++)
		{
//			std::cout<<"\tMapping pnt: "<<p1ID<<"/ "<<nPnts1<<std::endl;
			SurfaceSample samp1 = sampleSet1->GetSample(p1ID);
			std::map<int, double> corrs;
//			std::cout<<"\t\tsampConf="<<samp1.NearestVertex()<<std::endl;
			double c = confidence->ConfidenceAtSample(samp1);
//			std::cout<<"\t\tsampForwMap="<<samp1.NearestVertex()<<std::endl;
			SurfaceSample samp2 = currMap->ForwardMap(samp1);
			
//			std::cout<<"\t\tForwardMap:"<<samp2.NearestVertex()<<std::endl;
//			double mind=10000;
//			double maxc=0;
			for (int p2ID=0; p2ID < nPnts2; p2ID++)
			{
//				std::cout<<"\t\t\tExamine p2: "<<p2ID<<" / "<<nPnts2<<std::endl;
				int v2ID = sampleSet2->GetSample(p2ID).NearestVertex(&atVtx);
//				std::cout<<"\t\t\tExamine p2 : 2"<<std::endl;
				assert(atVtx);
				double d = dist2->Distance(samp2, sampleSet2->GetSample(p2ID));
//				std::cout<<"\t\t\tExamine p2 : 3"<<std::endl;				
				if (d < maxDist)
				{
					corrs[v2ID] = c * exp(-d / sigma2);
//					if (mind>d)
//					{
//						mind = d;
//						maxc = c;
//					}
				}
			}
			
//			double v = maxc * exp(-mind / sigma2);
			//std::cout<<"\t\tCorrVal: "<<mind<<" "<<maxc<<" "<<sigma2<<" = "<<v<<std::endl;
			
//			std::cout<<"\t\tFuzzyMap:"<<corrs.size()<<std::endl;
			
			int nCorrs = (int)corrs.size();
			int v1ID = samp1.NearestVertex(&atVtx);
			assert(atVtx);
			fwrite(&v1ID, sizeof(int), 1, fuzzyCorrFile);
			fwrite(&nCorrs, sizeof(int), 1, fuzzyCorrFile);
			for (std::map<int, double>::iterator iter = corrs.begin(); 
				 iter != corrs.end(); iter++)
			{
				fwrite(&model2ID, sizeof(int), 1, fuzzyCorrFile);
				fwrite(&iter->first, sizeof(int), 1, fuzzyCorrFile);
				fwrite(&iter->second, sizeof(double), 1, fuzzyCorrFile);
			}
//			std::cout<<"\t-----"<<std::endl;
		}
//		std::cout<<"fuzzy - ALLMAPS - DONE"<<std::endl;		
		int checksum = -29349234;
		fwrite(&checksum, sizeof(int), 1, fuzzyCorrFile);
	}
	fclose(fuzzyCorrFile);
	textStream.close();
	
	VKString doneFile = LocalFilename(meshNames[0], filebody, ".done.txt");
	std::ofstream doneFileStream(doneFile.c_str());
	assert(doneFileStream.is_open());
	doneFileStream<<"DONE!";
	doneFileStream.close();
	
}

void PipelineExploreConfmaps::WriteDenseSHREC10(const VKString & benchMeshFromNull,
												const VKString & benchMeshToXform,
												const VKString & outputCorrFile)
{
	assert(m_pTheFinalMap!=NULL);
	// load files
	int numFromVertices=-1;
	int numToVertices=-1;
	SHREC10CorrVertex * vertexFromArray=NULL;
	SHREC10CorrVertex * vertexToArray=NULL;
	
	ReadSHREC10CorrVerticesIntoKDTree(benchMeshFromNull, &numFromVertices, &vertexFromArray);
	
	R3Kdtree<SHREC10CorrVertex*> * toKdTree;
	ReadSHREC10CorrVerticesIntoKDTree(benchMeshToXform, &numToVertices, &vertexToArray, &toKdTree);
	
	std::ofstream textStream(outputCorrFile.c_str());
	// for every 'from' vertex: 
	for (int i=0; i<numFromVertices; i++)
	{
		// find an appropriate 'to' map: 		
		SHREC10CorrVertex & vFrom=vertexFromArray[i];
		if (vFrom.triID==-1)
			continue;	// ignore isolated vertices
		SurfaceSample s1 = GetSurface(0)->GetNearestSample(vFrom.pnt);
		assert(!s1.Invalid());
		SurfaceSample s2 = m_pTheFinalMap->ForwardMap(s1);
		assert(!s2.Invalid());
		SHREC10CorrVertex * vTo = toKdTree->FindClosest(s2.GetPosition());
		assert (vTo->triID!=-1);	// kd tree should not have isolated vertices
		textStream<<vFrom.triID<<" "<<vFrom.bary[0]<<" "<<vFrom.bary[1]<<" "<<vFrom.bary[2]<<" ";
		textStream<<vTo->triID<<" "<<vTo->bary[0]<<" "<<vTo->bary[1]<<" "<<vTo->bary[2]<<"\n";		
	}
	textStream.close();
}

void PipelineExploreConfmaps::WriteSparseSHREC10(const VKString & benchMeshFromNull,
												 const VKString & benchMeshToXform,
												 const VKString & outputCorrFile,
												 int numSparseCorrs,
												 const VKString & samplesetName,
												 const VKString & confidenceName)
{
	assert(m_pTheFinalMap!=NULL);
	
	// load files
	int numFromVertices=-1;
	int numToVertices=-1;
	SHREC10CorrVertex * vertexFromArray=NULL;
	SHREC10CorrVertex * vertexToArray=NULL;
	R3Kdtree<SHREC10CorrVertex*> * fromKdTree;	
	ReadSHREC10CorrVerticesIntoKDTree(benchMeshFromNull, &numFromVertices, &vertexFromArray, &fromKdTree);
	
	R3Kdtree<SHREC10CorrVertex*> * toKdTree;
	ReadSHREC10CorrVerticesIntoKDTree(benchMeshToXform, &numToVertices, &vertexToArray, &toKdTree);

	// find top sparse corrs
	SurfaceSampleSet * sampleSet = GetSurface(0)->GetSampleSet(samplesetName);
	assert(sampleSet!=NULL);
	SurfaceSampleSet * bestSparse;
	
	SurfaceMapConfidence * confidence = m_pTheFinalMap->GetMapConfidenceCalculator(confidenceName);
	MapCoarse * sparseMap = confidence->GetBestSparseCorrs(sampleSet, numSparseCorrs, &bestSparse);
	
	std::ofstream textStream(outputCorrFile.c_str());	
	// for each sparse corr - find correspondence	
	for (int i=0; i<bestSparse->NumSamples(); i++)
	{
		SurfaceSample s1 = bestSparse->GetSample(i);
		assert(!s1.Invalid());
		SurfaceSample s2 = sparseMap->ForwardMap(s1);
		assert(!s2.Invalid());
		
		SHREC10CorrVertex * vFrom = fromKdTree->FindClosest(s1.GetPosition());
		SHREC10CorrVertex * vTo = toKdTree->FindClosest(s2.GetPosition());

		assert (vFrom->triID!=-1);	// kd tree should not have isolated vertices		
		assert (vTo->triID!=-1);	// kd tree should not have isolated vertices
		textStream<<vFrom->triID<<" "<<vFrom->bary[0]<<" "<<vFrom->bary[1]<<" "<<vFrom->bary[2]<<" ";
		textStream<<vTo->triID<<" "<<vTo->bary[0]<<" "<<vTo->bary[1]<<" "<<vTo->bary[2]<<"\n";		
		
	}
}

void PipelineExploreConfmaps::ReadSHREC10CorrVerticesIntoKDTree(const VKString & benchMesh, 
																int * numVertices,
																SHREC10CorrVertex ** vertexArray,
																R3Kdtree<SHREC10CorrVertex*> ** kdtree)
{
	assert(numVertices!=NULL);
	assert(vertexArray!=NULL);
	std::ifstream textStream(benchMesh.c_str());
	std::string tempStr;
	textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "OFF")==0);
	int numFaces, numEdges;
	textStream>>(*numVertices)>>numFaces>>numEdges;
	*vertexArray = new SHREC10CorrVertex[*numVertices];
	double x, y, z;	
	for (int i=0; i<(*numVertices); i++)
	{
		SHREC10CorrVertex & v = (*vertexArray)[i];
		textStream>>x>>y>>z;
		v.pnt = R3Point(x, y, z);
		v.id = i;
		v.triID = -1;
		v.bary[0] = -1;
		v.bary[1] = -1;
		v.bary[2] = -1;
	}
	
	int v1, v2, v3;	
	for (int i=0; i<numFaces; i++)
	{
		textStream>>v1>>v2>>v3;
		SHREC10CorrVertex & v1SHREC = (*vertexArray)[v1-1];
		SHREC10CorrVertex & v2SHREC = (*vertexArray)[v2-1];
		SHREC10CorrVertex & v3SHREC = (*vertexArray)[v3-1];		
		int triID = i+1;	// NOTE: I am not sure if it's 1-based or 0-based (i vs i+1)
		v1SHREC.triID = triID;	
		v1SHREC.bary[0] = 1.;		v1SHREC.bary[1] = 0;		v1SHREC.bary[2] = 0;
		v2SHREC.triID = triID;
		v2SHREC.bary[0] = 0.;		v2SHREC.bary[1] = 1.;		v2SHREC.bary[2] = 0;
		v3SHREC.triID = triID;		
		v3SHREC.bary[0] = 0.;		v3SHREC.bary[1] = 0;		v3SHREC.bary[2] = 1.;
	}
	
	if (kdtree!=NULL)
	{
		RNArray<SHREC10CorrVertex*> vertexRNArray;
		for (int i=0; i<(*numVertices); i++)
		{
			SHREC10CorrVertex & v = (*vertexArray)[i];
			if (v.triID>=0)
				vertexRNArray.InsertTail(&v);
		}
		*kdtree = new R3Kdtree<SHREC10CorrVertex*>(vertexRNArray, &global_SHREC10CorrVertex, NULL);
	}
}

void PipelineExploreConfmaps::WriteSymmetryFunction(const VKString & symmFnName,
													const VKString & relationName)
{
	std::cout<<"Writing Symmetry Function"<<std::endl;
	assert(m_pFinalCoarseMap!=NULL);
	SurfacePerVertexValue vals(GetSurface(0), false); 
	SurfaceDistance * distance = GetSurface(0)->GetDistanceMetric("default");
	const SurfaceSampleSet * domain = m_pFinalCoarseMap->GetValidDomain();
	for (int i=0; i<domain->NumSamples(); i++)
	{
		SurfaceSample s1 = domain->GetSample(i);
		SurfaceSample s2 = m_pFinalCoarseMap->ForwardMap(s1);
		bool atVertex;
		int vertexID = s1.NearestVertex(&atVertex);
		assert(atVertex);
		if (!s2.Invalid())
			vals.SetVertexValue(vertexID, distance->Distance(s1, s2));
	}
	vals.SmoothGaussianFromKnown(-1);
	vals.WritePerVertexValuesInArff(symmFnName, relationName, "DistanceToSymmetricCorr");
}

void PipelineExploreConfmaps::RescoreCollection(MapScoredCollection * collection, 
												const VKString & confidenceName)
{
	if (confidenceName!="none")
		collection->RecalculateScores(confidenceName);

}

void PipelineExploreConfmaps::LoadCoarseCorrsFromTruth()
{
	m_pFinalCoarseMap = NULL;
	m_pTheFinalMap = NULL;
	m_pCollectionOfMaps = NULL;
	m_pCollectionOfConformalMaps = NULL;
	VKString truthMapName = m_Params.GetStrValue("Pipeline", "BenchmarkQuery", valid);		assert(valid);
	VKString mapType = m_Params.GetStrValue(truthMapName, "MapType", valid);				 assert(valid);
	VKStringList truthPaths = m_Params.GetStrValues(truthMapName, "Path", valid);					 assert(valid);
	if (truthPaths.count()==1)
		truthPaths.push_back(truthPaths[0]);
	VKStringList truthPnts = m_Params.GetStrValues(truthMapName, "FileVertexIDs", valid);	 assert(valid);
	VKString truthCorr = m_Params.GetStrValue(truthMapName, "FileFeaturesMap", valid);		

	assert(mapType=="CoarseMap_SameVertexIDs");	
	MapCoarse humanFeatureMap(GetSurface(0), GetSurface(1));

	if (GetSurface(0)->GetName()==GetSurface(1)->GetName() && truthPnts.count()==1)		// symmetry
		humanFeatureMap.LoadCorrFeaturePntSymmetryMap(truthPaths[0]+truthPnts[0], truthPaths[0]+truthCorr);
	else
		humanFeatureMap.LoadCorrSameFeaturePntIDs(truthPaths[0]+truthPnts[0], truthPaths[1]+truthPnts[1]);	
	
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);	
	int numPnts = m_Params.GetIntValue(algorithmName, "NumFeaturePnts", valid);
	if (!valid || numPnts==-1 || numPnts>humanFeatureMap.GetValidDomain()->NumSamples())
		numPnts = humanFeatureMap.GetValidDomain()->NumSamples();
	std::vector<int> featurePntOrder = m_Params.GetIntValues(algorithmName, "FeaturePntOrder", valid);
	int givenOrdPnts = (int)featurePntOrder.size();
	
	//	initialize list to check if a sample was given in featurePntOrder
	std::vector<bool> givenInOrderList;
	for (int i=0; i<(int)humanFeatureMap.GetValidDomain()->NumSamples(); i++)
		givenInOrderList.push_back(false);	
	for (int i=0; i<(int)featurePntOrder.size(); i++)
		givenInOrderList[featurePntOrder[i]] = true;
	
	// add uninitialized samples to make sure all samples are included
	for (int i=0; i<(int)humanFeatureMap.GetValidDomain()->NumSamples(); i++)
		if (!givenInOrderList[i])
			featurePntOrder.push_back(i);
//	std::cout<<featurePntOrder.size()<<" =? "<<humanFeatureMap.GetValidDomain()->NumSamples()<<std::endl;
	assert((int)featurePntOrder.size()==humanFeatureMap.GetValidDomain()->NumSamples());
	
	// randomize uninitialize samples
	for (int i=givenOrdPnts; i<(int)featurePntOrder.size(); i++)
	{
		int randVal=rand()%(featurePntOrder.size()-givenOrdPnts);
		if (randVal<0)
			randVal *= -1;
		int randSwap = givenOrdPnts + randVal;
		int tempVal = featurePntOrder[randSwap];
		featurePntOrder[randSwap] = featurePntOrder[i];
		featurePntOrder[i] = tempVal;
	}
			
	// verify loaded desired set
	m_sCoarseSetName = "BenchmarkQuery";
	std::vector<int> corrs; 	
	SurfaceSampleSet * benchmarkDomain = GetSurface(0)->GetSampleSet(m_sCoarseSetName);
	SurfaceSampleSet * benchmarkRange = GetSurface(1)->GetSampleSet(m_sCoarseSetName);			
	assert(benchmarkDomain!=NULL);
	assert(benchmarkRange!=NULL);	
	assert(benchmarkDomain->NumSamples()==humanFeatureMap.GetValidDomain()->NumSamples());
	for (int i=0; i<benchmarkDomain->NumSamples(); i++)
	{
		assert(humanFeatureMap.GetValidDomain()->GetSample(i)==benchmarkDomain->GetSample(i));
		corrs.push_back(-1);
	}
	
	// initialize new correspondence
	for (int i=0; i<numPnts; i++)
	{
		int s1ID = featurePntOrder[i];
		SurfaceSample s1 = humanFeatureMap.GetValidDomain()->GetSample(s1ID);
		SurfaceSample s2 = humanFeatureMap.ForwardMap(s1);
		
		int s2ID=-1;
		for (int j=0; j<benchmarkRange->NumSamples(); j++)
		{
			if (s2==benchmarkRange->GetSample(j))
				s2ID = j;
		}
		assert(s2ID!=-1);
		
		corrs[featurePntOrder[i]] = s2ID;
	}	
	m_pFinalCoarseMap = new MapCoarse(GetSurface(0), GetSurface(1), m_sCoarseSetName, m_sCoarseSetName);
	m_pFinalCoarseMap->SetFinalCorrMap(corrs);
}

void PipelineExploreConfmaps::ExploreConformalMaps()
{
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);
	VKString mapEvalName = m_Params.GetStrValue(algorithmName, "MapEvaluator", valid);
	VKString aggregationType = m_Params.GetStrValue(algorithmName, "AggregationMethod", valid);
	bool mobiusVoting = (aggregationType=="VoteInCorrMatrixGreedy");
	assert(GetSurface(0)->GetSurfaceType()=="SurfaceMidEdgeConf" 
		   && GetSurface(1)->GetSurfaceType()=="SurfaceMidEdgeConf");
	SurfaceMidEdgeConf * M1 = (SurfaceMidEdgeConf *) GetSurface(0);
	SurfaceMidEdgeConf * M2 = (SurfaceMidEdgeConf *) GetSurface(1);	
	
	assert(m_pCollectionOfConformalMaps!=NULL);
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("ExploringConformalMaps");							
	if (Verb("Min"))
		std::cout<<"\tExploring Conformal Maps:"<<std::endl;
	
	int progress=0;
	int outOf = 0;
	int everyN = 5;
	while (PickNextNplet(&progress, &outOf))
	{
		if (progress<10)
		{
			if (outOf > 1000000)
				everyN = 100000;
			else if (outOf > 100000)
				everyN = 10000;
		}
		if (Verb("Min") && (progress%everyN==0))
			AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("ExploringConformalMaps", progress, outOf);

		VoteForCurrentConfMap(M1, M2, mapEvalName, mobiusVoting);
	}

	if (aggregationType=="VoteInCorrMatrixGreedy")
	{
		assert(m_pFinalCoarseMap!=NULL);
		double fracToKeep = m_Params.GetDoubleValue(algorithmName, "PickBestConf", valid);
		bool furtherPairs = m_Params.GetStrValue(algorithmName, "FurtherPairsSearch", valid)=="true" && valid;
		m_pFinalCoarseMap->AssignGreedyCorr(fracToKeep, furtherPairs);
	}
	
	if (Verb("Min"))
		std::cout<<std::endl;
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("ExploringConformalMaps");	
}

bool PipelineExploreConfmaps::PickNextNplet(int * curr, int * total)
{
	VKStringList retVal;
	bool valid=false;
	for (int i=0; i<(int)m_CandidateCorrGenerators.size(); i++)
	{
		m_CurrentSequence = &(m_CandidateCorrGenerators[i]->GenerateNextMap(&valid));
		if (valid)
			break;
		else
			m_CurrentSequence = NULL;
	}
	
	if (curr!=NULL || total!=NULL)
	{
		if (curr!=NULL)
			*curr = 0;
		if (total!=NULL)
			*total = 0;
		
		int currVal, currTotal;
		for (int i=0; i<(int)m_CandidateCorrGenerators.size(); i++)
		{
			m_CandidateCorrGenerators[i]->GetProgress(currVal, currTotal);
			if (curr!=NULL)
				*curr += currVal;
			if (total!=NULL)
				*total += currTotal;
		}
	}
	
	return valid;
}

void PipelineExploreConfmaps::VoteForCurrentConfMap(SurfaceMidEdgeConf * M1, 
													SurfaceMidEdgeConf * M2,
													const VKString & confidenceName,
													bool mobiusVoting)
{
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("ExplrVotePrefix");
	if (m_CurrentSequence==NULL)
		return;
	
	// create a conformal map
	MapConformal * newMap=NULL;

//	assert(GetSurface(0)->GetSurfaceType()=="SurfaceMidEdgeConf" 
//		   && GetSurface(1)->GetSurfaceType()=="SurfaceMidEdgeConf");
//	SurfaceMidEdgeConf * M1 = (SurfaceMidEdgeConf *) GetSurface(0);
//	SurfaceMidEdgeConf * M2 = (SurfaceMidEdgeConf *) GetSurface(1);	
	//assert(m_GeneratorSetName1==m_GeneratorSetName2);	
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("ExplrVotePrefix");			
	
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("ExplrCreateConf");	
	newMap = new MapConformal(M1, M2, *m_CurrentSequence, m_GeneratorSetName1);
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("ExplrCreateConf");	
	if (mobiusVoting)
	{
		if (m_pFinalCoarseMap==NULL)
		{
			VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);
			VKString samplesName = m_Params.GetStrValue(algorithmName, "SamplesForErr", valid);
			m_pFinalCoarseMap = new MapCoarse(GetSurface(0), GetSurface(1),
											  M1->GetSampleSet(samplesName),
											  M2->GetSampleSet(samplesName));
			m_pFinalCoarseMap->InitializeVoting();
		}	
//		AnalysisStats::m_GlobalStats.m_Timing.startedProcess("ExplrVoteInMatrix");
		
		double val = m_pFinalCoarseMap->CastVote(newMap, MapCoarse::F2C_MCN_CONFORMAL_EUCLIDEAN);
//		AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("ExplrVoteInMatrix");				
		if (m_pCollectionOfConformalMaps->GetNumMaps()<20)
			m_pCollectionOfConformalMaps->AddMap(newMap, val);
		else
			delete newMap;
	}
	else
	{
		// add a conformal map to our collection
		assert(newMap!=NULL && m_pCollectionOfConformalMaps!=NULL);
		m_pCollectionOfConformalMaps->AddMap(newMap, confidenceName);
		//newMap->SetAdditionalData(GetCurrentVoteDescriptor());
	}
	
	AnalysisStats::m_GlobalStats.Increase("Voting::NumConformalMapsExplored");	
}

VKString PipelineExploreConfmaps::GetCurrentVoteDescriptor()
{
	VKString voteDescr="";
	if (m_CurrentSequence==NULL)
	{
		voteDescr+="Exhausted all votes\n";
	}
	else
	{
		voteDescr+="Map=[";
		for (int i=0; i<(int)m_CurrentSequence->size(); i+=2)
		{
			voteDescr+=VKString::number((*m_CurrentSequence)[i])+" -> ";
			voteDescr+=VKString::number((*m_CurrentSequence)[i+1])+", ";
		}
		voteDescr+=VKString("] Value=")+VKString::number(m_pCollectionOfConformalMaps->GetMapValue(m_pCollectionOfConformalMaps->GetNumMaps()-1))+"\n";
	}
	return voteDescr;
}

// RENDERING
void PipelineExploreConfmaps::Draw(ParamParser & drawParams)
{		
	SurfaceMap * currMap=NULL;
	SurfaceMap * otherMap=NULL;
	
	int drawConfMap = drawParams.GetIntValue("RendererDefault", "SelectedConfMap", valid);
	int distFromConfMap = drawParams.GetIntValue("RendererDefault", "DistFromConfMap", valid);

	if (valid && distFromConfMap!=-1 && m_pCollectionOfConformalMaps!=NULL)
	{
		assert(m_pCollectionOfConformalMaps->GetMapByID(distFromConfMap)->GetSurfaceMapType()=="MapConformal");
		((MapConformal*)m_pCollectionOfConformalMaps->GetMapByID(distFromConfMap))->SetDrawOnlyGenerators(0, 0, 1);
		m_pCollectionOfConformalMaps->GetMapByID(distFromConfMap)->Draw(m_pWindow, &drawParams);
		currMap = m_pCollectionOfConformalMaps->GetMapByID(distFromConfMap);	
		otherMap = currMap;
	}
	
	if (valid && drawConfMap!=-1 && m_pCollectionOfConformalMaps!=NULL)
	{
		assert(m_pCollectionOfConformalMaps->GetMapByID(drawConfMap)->GetSurfaceMapType()=="MapConformal");
		
		if (distFromConfMap!=-1)	// other map exists
			((MapConformal*)m_pCollectionOfConformalMaps->GetMapByID(drawConfMap))->SetDrawOnlyGenerators(1., .2, .2);
		else	// other map does not exist
			((MapConformal*)m_pCollectionOfConformalMaps->GetMapByID(drawConfMap))->SetDrawOnlyGenerators(-1, 0, 0);
		
		m_pCollectionOfConformalMaps->GetMapByID(drawConfMap)->Draw(m_pWindow, &drawParams);
		currMap = m_pCollectionOfConformalMaps->GetMapByID(drawConfMap);
	}
	else
		otherMap = NULL;
	
	VKStringList mapsToDraw = drawParams.GetStrValues("RendererDefault", "RenderMap", valid);
	if (mapsToDraw.contains("BestConformal") && m_pFinalCoarseMap!=NULL)
	{
		m_pFinalCoarseMap->Draw(m_pWindow, &drawParams);
		currMap = m_pFinalCoarseMap;
	}
	
	if (mapsToDraw.contains("BestExtrapolated") && m_pTheFinalMap!=NULL)
	{
		m_pTheFinalMap->Draw(m_pWindow, &drawParams);
		currMap = m_pTheFinalMap;
	}
	
	if (mapsToDraw.contains("Truth") && m_pTruthMap!=NULL)
	{
		m_pTruthMap->Draw(m_pWindow, &drawParams);
		currMap = m_pTruthMap;
	}
	
	if (mapsToDraw.contains("Collection") && m_pCollectionOfMaps!=NULL)
	{
		m_pCollectionOfMaps->Draw(m_pWindow, &drawParams);
		currMap = m_pCollectionOfMaps->GetMapByID(m_pCollectionOfMaps->GetSelectedMapID());
	}
	
	if (mapsToDraw.contains("BenchmarkQuery") && m_pBenchmarkQuery!=NULL)
	{
		m_pBenchmarkQuery->Draw(m_pWindow, &drawParams);	
		currMap = m_pBenchmarkQuery;
	}
	
	VKString drawType = drawParams.GetStrValue("RendererDefault", "DrawMapColorsOnSurf", valid);
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);
	VKString sampleSet = m_Params.GetStrValue(algorithmName, "SamplesMobVoteCast", valid);

	VKString similarityName = "none";
	VKString confidenceName = m_Params.GetStrValue(algorithmName, "MapEvaluator", valid);
	if (drawType=="MapDiscrepancy")
		similarityName = m_Params.GetStrValue(algorithmName, "MapSimilarity", valid);
	
	for (int i=0; i<(int)m_Surfaces.size(); i++)
	{		
		GetSurface(i)->SetDrawSurfaceColors(&drawParams, sampleSet, currMap, otherMap, 
											similarityName, confidenceName);
		GetSurface(i)->Draw(&drawParams);
	}
	
	int numSurfacesToRender = drawParams.GetIntValue("RendererDefault", "NumSurfaces", valid);
	assert(valid);
	if (numSurfacesToRender==1)
	{
		DrawSymmetricGeneratorSet(GetSurface(0), &drawParams);
	}
	else if (numSurfacesToRender==2)
	{
		for (int i=0; i<NumSurfaces(); i++)
			if (i==0)
				DrawGeneratorSet(GetSurface(i), i, 0, &drawParams);
			else
				DrawGeneratorSet(GetSurface(i), i, '\'', &drawParams);		
	}
	else 
		assert(false);
	
	
//	// TEMP STUFF FOR PAPER Dec 30
//	if (m_MapsDueToTopEigenvectors.size()>0)
//	{
//		double rad = GetSurface(1)->GetStandardRadius(&drawParams);
//		int v = GetSurface(0)->GetSelectedVertex(&drawParams);
//		if (v==-1)
//		{
//			std::cout<<"\tVertex Not Selected"<<std::endl;
//			return;
//		}
//		glDisable(GL_LIGHTING);		
//		SurfaceSample sampV(v, GetSurface(0)->GetMesh());
//		R3Point p1 = GetSurface(0)->GetDrawablePosition(sampV, &drawParams);
//		for (int i=0; i<(int)m_MapsDueToTopEigenvectors.size(); i++)
//		{
//			R3Point p2 = GetSurface(1)->GetDrawablePosition(m_MapsDueToTopEigenvectors[i], &drawParams);
//			
//			double err = 1. - m_EigenValues[i] / m_EigenValues[0];
//			glColor3d(1, 1-err, 0);
//			R3Sphere(p2, rad).Draw();
//		}
//	}
}

void PipelineExploreConfmaps::DrawGeneratorSet(SampledSurface * surface, int surfID, char secondChar, 
										ParamParser * params,
										const VKString & renderingParams,
										const VKString & surfaceName)
{	
	if (m_CurrentSequence == NULL)
		return;
	assert(surfID<2 && surfID>=0);
	VKString mySet = params->GetStrValue(renderingParams, "SeqFromSet", valid);
	SurfaceSampleSet * set1 = surface->GetSampleSet(mySet);
	if (!valid || set1==NULL)
	{
		std::cout<<"[WARNING] Cannot find generator set for a sequence. Ignore rendering..."<<std::endl;
		return;
	}
	
	for (int i=0; i<(int)m_CurrentSequence->size()-1; i+=2)
	{
		int id1 = (*m_CurrentSequence)[i+surfID];
		
		assert(id1>=0 && id1<set1->NumSamples() );
		
		R3Point p1 = surface->GetDrawablePosition(set1->GetSample(id1), params, renderingParams, surfaceName);
		
		double r, g, b;
		char buf[3];
		buf[1]=secondChar;
		buf[2]=0;		
		buf[0] = 'A'+(char)i/2;
		r = 1;		g = 0;		b = 0;
		m_pWindow->RenderText(p1, VKString(buf), r, g, b, AnalysisWindow::TEXT_SIZE_LARGE);
	}
}


void PipelineExploreConfmaps::DrawSymmetricGeneratorSet(SampledSurface * surface, ParamParser * params,
												 const VKString & renderingParams,
												 const VKString & surfaceName)
{	
	if (m_CurrentSequence == NULL)
		return;
	VKString mySet = params->GetStrValue(renderingParams, "SeqFromSet", valid);
	SurfaceSampleSet * set1 = surface->GetSampleSet(mySet);
	std::map<int, int> assigned;
	if (!valid || set1==NULL)
	{
		std::cout<<"[WARNING] Cannot find symmetric set for a sequence. Ignore rendering..."<<std::endl;
		return;
	}
	
	int letterID = 0;
	for (int i=0; i<(int)m_CurrentSequence->size()-1; i+=2)
	{
		int id1 = (*m_CurrentSequence)[i];
		int id2 = (*m_CurrentSequence)[i+1];
		
		assert(id1>=0 && id2>=0 && id1<set1->NumSamples() && id2<set1->NumSamples());
		
		//		R3Point p1 = set1->GetSample(id1).GetPosition();
		//		R3Point p2 = set1->GetSample(id2).GetPosition();
		R3Point p1 = surface->GetDrawablePosition(set1->GetSample(id1), params, renderingParams, surfaceName);
		R3Point p2 = surface->GetDrawablePosition(set1->GetSample(id2), params, renderingParams, surfaceName);		
		
		double r, g, b;
		char buf[2];
		buf[1]=0;
		if (id1==id2)	//stationary points
		{
			if (assigned.find(id1)==assigned.end())
			{
				assigned[id1] = letterID;
				buf[0] = 'A'+(char)letterID;
				//m_pWindow->MapIntToColor(letterID, r, g, b);
				r = 1;		g = 0;		b = 0;
				m_pWindow->RenderText(p1, VKString(buf)+"*", r, g, b, AnalysisWindow::TEXT_SIZE_LARGE);
				letterID++;
			}
		}
		else
		{
			if (assigned.find(id1)==assigned.end())
			{
				assigned[id1] = letterID;
				buf[0] = 'A'+(char)letterID;
				//m_pWindow->MapIntToColor(letterID, r, g, b);
				r = 1;		g = 0;		b = 0;
				m_pWindow->RenderText(p1, VKString(buf), r, g, b, AnalysisWindow::TEXT_SIZE_LARGE);			
				letterID++;
			}
			
			if (assigned.find(id2)==assigned.end())
			{
				assigned[id2] = letterID;
				buf[0] = 'A'+(char)letterID;
				//m_pWindow->MapIntToColor(letterID-1, r, g, b);
				r = 1;		g = 0;		b = 0;
				m_pWindow->RenderText(p2, VKString(buf), r, g, b, AnalysisWindow::TEXT_SIZE_LARGE);
				letterID++;
			}
		}
	}
}


