#include "AnalysisPipeline.h"
#include "PipelineExploreConfmaps.h"
#include "FeatureCurvature.h"
#include "PipelineGeneral.h"

VKString AnalysisPipeline::s_Verbous = "Debug";
bool AnalysisPipeline::Verb(const VKString & v)
{
	if (v=="Debug")
		return s_Verbous=="Debug";
	else if (v=="All")
		return s_Verbous=="All"||s_Verbous=="Debug";
	else if (v=="Min")
		return (s_Verbous=="Debug"||s_Verbous=="All"||s_Verbous=="Min");
	else
		return false;
}

AnalysisPipeline::AnalysisPipeline(ParamParser & params)
: m_Params(params)
{
	m_pTheFinalMap = NULL;
	m_pTruthMap = NULL;
	m_pCollectionOfMaps = NULL;
	m_pBenchmarkQuery = NULL;
	m_pWindow = NULL;	
	m_pExternalAlgorithmTime = NULL;
}

AnalysisPipeline::~AnalysisPipeline()
{
}

AnalysisPipeline * AnalysisPipeline::CreatePipeline(ParamParser & params)
{
//	params.PrintParams();
	bool valid;
	VKString algorithmName = params.GetStrValue("Pipeline", "Algorithm", valid);
	assert(valid);
	VKString algorithmType = params.GetStrValue(algorithmName, "Type", valid);
	if (algorithmName=="PipelineGeneral")
		algorithmType = "PipelineGeneral";
	else
		assert(valid);
	
	if (algorithmType=="PipelineExploreConfmaps" || algorithmType=="OneStepSpectral")
		return new PipelineExploreConfmaps(params);
	else if (algorithmType=="PipelineGeneral")
		return new PipelineGeneral(params);
	else
		assert(false);
}

VKString AnalysisPipeline::GetPipelineType()
{
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);
	return m_Params.GetStrValue(algorithmName, "Type", valid);
}

void AnalysisPipeline::TryLoadCoarseMap(const VKString & filename, MapCoarse ** saveTo)
{
	*saveTo = new MapCoarse(GetSurface(0), GetSurface(1));
	if (!((SurfaceMap*)*saveTo)->LoadMap(filename))
	{
		delete *saveTo;
		*saveTo = NULL;
	}
}

void AnalysisPipeline::LoadAndProcessSurfaces(const VKString & algorithmName)
{
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);	
	VKStringList surfaces = m_Params.GetStrValues(algorithmName, "Surface", valid);
	
	std::vector<SampledSurface *> surfaceSymmetries;
	if (loads.contains("Surface") || process.contains("Surface"))	// also might load Distances and Features
	{
		for (int i=0; i<surfaces.count(); i++)
		{
			VKStringList meshes = m_Params.GetStrValues(surfaces[i], "MeshName", valid);
			assert(meshes.count()>0);
			
			for (int j=0; j<meshes.count(); j++)
			{
				surfaceSymmetries.clear();
				m_Surfaces.push_back(CreateSurface(surfaces[i], surfaceSymmetries, j));
				for (int s=0; s<(int)surfaceSymmetries.size(); s++)
					m_Surfaces.push_back(surfaceSymmetries[s]);
			}
		}
	}
	std::cout<<"Loaded "<<m_Surfaces.size()<<" surfaces"<<std::endl;
}

void AnalysisPipeline::LoadSetSamplers(const VKString & algorithmName)
{
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);	
	m_GeneratorSetName1 =  m_Params.GetStrValue(algorithmName, "SamplesMobVoteCast", valid);
	m_GeneratorSetName2 = m_GeneratorSetName1;
	assert(valid);
	
	// some more stuff to load that is not
	VKStringList samplers = m_Params.GetStrValues(algorithmName, "GeneratorSets", valid);	
	
	if (loads.contains("Surface") && loads.contains("Samples") 
		&& loads.contains("Features") && loads.contains("Distances"))
	{		
		if (Verb("Min"))
			std::cout<<"Initializing Generator sets"<<std::endl;	
		AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CreatingGeneratorSets");

		// initialize generator sets
		for (int i=0; i<samplers.count(); i++)
		{
			if (samplers[i]!="none")
				m_CandidateCorrGenerators.push_back(CreateSetSampler(samplers[i], 
																	 m_GeneratorSetName1, m_Surfaces));
		}
		
		AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("CreatingGeneratorSets");				
	}
}

void AnalysisPipeline::LoadTruthMap()
{
	if (m_pTruthMap!=NULL)
		return;
	
	VKString truthMapName = m_Params.GetStrValue("Pipeline", "BenchmarkQuery", valid);		assert(valid);
	VKString mapType = m_Params.GetStrValue(truthMapName, "MapType", valid);				 assert(valid);
	VKStringList truthPaths = m_Params.GetStrValues(truthMapName, "Path", valid);			 assert(valid);
	if (truthPaths.count()==1)
		truthPaths.push_back(truthPaths[0]);
	VKStringList truthPnts = m_Params.GetStrValues(truthMapName, "FileVertexIDs", valid);	 assert(valid);
	VKString truthCorr = m_Params.GetStrValue(truthMapName, "FileFeaturesMap", valid);	
	
	m_pTruthMap = new MapCoarse(GetSurface(0), GetSurface(1));
	
	if (mapType=="CoarseMap_FeaturePointsSymmetry")
		m_pTruthMap->LoadCorrFeaturePntSymmetryMap(truthPaths[0]+truthPnts[0], truthPaths[0]+truthCorr);
	else if (mapType=="CoarseMap_FeaturePointsCorrespondence")
		m_pTruthMap->LoadCorrSameFeaturePntIDs(truthPaths[0]+truthPnts[0], truthPaths[1]+truthPnts[1]);	
	else if (mapType=="CoarseMap_SameVertexIDs")
		m_pTruthMap->LoadCorrSameVertexIDs();
	else if (mapType=="CoarseMap_FunkCor")
		m_pTruthMap->LoadMapInVertexByVertexFormat(truthPaths[0]+truthPnts[0]);
	else
		assert(false);
	
	truthMapName = m_Params.GetStrValue("Pipeline", "SecondBenchmark", valid);	
	if (valid && truthMapName!="none")
	{
		mapType = m_Params.GetStrValue(truthMapName, "MapType", valid);			assert(valid);
		truthPaths = m_Params.GetStrValues(truthMapName, "Path", valid);		assert(valid);
		if (truthPaths.count()==1)
			truthPaths.push_back(truthPaths[0]);
		truthPnts = m_Params.GetStrValues(truthMapName, "FileVertexIDs", valid);	 
		assert(valid);
		truthCorr = m_Params.GetStrValue(truthMapName, "FileFeaturesMap", valid);	
		assert(valid);
		
		m_pAlternativeTruth = new MapCoarse(GetSurface(0), GetSurface(1));
		m_pAlternativeSymmetricTruth = new MapCoarse(GetSurface(0), GetSurface(1));
		
		if (mapType=="CoarseMap_CreateTrueAndSymmetric")
		{
			// ground truth
			m_pAlternativeTruth->LoadCorrSameFeaturePntIDs(truthPaths[0]+truthPnts[0], 
														   truthPaths[1]+truthPnts[1]);	
			m_pAlternativeSymmetricTruth->LoadCorrSameFeaturePntIDsFlipSymmetry(truthPaths[0]+truthPnts[0], 
																				truthPaths[1]+truthPnts[1],
																				truthPaths[0]+truthCorr);
		}
		else
			assert(false);
	}	
}

void AnalysisPipeline::AddBenchmark(const VKString & confidenceName, SurfaceMap * predictionMap, 
									ParamParser * benchParams, int statID, double score,
									bool includesSymFlip)
{
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);
	assert(valid);
	VKString genName = m_Params.GetStrValue(algorithmName, "SamplesMobVoteCast", valid);
	
	// create error evaluator	
	bool flipped = false;
	SurfaceMapConfidence * confidence = NULL;
	if (!includesSymFlip)
		confidence = predictionMap->GetMapConfidenceCalculator(confidenceName);
	else 
	{
		SurfaceMapConfidence * confidenceNotFlipped;
		SurfaceMapConfidence * confidenceFlipped;		
		confidenceNotFlipped = predictionMap->GetMapConfidenceCalculator(confidenceName+"_NotFlipped");
		confidenceFlipped = predictionMap->GetMapConfidenceCalculator(confidenceName+"_Flipped");
		if (confidenceNotFlipped->Error() > confidenceFlipped->Error())
		{
			confidence = confidenceFlipped;
			flipped = true;
		}
		else
		{
			confidence = confidenceNotFlipped;
			flipped = false;
		}
	}
	
	std::vector<double> perSampleErrors;
	confidence->GetPerSampleErrors(GetSurface(1), perSampleErrors, false);
	double corrRate = confidence->Confidence();
	double aveGeo = confidence->Error();
	
	if (Verb("Min"))
	{
		std::cout<<"\t -> Error: "<<aveGeo<<std::endl;
		std::cout<<"\t -> Matched: "<<corrRate<<std::endl;
		std::cout<<"\t -> ObjFn: "<<score<<std::endl;		
		if (includesSymFlip)
			std::cout<<"\t -> "<<(flipped ? ("Flipped!") : ("Not Flipped"))<<std::endl;
	}
		
	VKString statsStr=((statID==0) ? "Stats" : (VKString("Stats")+VKString::number(statID)));

	if (includesSymFlip)
	{
		statsStr=((statID==0) ? "StatsWithFlip" : (VKString("StatsWithFlip")+VKString::number(statID)));
		if (flipped)
			(*benchParams).m_ParamModules[statsStr]["Flipped"].push_back("true");
		else
			(*benchParams).m_ParamModules[statsStr]["Flipped"].push_back("false");
	}
		
	// print error stats
	(*benchParams).m_ParamModules[statsStr]["Average"].push_back(VKString::number(aveGeo));
	(*benchParams).m_ParamModules[statsStr]["Guessed"].push_back(VKString::number(corrRate));	
	(*benchParams).m_ParamModules[statsStr]["ObjFn"].push_back(VKString::number(score));
	SurfaceSampleSet * set1 = GetSurface(0)->GetSampleSet(genName);
	SurfaceSampleSet * set2 = GetSurface(1)->GetSampleSet(genName);	
	if (set1!=NULL && set2!=NULL)
	{
		int G1 = set1->NumSamples();
		int G2 = set2->NumSamples();		
		(*benchParams).m_ParamModules[statsStr]["NumGen"].push_back(VKString::number(G1));
		(*benchParams).m_ParamModules[statsStr]["NumGen"].push_back(VKString::number(G2));
	}
	
	for (int i=0; i<(int)perSampleErrors.size(); i++)
	{
		VKString errAtSampStr = VKString::number(perSampleErrors[i]);
		(*benchParams).m_ParamModules[statsStr]["PerSampleErrors"].push_back(errAtSampStr);
	}
}

void AnalysisPipeline::LoadSamplesFromCoarseAlgorithm(const VKString & otherMapPrefix)
{
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);
	VKString surfDescr = m_Params.GetStrValue(algorithmName, "Surface", valid);	
	VKString confMapExt = m_Params.GetStrValue(algorithmName, "ConfMapExt", valid);	
	VKStringList meshNames = m_Params.GetStrValues(surfDescr, "MeshName", valid);

	assert(valid && meshNames.count()>0);
	if (meshNames.count()<2)
		meshNames.push_back(meshNames[0]);

	VKString filename=otherMapPrefix+"_Map_"+meshNames[0]+"_to_"+meshNames[1];
	VKString coarseMapName = LocalFilename(meshNames[0], filename, confMapExt);

	SurfaceSampleSet * fromSet=new SurfaceSampleSet();
	SurfaceSampleSet * toSet=new SurfaceSampleSet();	
	SurfaceMap * coarseMap = SurfaceMap::CreateMap(coarseMapName, GetSurface(0), GetSurface(1));
	const SurfaceSampleSet * validDomain = coarseMap->GetValidDomain();
	for (int i=0; i<validDomain->NumSamples(); i++)
	{
		const SurfaceSample & s1 = validDomain->GetSample(i);
		SurfaceSample s2;
		if (!s1.Invalid())
		{
			s2 = coarseMap->ForwardMap(s1);
			if (!s2.Invalid())
			{
				fromSet->AddSample(s1);
				toSet->AddSample(s2);				
			}
		}
	}
	
	GetSurface(0)->AddSampleSet("FinalConformalMap", fromSet);
	GetSurface(1)->AddSampleSet("FinalConformalMap", toSet);	
}

void AnalysisPipeline::LoadAndProcessBenchmark()
{
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("AABasicBenchmarking");			
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);
	VKStringList processes = m_Params.GetStrValues("Pipeline", "Process", valid);
	if (!loads.contains("BenchmarkQuery") && !processes.contains("BenchmarkQuery"))
		return;

	if (Verb("Min"))
		std::cout<<"Benchmarking"<<std::flush;
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("Benchmarking");	
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);	
	VKString algNameOut = m_Params.GetStrValue(algorithmName, "OutputName", valid);
	if (!valid)
		algNameOut = algorithmName;
	
	VKString surfDescr = m_Params.GetStrValue(algorithmName, "Surface", valid);
	assert(valid);
	VKStringList meshNames = m_Params.GetStrValues(surfDescr, "MeshName", valid);
	assert(valid && meshNames.count()>0);
	
	ParamParser * benchmarkParams=NULL;
	
	VKString coarseMapFilename;
	if (meshNames.count()==1)
		coarseMapFilename = algNameOut+"_Map_"+meshNames[0]+"_to_"+meshNames[0];
	else if (meshNames.count()==2)
		coarseMapFilename = algNameOut+"_Map_"+meshNames[0]+"_to_"+meshNames[1];	
	else
		assert(false);	
	
	VKString workFolder = m_Params.GetStrValue("Pipeline", "WorkFolder", valid);
	VKString filenameBenchParams = workFolder+"/BenchResult_"+algNameOut+"_"+meshNames[0]+"_to_";
	if (meshNames.count()==1)
		filenameBenchParams += meshNames[0]+".txt";
	else
		filenameBenchParams += meshNames[1]+".txt";			

	m_pBenchmarkQuery = NULL;
	if (!processes.contains("BenchmarkQuery"))
	{
		std::cout<<"Loading benchmark"<<std::endl;
		std::ifstream textStream(filenameBenchParams.c_str());
		if (textStream.is_open())
		{
			textStream.close();
			benchmarkParams = new ParamParser(filenameBenchParams);
		}
		
		SurfaceMap * queryMap = SurfaceMap::CreateMap(coarseMapFilename, 
													  GetSurface(0), GetSurface(1));
		if (queryMap!=NULL)
		{
			assert(queryMap!=NULL && queryMap->GetSurfaceMapType()=="MapCoarse");
			m_pBenchmarkQuery = (MapCoarse*)queryMap;
		}
	}

	VKString confidenceName = m_Params.GetStrValue("Pipeline", "BenchmarkQuery", valid);
	assert(valid);
	VKString alternativeConfidenceName = m_Params.GetStrValue("Pipeline", "SecondBenchmark", valid);	
	if (!valid)
		alternativeConfidenceName = "none";
	VKString outCorr = m_Params.GetStrValue(confidenceName, "Extension", valid);
	assert(valid);			
	LoadTruthMap();
	
	assert(m_pTheFinalMap!=NULL);
	assert(m_pTruthMap!=NULL);
	
	if (m_pBenchmarkQuery==NULL)
	{
		VKString atSamplesRange = m_Params.GetStrValue(confidenceName, "AtSamplesRange", valid);
		VKString atSamples = m_Params.GetStrValue(confidenceName, "AtSamples", valid);	
		VKString otherAlgName = m_Params.GetStrValue(confidenceName, "RangeAlgPrefix", valid);
		if (atSamplesRange=="FinalConformalMap" || atSamples=="FinalConformalMap")
			LoadSamplesFromCoarseAlgorithm(otherAlgName);

		SurfaceSampleSet * domainSet = new SurfaceSampleSet();		
		if (atSamplesRange!="none")
		{
			SurfaceSampleSet * sampleSet = GetSurface(1)->GetSampleSet(atSamplesRange);
			assert(sampleSet!=NULL);

			for (int i=0; i<sampleSet->NumSamples(); i++)
			{
				// verify that samples at which truth is to be evaluated exist / known
				assert(!sampleSet->GetSample(i).Invalid());	
				const SurfaceSample & invSamp = m_pTruthMap->InverseMap(sampleSet->GetSample(i));
				assert(!invSamp.Invalid());	
				domainSet->AddSample(invSamp);
			}
		}
		if (atSamples!="none")
		{
			SurfaceSampleSet * sampleSet = GetSurface(0)->GetSampleSet(atSamples);
			assert(sampleSet!=NULL);
			for (int i=0; i<sampleSet->NumSamples(); i++)			
				domainSet->AddSample(sampleSet->GetSample(i));
		}		
		
		assert(domainSet!=NULL && domainSet->NumSamples()>0);

		if (Verb("Min"))
			std::cout<<" "<<domainSet->NumSamples()<<" Query Points"<<std::endl;
			
		m_pBenchmarkQuery = new MapCoarse(domainSet, NULL, 
										  m_pTheFinalMap, MapCoarse::F2C_FORWARD_ONLY);
		((SurfaceMap*)m_pBenchmarkQuery)->SaveMap(LocalFilename(meshNames[0], 
																coarseMapFilename, outCorr));		
	}
	m_pBenchmarkQuery->SetCompareToTruthRendering(confidenceName, m_pTruthMap);	
	
	if (benchmarkParams==NULL)
	{		
		// General log stuff
		benchmarkParams = new ParamParser();
		
		int tryTopClusters=ClustersToAnalyze();
		
		if (Verb("Min"))
			std::cout<<"\tClusters to Analyze: "<<tryTopClusters<<std::endl;
			
		// NOTE: assume collection is sorted
		double score=0;
		if (tryTopClusters==1)
		{
			std::cout<<"Analyzing Final Cluster: "<<std::endl;			
			std::cout<<"\tGround Truth"<<std::endl;
			AddBenchmark(confidenceName, m_pTheFinalMap, benchmarkParams, 0, 0, false);
			if (alternativeConfidenceName!="none")
			{
				std::cout<<"\tAllowing Flips"<<std::endl;
				AddBenchmark(alternativeConfidenceName, m_pTheFinalMap, benchmarkParams, 0, 0, true);
			}
		}
		else
		{
			for (int i=0; i<tryTopClusters; i++)
			{
				std::cout<<"Analyzing Cluster: "<<i<<std::endl;
				SurfaceMap * currMap = m_pCollectionOfMaps->GetMapByID(i, &score);
				if (i==0)
					currMap = m_pTheFinalMap;
				std::cout<<"\tGround Truth"<<std::endl;
				AddBenchmark(confidenceName, currMap, benchmarkParams, i, score, false);
				if (alternativeConfidenceName!="none")
				{
					std::cout<<"\tAllowing Flips"<<std::endl;
					AddBenchmark(alternativeConfidenceName, currMap, benchmarkParams, i, score, true);
				}
			}
		}

		AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("AABasicCreateMaps");			
		TimeProfiler & timing = AnalysisStats::m_GlobalStats.m_Timing;
		int timeVal = timing.getTotalProcessLengthInSecondsInt("SurfaceCorrespondenceConfMap");
		int preproc = timing.getTotalProcessLengthInSecondsInt("AABasicCreateSurfaces");
		int findMaps = timing.getTotalProcessLengthInSecondsInt("AABasicCreateMaps");
		int findMapsSrch = timing.getTotalProcessLengthInSecondsInt("AABasicCreateMapsSearch");
		int findMapsCluster = timing.getTotalProcessLengthInSecondsInt("AABasicCreateMapsCluster");
		int findMapsWeights = timing.getTotalProcessLengthInSecondsInt("AABasicCreateMapsWeights");		
		int benchmark = timing.getTotalProcessLengthInSecondsInt("AABasicBenchmarking");
		//		int addProc = timing.getTotalProcessLengthInSecondsInt("AABasicAdditionalOutput");		
		(*benchmarkParams).m_ParamModules["Stats"]["Folder"].push_back(meshNames[0]);
		if (m_pExternalAlgorithmTime==NULL)
		{
			(*benchmarkParams).m_ParamModules["Stats"]["Time"].push_back(VKString::number(timeVal));
			(*benchmarkParams).m_ParamModules["Stats"]["TimeInit"].push_back(VKString::number(preproc));
			(*benchmarkParams).m_ParamModules["Stats"]["TimeMaps"].push_back(VKString::number(findMaps));
			(*benchmarkParams).m_ParamModules["Stats"]["TimeMapsSearch"].push_back(VKString::number(findMapsSrch));	
			(*benchmarkParams).m_ParamModules["Stats"]["TimeMapsCluster"].push_back(VKString::number(findMapsCluster));	
			(*benchmarkParams).m_ParamModules["Stats"]["TimeMapsWeights"].push_back(VKString::number(findMapsWeights));
		}
		else	
		{
			VKString initTime = m_pExternalAlgorithmTime->GetStrValue("ExternalTiming", "InitTime", valid);
			if (!valid)	initTime = VKString::number(0);
			VKString findMaps = m_pExternalAlgorithmTime->GetStrValue("ExternalTiming", "FindMaps", valid);
			if (!valid)	findMaps = VKString::number(0);			
			bool ok1, ok2;
			(*benchmarkParams).m_ParamModules["Stats"]["TimeInit"].push_back(initTime);			
			(*benchmarkParams).m_ParamModules["Stats"]["TimeMaps"].push_back(findMaps);
			int totalTime = initTime.toInt(&ok1) + findMaps.toInt(&ok2) + benchmark;
			assert(ok1 && ok2);
			(*benchmarkParams).m_ParamModules["Stats"]["Time"].push_back(VKString::number(totalTime));
		}
		(*benchmarkParams).m_ParamModules["Stats"]["TimeBench"].push_back(VKString::number(benchmark));
		//		(*benchmarkParams).m_ParamModules["Stats"]["TimeAdd"].push_back(VKString::number(addProc));
		
		benchmarkParams->WriteParams(filenameBenchParams);
	}
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("Benchmarking");	
}

void AnalysisPipeline::WriteLog()
{
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);
	assert(valid);
	VKString outputDir = m_Params.GetStrValue(algorithmName, "OutputName", valid);	
	assert(valid);
	VKString surfDescr = m_Params.GetStrValue(algorithmName, "Surface", valid);
	assert(valid);
	VKStringList meshNames = m_Params.GetStrValues(surfDescr, "MeshName", valid);
	assert(valid && meshNames.count()>0);
	VKString mesh1 = meshNames[0];
	VKString mesh2 = meshNames[(meshNames.count()==2) ? 1 : 0];	
	VKString logName = LocalFilename(".log.txt");	
	VKString timingName = LocalFilename(".time.txt");			
	
	std::ofstream textStream(logName.c_str(), std::ios::out);
	assert(textStream.is_open());
	textStream<<"[WARNING] Log is empty..."<<std::endl;
	textStream.close();	

	if (m_Params.GetStrValue("Pipeline", "RecordTiming", valid)=="true" && valid)
	{
		std::ofstream textStreamTime(timingName.c_str(), std::ios::out);
		assert(textStreamTime.is_open());
		textStreamTime<<AnalysisStats::m_GlobalStats.m_Timing.profile().c_str()<<std::endl;
		textStreamTime.close();
	}
}

int AnalysisPipeline::ClustersToAnalyze()
{
	double fracToAnalyze = m_Params.GetDoubleValue("Pipeline", "AnalyzeCollFrac", valid);
	if (!valid)
		fracToAnalyze=0;
	int numToAnalyze = m_Params.GetIntValue("Pipeline", "AnalyzeCollNum", valid);	
	if (!valid)
		numToAnalyze = -1;
	std::cout<<"fracToAnalyze="<<fracToAnalyze<<" num="<<numToAnalyze<<std::endl;
	
	if (numToAnalyze<0)
	{
		if (m_pCollectionOfMaps!=NULL && m_pCollectionOfMaps->GetNumMaps()!=0)
			return m_pCollectionOfMaps->GetNumMaps();
		else
			return 1;			
	}
	if (m_pCollectionOfMaps==NULL || m_pCollectionOfMaps->GetNumMaps()==0)
		return 1;
	m_pCollectionOfMaps->SortMaps();
	int id;		m_pCollectionOfMaps->GetBestMap(&id);
	double threshold = m_pCollectionOfMaps->GetMapValue(id) * fracToAnalyze;
	for (int i=0; i<m_pCollectionOfMaps->GetNumMaps(); i++)
		if (m_pCollectionOfMaps->GetMapValue(i) < threshold)
			return vkMin(i, numToAnalyze);
	return vkMin(numToAnalyze, m_pCollectionOfMaps->GetNumMaps());
}

void AnalysisPipeline::FillSymmetricSurfacesConf(VKStringList & typeSymmetries, 
												 R3Mesh * mesh, R3Mesh * flatMesh,
												 std::vector<SampledSurface*> & fillMe)
{
	for (int i=0; i<typeSymmetries.count(); i++)
	{
		if (typeSymmetries[i]=="Reflection")
			fillMe.push_back(new SurfaceMidEdgeConf(CopyMesh(mesh), flatMesh, true));
		else if (typeSymmetries[i]=="Identity")
			fillMe.push_back(new SurfaceMidEdgeConf(CopyMesh(mesh), flatMesh, false));
		else if (typeSymmetries[i]!="none")
			assert(false);
	}
}

R3Mesh * AnalysisPipeline::CopyMesh(R3Mesh * mesh)
{
	R3Mesh * newMesh = new R3Mesh();
	for (int i=0; i<mesh->NVertices(); i++)
		newMesh->CreateVertex(mesh->VertexPosition(mesh->Vertex(i)));
	for (int i=0; i<mesh->NFaces(); i++)	
	{
		R3MeshFace * face = mesh->Face(i);
		int v1 = mesh->VertexID(mesh->VertexOnFace(face, 0));
		int v2 = mesh->VertexID(mesh->VertexOnFace(face, 1));
		int v3 = mesh->VertexID(mesh->VertexOnFace(face, 2));		
		newMesh->CreateFace(newMesh->Vertex(v1), newMesh->Vertex(v2), newMesh->Vertex(v3));
	}
	return newMesh;
}

SampledSurface * AnalysisPipeline::CreateSurface(const VKString & surfaceDescriptor, 
												 std::vector<SampledSurface*> &symmetricSurfaces, 
												 int meshID)
{
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CreatingSurface");	
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);
	
	VKStringList typesOfSymmetries = m_Params.GetStrValues(surfaceDescriptor, "LoadSymmetry", valid);
	
	VKString surfaceType = m_Params.GetStrValue(surfaceDescriptor, "Type", valid);
	assert(valid);
	SampledSurface * surface=NULL;
	SurfaceMidEdgeConf * needFlattening=NULL;
	
	VKStringList meshNames = m_Params.GetStrValues(surfaceDescriptor, "MeshName", valid);
	assert(valid);
	VKString meshName;

	if (meshNames.count()==1 || meshID<0)
		meshName = meshNames[0];
	else if (meshNames.count() > meshID)
		meshName = meshNames[meshID];
	else
		assert(false);
		
	VKString meshExt = m_Params.GetStrValue(surfaceDescriptor, "MeshExt", valid);	
	assert(valid);
	VKString flatExt = m_Params.GetStrValue(surfaceDescriptor, "FlatExt", valid);
	assert(valid);
	
	R3Mesh * mesh = new R3Mesh();
	mesh->ReadFile(LocalFilename(meshName, meshName, meshExt).c_str());
	assert(mesh!=NULL);
	if (Verb("Min"))
		std::cout<<"Loading Surface: "<<LocalFilename(meshName, meshName, meshExt).c_str()<<std::endl;
	if (surfaceType=="SurfaceMidEdgeConf")
	{
		VKString flatOff = LocalFilename(meshName, flatExt);
		FILE * checkExist = fopen(flatOff.c_str(), "r");
		if (checkExist==NULL || process.contains("Surface"))
		{
			needFlattening = new SurfaceMidEdgeConf(mesh, false);
			surface = needFlattening;
		}
		else
		{
			R3Mesh * flatMesh = new R3Mesh();
			flatMesh->ReadFile(flatOff.c_str());
			surface = new SurfaceMidEdgeConf(mesh, flatMesh, false);
			delete flatMesh;
			flatMesh = new R3Mesh();
			flatMesh->ReadFile(flatOff.c_str());
			FillSymmetricSurfacesConf(typesOfSymmetries, mesh, flatMesh, symmetricSurfaces);
			delete flatMesh;
		}
	}
	else if (surfaceType=="Surface")
	{
		surface = new SampledSurface(mesh);
		for (int i=0; i<(int)typesOfSymmetries.count(); i++)
			if (typesOfSymmetries[i]!="none")
				symmetricSurfaces.push_back(new SampledSurface(mesh));
	}
	assert(surface!=NULL);
	surface->SetName(meshName);
	
	// Load Camera
	VKStringList cameraFiles = m_Params.GetStrValues(surfaceDescriptor, "CameraFile", valid);
	VKStringList cameraPositions = m_Params.GetStrValues(surfaceDescriptor, "CameraOrigin", valid);
	VKStringList cameraOrientations = m_Params.GetStrValues(surfaceDescriptor, "CameraTowards", valid);
	VKStringList cameraUps = m_Params.GetStrValues(surfaceDescriptor, "CameraUp", valid);	
	VKStringList cameraLines = m_Params.GetStrValues(surfaceDescriptor, "CameraLine", valid);		
		
	AddSurfaceCamera(cameraFiles, cameraPositions, cameraOrientations, cameraUps, cameraLines, meshID);

	// add camera due symmetric surfaces
	for (int i=0; i<typesOfSymmetries.count(); i++)
	{
		if (typesOfSymmetries[i]=="none")
			continue;
		int l = (int)m_InitCameraPositions.size();
		if (m_InitCameraPositions[l-1]==NULL)
			m_InitCameraPositions.push_back(NULL);
		else
		{
			m_InitCameraPositions.push_back(new double[9]);
			for (int i=0; i<9; i++)
				m_InitCameraPositions[l][i] = m_InitCameraPositions[l-1][i];
		}
	}	
	
	
	if (loads[1]=="Distances")
	{
		LoadAndProcessDistances(surface, surfaceDescriptor);
		LoadAndProcessFeatures(surface, surfaceDescriptor, &needFlattening, symmetricSurfaces, meshName);
		LoadAndProcessSampleSets(surface, surfaceDescriptor);
	}
	else if (loads[3]=="Distances")
	{
		LoadAndProcessFeatures(surface, surfaceDescriptor, &needFlattening, symmetricSurfaces, meshName);
		LoadAndProcessSampleSets(surface, surfaceDescriptor);
		LoadAndProcessDistances(surface, surfaceDescriptor);
	}
	else 
	{
		assert(!loads.contains("Distances"));
		LoadAndProcessFeatures(surface, surfaceDescriptor, &needFlattening, symmetricSurfaces, meshName);
		LoadAndProcessSampleSets(surface, surfaceDescriptor);
	}
	
	assert(needFlattening==NULL);
	
	if (symmetricSurfaces.size()!=0)
	{
		assert((int)symmetricSurfaces.size()==typesOfSymmetries.count());
		for (int i=0; i<(int)symmetricSurfaces.size(); i++)
			symmetricSurfaces[i]->CopyFromSameSurface(surface);
	}
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("CreatingSurface");	
	return surface;
}

void AnalysisPipeline::AddSurfaceCamera(const VKStringList & cameraFiles, 
										const VKStringList & cameraPositions,
										const VKStringList & cameraOrientations,
										const VKStringList & cameraUps,
										const VKStringList & cameraLines,
										int meshID)
{
	VKString cameraFile = "none";
	int posID = -1;
	int orientID = -1;
	int upID = -1;
	int cameraLine = -1;
	
	if (cameraFiles.count()==1 || meshID<0)
		cameraFile = ((VKStringList &)cameraFiles)[0];
	else if (cameraFiles.count() > meshID)
		cameraFile = ((VKStringList &)cameraFiles)[meshID];
	
	if (cameraPositions.count()==3 || meshID<0)
		posID = (((VKStringList &)cameraPositions)[0]!="none") ? 0 : -1;
	else if (cameraPositions.count() > meshID*3+2)
		posID = (((VKStringList &)cameraPositions)[meshID*3]!="none") ? (meshID*3) : (-1);
	
	if (cameraOrientations.count()==3 || meshID<0)
		orientID = (((VKStringList &)cameraOrientations)[0]!="none") ? 0 : (-1);
	else if (cameraOrientations.count() > meshID*3+2)
		orientID = (((VKStringList &)cameraOrientations)[meshID*3]!="none") ? (meshID*3) : (-1);
	
	if (cameraUps.count()==3 || meshID<0)
		upID = (((VKStringList &)cameraUps)[0]!="none") ? 0 : (-1);
	else if (cameraUps.count() > meshID*3+2)
		upID = (((VKStringList &)cameraUps)[meshID*3]!="none") ? (meshID*3) : (-1);
	
	if(cameraFile!="none")
	{
		if (cameraLines.count()==1 || meshID<0)
			cameraLine = ((VKStringList &)cameraLines)[0].toInt();
		else if (cameraLines.count() > meshID)
			cameraLine = ((VKStringList &)cameraLines)[meshID].toInt();
		else
			cameraLine = -1;
	}
	
	if (posID!=-1 && orientID!=-1 && upID!=-1)	// has camera pos in settings file
	{
		m_InitCameraPositions.push_back(new double[9]);
		int l = m_InitCameraPositions.size()-1;
		
		m_InitCameraPositions[l][0] = ((VKStringList &)cameraPositions)[posID+0].toDouble();
		m_InitCameraPositions[l][1] = ((VKStringList &)cameraPositions)[posID+1].toDouble();
		m_InitCameraPositions[l][2] = ((VKStringList &)cameraPositions)[posID+2].toDouble();
		m_InitCameraPositions[l][3] = ((VKStringList &)cameraOrientations)[orientID+0].toDouble();
		m_InitCameraPositions[l][4] = ((VKStringList &)cameraOrientations)[orientID+1].toDouble();
		m_InitCameraPositions[l][5] = ((VKStringList &)cameraOrientations)[orientID+2].toDouble();
		m_InitCameraPositions[l][6] = ((VKStringList &)cameraUps)[upID+0].toDouble();
		m_InitCameraPositions[l][7] = ((VKStringList &)cameraUps)[upID+1].toDouble();
		m_InitCameraPositions[l][8] = ((VKStringList &)cameraUps)[upID+2].toDouble();
	}
	else if (cameraLine>=0)	// try loading camera from a camera file
	{
		ifstream dataStream(cameraFile.c_str(), ios::in);
		if(!dataStream.is_open())
		{
			std::cout<<"[WARNING] Could not open camera file: "<<cameraFile.c_str()<<std::endl;
			m_InitCameraPositions.push_back(NULL);
		}
		else
		{
			VKString tempStr;
			for (int i=0; i<=cameraLine; i++)
				tempStr.readLine(dataStream);
			VKStringList cameraParams = tempStr.split(" ");
			assert(cameraParams.count()>=9);
			m_InitCameraPositions.push_back(new double[9]);		
			for (int i=0; i<9; i++)
			{
				m_InitCameraPositions[m_InitCameraPositions.size()-1][i] = ((VKStringList &)cameraParams)[i].toDouble(&valid);
				assert(valid);
			}
		}
	}
	else
	{
		std::cout<<"[WARNING] Camera not specified"<<std::endl;
		m_InitCameraPositions.push_back(NULL);	
	}
}

void AnalysisPipeline::LoadAndProcessDistances(SampledSurface * surface, const VKString & surfaceDescriptor)
{
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);
	
	if (loads.contains("Distances") || process.contains("Distances"))
	{
		AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CreatingSurface_Distances");		
		if (Verb("Min"))
			std::cout<<"Loading Distances"<<std::endl;	
		
		assert(surface!=NULL);
		VKString defaultMetric = m_Params.GetStrValue(surfaceDescriptor, "DefaultMetric", valid);
		assert(valid);
		SurfaceDistance * distance = CreateDistanceMetric(defaultMetric, surface);
		assert(distance!=NULL);
		surface->AddDistanceMetric("default", distance);
		surface->AddDistanceMetric(defaultMetric, distance);	
		
		VKStringList additionalMetrics = m_Params.GetStrValues(surfaceDescriptor, "AdditionalMetrics", valid);
		for (int i=0; i<additionalMetrics.count(); i++)
		{
			SurfaceDistance * metric = CreateDistanceMetric(additionalMetrics[i], surface);
			assert(metric!=NULL);
			surface->AddDistanceMetric(additionalMetrics[i], metric);
		}
		
		AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("CreatingSurface_Distances");
	}
}

void AnalysisPipeline::LoadAndProcessFeatures(SampledSurface * surface, const VKString & surfaceDescriptor,
											  SurfaceMidEdgeConf ** needFlattening, 
											  std::vector<SampledSurface*> &symmetricSurfaces,
											  const VKString & meshName )	
{
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);
	VKStringList typesOfSymmetries = m_Params.GetStrValues(surfaceDescriptor, "LoadSymmetry", valid);
	
	if (loads.contains("Features") || process.contains("Features"))
	{
		if (Verb("Min"))
			std::cout<<"Loading Features"<<std::endl;
		
		assert(surface!=NULL);
		VKStringList features = m_Params.GetStrValues(surfaceDescriptor, "Features", valid);

		for (int i=0; i<features.count(); i++)
		{
			if (features[i]=="none")
				continue;
			AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CreatingSurface_Features");
			SurfaceFeature * feature = CreateFeature(features[i], surface);
			assert(feature!=NULL);
			surface->AddFeature(features[i], feature);
			
			// NOTE: to flatten mesh more robustly we use some geodesic feature points 
			//			(to cut and make into canonical frame). They are calculated on the fly.
			VKString featureType = m_Params.GetStrValue(features[i], "Type", valid);
			AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("CreatingSurface_Features");
			
			if ((featureType=="FeatureAGD" || featureType=="FeatureFastAGD")  && (*needFlattening)!=NULL)
			{
				VKString flatExt = m_Params.GetStrValue(surfaceDescriptor, "FlatExt", valid);
				int numSmoothIter = m_Params.GetIntValue(surfaceDescriptor, "FlatSmoothIters", valid);				
				VKString flatOff = LocalFilename(meshName, flatExt);
				AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CreatingSurface_Flattening");
				R3Mesh * flatMesh = (*needFlattening)->Flatten(numSmoothIter, feature);
				AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("CreatingSurface_Flattening");										
				flatMesh->WriteFile(flatOff.c_str());
				*needFlattening = NULL;
				FillSymmetricSurfacesConf(typesOfSymmetries, surface->GetMesh(), flatMesh, symmetricSurfaces);
			}
		}
	}	
}

void AnalysisPipeline::LoadAndProcessSampleSets(SampledSurface * surface, const VKString & surfaceDescriptor)
{
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);
	
	if (loads.contains("Samples") || process.contains("Samples"))
	{
		AnalysisStats::m_GlobalStats.m_Timing.startedProcess("CreatingSurface_Samples");		
		if (Verb("Min"))
			std::cout<<"Loading Samples"<<std::endl;	
		
		assert(surface!=NULL);
		VKStringList sampleSets = m_Params.GetStrValues(surfaceDescriptor, "Samples", valid);
		assert(valid);
		for (int i=0; i<sampleSets.count(); i++)
		{
			SurfaceSampleSet* sampleSet = CreateSampleSet(sampleSets[i], surface);
			assert (sampleSet!=NULL);
			surface->AddSampleSet(sampleSets[i], sampleSet);
			std::cout<<"\tSample Set: "<<sampleSets[i].c_str()<<" Num="<<sampleSet->NumSamples()<<std::endl;
		}		
		AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("CreatingSurface_Samples");				
	}	
}

SurfaceFeature * AnalysisPipeline::CreateFeature(const VKString & featureDescriptor, 
												 SampledSurface * surface)
{
	VKString meshName = surface->GetName();
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);	
	VKString featureType = m_Params.GetStrValue(featureDescriptor, "Type", valid);
	VKString featureExt = m_Params.GetStrValue(featureDescriptor, "Extension", valid);
	VKString featureFile = LocalFilename(meshName, featureExt);
	
	assert(valid);
	SurfaceFeature * feature = NULL;
	int nInit = m_Params.GetIntValue(featureDescriptor, "NInit", valid);
	if (!valid)
		nInit = 256;
	if (featureType=="FeatureAGD")
		feature = new FeatureAGD(surface, "default", true, featureDescriptor);
	else if (featureType=="FeatureFastAGD")
		feature = new FeatureFastAGD(surface, nInit, featureDescriptor);
	else if (featureType=="FeatureCurvature")
		feature = new FeatureCurvature(surface, nInit, featureDescriptor);
	
	assert(feature!=NULL);
	if (process.contains("Features") || !feature->LoadFeature(featureFile))
	{
		feature->PrecomputeVertexFeatures();
	//	if (process.contains("Features"))
			feature->WriteFeature(featureFile);
	}
	
	return feature;
}


SurfaceSampleSet * AnalysisPipeline::CreateSampleSet(const VKString & setName, SampledSurface * surface)
{
	VKString meshName = surface->GetName();
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);
	
	VKString type = m_Params.GetStrValue(setName, "Type", valid);
	VKString setExt = m_Params.GetStrValue(setName, "Extension", valid);
	if (setName=="BenchmarkQuery")
		setExt=".smplset.benchmark";
	VKString setFile = LocalFilename(meshName, setExt);
	
	
	//assert(valid);
	SurfaceSampleSet * sampleSet = NULL;
	
	if (!process.contains("Samples") && setName!="AllVertices")
	{
		sampleSet = new SurfaceSampleSet();
		if (!sampleSet->LoadSet(setFile, surface->GetMesh()))
		{
			delete sampleSet;
			sampleSet = NULL;
		}
	}
	
	if (sampleSet==NULL)
	{
		if (type=="SurfaceSampleSet_Extrema"
			|| type=="SurfaceSampleSet_Max"
			|| type=="SurfaceSampleSet_Min")
		{
			VKString featureName = m_Params.GetStrValue(setName, "CorrFeature", valid);
			assert(valid);
			SurfaceFeature * feature = surface->GetFeature(featureName);
			assert(feature!=NULL);
			
			//int numSmoothing = m_Params.GetIntValue(setName, "SmoothingIters", valid);
			assert(valid);
			double sizeRing = surface->AdjustedRadius(m_Params.GetDoubleValue(setName, "ExtremaRingRadius", valid));
			assert(valid);
			int minNumSamples = m_Params.GetIntValue(setName, "MinSamples", valid);
			assert(valid);
			int maxNumSamples = m_Params.GetIntValue(setName, "MaxSamples", valid);
			assert(valid);			
			
			if (type.endsWith("_Extrema"))
				sampleSet = new SurfaceSampleSet(feature, sizeRing, minNumSamples, maxNumSamples,
												 SurfaceSampleSet::FEAT_ALL_EXTREMA);
			else if (type.endsWith("_Max"))
				sampleSet = new SurfaceSampleSet(feature, sizeRing, minNumSamples, maxNumSamples,
												 SurfaceSampleSet::FEAT_MAXIMA_ONLY);
			else if (type.endsWith("_Min"))
				sampleSet = new SurfaceSampleSet(feature, sizeRing, minNumSamples, maxNumSamples,
												 SurfaceSampleSet::FEAT_MINIMA_ONLY);
			else
				assert(false);
		}
		else if (type=="SurfaceSampleSet_IterFurthest")
		{
			int numSamples = m_Params.GetIntValue(setName, "NumSamples", valid);
			VKString sampleSetName = m_Params.GetStrValue(setName, "InitialSet", valid);
			if (!valid)
				sampleSetName = "none";
			assert(sampleSetName=="none");	// TODO: allow other sets
			std::vector<int> evenVertices;
			FeatureMGD::EvenlySpreadVertices(surface, numSamples, evenVertices);
			sampleSet = new SurfaceSampleSet(surface->GetMesh(), evenVertices);			
		}
		else if (type=="SurfaceSampleSet_MGD")
		{
			int maxIterations = m_Params.GetIntValue(setName, "MaxIterations", valid);
			double tau = m_Params.GetDoubleValue(setName, "IterationTau", valid);
			int maxSamples = m_Params.GetIntValue(setName, "MaxSamples", valid);		
			VKString sampleSetName = m_Params.GetStrValue(setName, "InitialSet", valid);
			
			sampleSet = FeatureMGD::ConstructSymmetryInvariantSet(surface, maxIterations, tau, 
																  maxSamples, sampleSetName);			
		}
		else if (type=="External")
		{
			VKStringList externalStrList = m_Params.GetStrValues(setName, "Filenames", valid);
			assert(valid);
			VKString externalStr;
			if (NumSurfaces()<=1 || surface==GetSurface(0))
				externalStr = externalStrList[0];
			else 
				externalStr = externalStrList[1];
			
			sampleSet = new SurfaceSampleSet();
			std::ifstream textStream(externalStr.c_str());
			assert(textStream.is_open());
			while(true)
			{
//				std::string vID;
//				textStream>>vID;
//				VKString tempStr(vID);
				VKString tempStr;
				tempStr.readLine(textStream);
				bool ok;
				int intID = tempStr.toInt(&ok);
				if (ok)
					sampleSet->AddSample(SurfaceSample(intID, surface->GetMesh()));
				if (textStream.eof())
					break;
			}
		}		
		else if (setName=="BenchmarkQuery")
		{
			VKString truthDescriptor = m_Params.GetStrValue("Pipeline", "BenchmarkQuery", valid);
			assert(valid);
			VKString truthType = m_Params.GetStrValue(truthDescriptor, "MapType", valid);
			if (truthType=="CoarseMap_FeaturePointsCorrespondence" 
				|| truthType=="CoarseMap_FeaturePointsSymmetry"
				|| truthType=="CoarseMap_SameVertexIDs"
				|| truthType=="CoarseMap_CreateTrueAndSymmetric")
			{
				int surfaceID = NumSurfaces();
				VKStringList featureVerticesNames = m_Params.GetStrValues(truthDescriptor, "FileVertexIDs", valid);
				assert(valid);
				VKStringList featureVerticesPaths = m_Params.GetStrValues(truthDescriptor, "Path", valid);
				assert(valid);
				if (featureVerticesPaths.count()==1)
					featureVerticesPaths.push_back(featureVerticesPaths[0]);
	//			std::cout<<surfaceID<<" <? "<<featureVerticesPaths.count()<<std::endl;
				
				assert(surfaceID < featureVerticesPaths.count());
				assert(surfaceID < featureVerticesNames.count());				
				VKString truthFile=featureVerticesPaths[surfaceID]+"/"+featureVerticesNames[surfaceID];
				assert(featureVerticesNames.count()>surfaceID);
				sampleSet = MapCoarse::GetSetFeaturePoints(surface, truthFile);
			}
			else
				assert(false);
		}
		else if (setName=="AllVertices")
		{
			sampleSet = new SurfaceSampleSet(surface->GetMesh());
		}
		else
		{
			std::cout<<"[ERROR] Unknown Set Type: "<<type.c_str()<<" Name="<<setName.c_str()<<std::endl;
			assert(false);
		}
		
		assert(sampleSet!=NULL);
		if (setName!="AllVertices")
			sampleSet->WriteSet(setFile);
	}
	
	surface->AddSampleSet(setName, sampleSet);
	
	return sampleSet;	
}

SurfaceDistance * AnalysisPipeline::CreateDistanceMetric(const VKString & metricDescriptor, 
														 SampledSurface * surface)
{
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);	
	SurfaceDistance * newDistance = NULL;
	VKString meshName = surface->GetName();
	VKString metricType = m_Params.GetStrValue(metricDescriptor, "Type", valid);
	
	VKString distanceExt = m_Params.GetStrValue(metricDescriptor, "Extension", valid);
	VKString distanceFile = LocalFilename(meshName, distanceExt);
	
	assert(valid);
	if (metricType=="DistanceGeodesic")
	{
		VKString calcType = m_Params.GetStrValue(metricDescriptor, "Calculation", valid);
		newDistance = new DistanceGeodesic(surface, DistanceGeodesic::GetTypeFromStr(calcType));
	}
	else if (metricType=="DistanceOnTheFly")
	{
		VKString calcType = m_Params.GetStrValue(metricDescriptor, "Calculation", valid);
		bool tryDefault = m_Params.GetStrValue(metricDescriptor, "TryDefault", valid)=="true" && valid;
		int numRowsCached = m_Params.GetIntValue(metricDescriptor, "NumCachedVertices", valid);
		
		SurfaceDistance * tryFirst=NULL;
		if (tryDefault)
			tryFirst = surface->GetDistanceMetric("default");
		newDistance = new DistanceOnTheFly(surface, DistanceGeodesic::GetTypeFromStr(calcType),
										   tryFirst, numRowsCached);
	}
	else if (metricType=="DistanceLazySubsets")
	{
		VKString calcType = m_Params.GetStrValue(metricDescriptor, "Calculation", valid);
		VKStringList setFromNames = m_Params.GetStrValues(metricDescriptor, "SetsFrom", valid);
		VKStringList setToNames = m_Params.GetStrValues(metricDescriptor, "SetsTo", valid);
		SurfaceSampleSet * setsFrom = NULL;
		SurfaceSampleSet * setsTo = NULL;
		if (setFromNames.count()==0 || setFromNames[0]=="all" || setFromNames[0]=="none")
			setsFrom = new SurfaceSampleSet(surface->GetMesh());
		else
		{
			setsFrom = new SurfaceSampleSet();
			for (int i=0; i<setFromNames.count(); i++)
			{
				SurfaceSampleSet * currSet = surface->GetSampleSet(setFromNames[i]);
				
				if (currSet==NULL)
				{
					std::cout<<"[ERROR] Could not find set (AnalysisPipeline): "<<setFromNames[i].c_str()<<std::endl;
					assert(currSet!=NULL);
				}
				setsFrom->Union(*currSet);
			}
		}
		
		if (setToNames.count()==0 || setToNames[0]=="all" || setToNames[0]=="none")
			setsTo = new SurfaceSampleSet(surface->GetMesh());
		else 
		{
			setsTo = new SurfaceSampleSet();
			for (int i=0; i<setToNames.count(); i++)
			{
				SurfaceSampleSet * currSet = surface->GetSampleSet(setToNames[i]);
				assert(currSet!=NULL);				
				setsTo->Union(*currSet);
			}
		}
		newDistance = new DistanceLazySubsets(surface, setsFrom, setsTo,
											  DistanceGeodesic::GetTypeFromStr(calcType));
	}
	else
		assert(false);
	
	if (process.contains("Distances") || !newDistance->LoadPrecomputed(distanceFile))
	{
		newDistance->PrecomputeDistances();
	//	if (process.contains("Distances"))
	//		newDistance->WritePrecomputed(distanceFile);
	}

	return newDistance;
}

SampleSetMapGenerator * AnalysisPipeline::CreateSetSampler(const VKString & samplerDescriptor, 
														   const VKString & sampleSetName,
														   std::vector<SampledSurface *> & surfaces)
{
//	m_Params.PrintParams();
	VKString type = m_Params.GetStrValue(samplerDescriptor, "Type", valid);
//	std::cout<<"Getting: "<<samplerDescriptor.c_str()<<"->Type valid="<<valid<<std::endl;
	SampleSetMapGenerator * sampler=NULL;
	if (type=="SymmetricSetSampler" || type=="SetSampler" || type=="SetSamplerAdaptive")
	{
		assert(surfaces.size()==2);
		SampledSurface * surface1 = surfaces[0];
		SampledSurface * surface2 = surfaces[1];
		
		int numCorr = m_Params.GetIntValue(samplerDescriptor, "NumCorrespondences", valid);
		assert(valid || type=="SetSamplerAdaptive");
		VKString iterType = m_Params.GetStrValue(samplerDescriptor, "IterationType", valid);
		assert(valid || type=="SetSamplerAdaptive");
		int maxNumMaps = m_Params.GetIntValue(samplerDescriptor, "MaxNumSamples", valid);
//		if (!valid)
//			maxNumMaps = 100000;
		SurfaceSampleSet * set1 = surface1->GetSampleSet(sampleSetName);
		SurfaceSampleSet * set2 = surface2->GetSampleSet(sampleSetName);
		
		// Initialize pruning
		double minDistThresh = m_Params.GetDoubleValue(samplerDescriptor, "MinDistThresh", valid);		
		if (!valid)
			minDistThresh = 0;
		
		double distToOtherThresh = m_Params.GetDoubleValue(samplerDescriptor, "DistToOtherThresh", valid);				
		if (!valid)
			distToOtherThresh = 0;

		SurfaceDistance * distance1 = surface1->GetDistanceMetric("default");		
		SurfaceDistance * distance2 = surface2->GetDistanceMetric("default");		
		assert(distance1!=NULL && distance2!=NULL);
		
		VKStringList featureThreshStr = m_Params.GetStrValues(samplerDescriptor, "FeatureThresh", valid);
		SurfaceFeature * invariantFeature1 = NULL;
		SurfaceFeature * invariantFeature2 = NULL;
		double invariantFeatureThreshVal = 0;
		if (valid && featureThreshStr.count()==2)
		{	
			invariantFeature1 = surface1->GetFeature(featureThreshStr[0]);
			invariantFeature2 = surface2->GetFeature(featureThreshStr[0]);
			invariantFeatureThreshVal = featureThreshStr[1].toDouble(&valid);
		}

		SampleSetMapGenerator::GenerationMethod method = SampleSetMapGenerator::GEN_EXHAUSTIVE_SEARCH;
		if (iterType=="Exhaustive")
			method = SampleSetMapGenerator::GEN_EXHAUSTIVE_SEARCH;
		else if (iterType=="Random")
			method = SampleSetMapGenerator::GEN_RANDOM_SEARCH;	
		else
			assert(type=="SetSamplerAdaptive");
		if (type=="SymmetricSetSampler")
		{
			sampler = new MapGeneratorSymmetric(set1, numCorr, maxNumMaps, method, surface1->Area());
			if (invariantFeature1!=NULL)
				((MapGeneratorSymmetric*)sampler)->AddInvariantFunctionThreshold(invariantFeature1, 
																				 invariantFeatureThreshVal);
			((MapGeneratorSymmetric*)sampler)->AddDistancetoToOtherPoints(distance1, distToOtherThresh);
			((MapGeneratorSymmetric*)sampler)->AddMinDistanceThreshold(distance1, minDistThresh);			
		}
		else if (type=="SetSampler")
		{
			sampler = new SampleSetMapGenerator(set1, set2, numCorr, maxNumMaps, method, 
												surface1->Area(), surface2->Area());
			if (invariantFeature1!=NULL && invariantFeature2!=NULL)
				sampler->AddInvariantFunctionThreshold(invariantFeature1, 
													   invariantFeature2, invariantFeatureThreshVal);
			sampler->AddDistancetoToOtherPoints(distance1, -1, distance2, -1, distToOtherThresh);
			sampler->AddMinDistanceThreshold(distance1, -1, distance2, -1, minDistThresh);			
		}
		else if (type=="SetSamplerAdaptive")
		{
			VKString fullSampler = m_Params.GetStrValue(samplerDescriptor, "FullSampler", valid);			
			assert(valid);
			SampleSetMapGenerator * fullGen = CreateSetSampler(fullSampler, sampleSetName, surfaces);
			VKString otherType = m_Params.GetStrValue(fullSampler, "Type", valid);
			assert(valid);
			int K = m_Params.GetIntValue(samplerDescriptor, "SuggestedNumSampC", valid); 
			if (valid && otherType!="SymmetricSetSampler" && K>0)
			{
				int M = set1->NumSamples();
				int N = set2->NumSamples();			
				int minN = vkMin(M, N);
				double numSugg = (double)K * pow((double)N*M, 3.) / pow((double)minN, 3.);
				maxNumMaps = vkMin((int)ceil(numSugg), (int)maxNumMaps);
			}
			
			sampler = new MapGeneratorAdaptive(fullGen, maxNumMaps);
		}
		else
			assert(false);
	}
	else
	{
		std::cout<<"[ERROR] Unknown set sampler type "<<type.c_str()<<std::endl;
		assert(false);
	}

//	m_Params.PrintParams();
	sampler->StartGeneration();
	
	return sampler;
}


void AnalysisPipeline::AddClonableConfidence(const VKString & mapEvalDescriptor)
{	
//	std::cout<<"Add Clonable Confidence: "<<mapEvalDescriptor.c_str()<<std::endl;
	if (mapEvalDescriptor=="none")
		return;

	if (SurfaceMapConfidence::GetClonableSeed(mapEvalDescriptor)!=NULL)
		return;
		
	VKString type = m_Params.GetStrValue(mapEvalDescriptor, "Type", valid);
	VKString atSamples = m_Params.GetStrValue(mapEvalDescriptor, "AtSamples", valid);
	if (!valid)
	{
		atSamples = m_Params.GetStrValue(mapEvalDescriptor, "INPUT_AtSamples", valid);
		if (!valid)
		{
			atSamples = "none";
			std::cout<<"[WARNING] AtSamples not defined for "<<mapEvalDescriptor.c_str()<<std::endl;
		}
	}
		
	
	VKString atSamplesRange = m_Params.GetStrValue(mapEvalDescriptor, "AtSamplesRange", valid); 
	if (!valid)
	{
		atSamplesRange = m_Params.GetStrValue(mapEvalDescriptor, "INPUT_AtSamplesRange", valid); 
		if (!valid)
		{
			atSamplesRange = "none";
			std::cout<<"[WARNING] AtSamplesRange not defined for "<<mapEvalDescriptor.c_str()<<std::endl;	
		}
	}
	
	SurfaceMapConfidence * mapEvaluator=NULL;
	if (type=="MapConfidenceMutuallyClosest")
	{
		MapConfidenceMutuallyClosest::ValueType mapVal = MapConfidenceMutuallyClosest::CMV_FRAC_MCN;
		VKString mapValStr = m_Params.GetStrValue(mapEvalDescriptor, "MapValue", valid);
		if (mapValStr=="MCN_Fraction")
			mapVal = MapConfidenceMutuallyClosest::CMV_FRAC_MCN;			
		else if (mapValStr=="INV_Distance")
			mapVal = MapConfidenceMutuallyClosest::CMV_INV_DIST;
		else
			std::cout<<"[WARNING] AnalysisPipeline.cpp: Map Value is not defined!"<<std::endl;
		
		mapEvaluator = new MapConfidenceMutuallyClosest(NULL, atSamples, atSamples, mapVal);
	}
	else if (type=="SameSizeLocalNhdsValue")
	{
		assert(false);
//		SameSizeLocalNhdsValue::ValueType valueType;
//		
//		if (m_Params.GetStrValue(mapEvalDescriptor, "ValueType", valid)=="NHD_VAL_DIFF")
//			valueType = SameSizeLocalNhdsValue::NHD_VAL_DIFF;
//		else if (m_Params.GetStrValue(mapEvalDescriptor, "ValueType", valid)=="NHD_VAL_MIN_OVER_MAX")
//			valueType = SameSizeLocalNhdsValue::NHD_VAL_MIN_OVER_MAX;
//		else
//			assert(false);
//		
//		mapEvaluator = new SameSizeLocalNhdsValue(m_Params.GetDoubleValue(mapEvalDescriptor, 
//							"NhdRadius", valid), set1, set2, fine2Coarse, valueType, dist1, dist2);
	}
	else if (type=="MapConfidenceFaceArea")
	{
		MapConfidenceFaceArea::ValueType type;
		VKString valueType = m_Params.GetStrValue(mapEvalDescriptor, "ValueType", valid);
		if (!valid)
			valueType = m_Params.GetStrValue(mapEvalDescriptor, "INPUT_ValueType", valid);
		if (valueType=="TRIAREA_VAL_DIFF")
			type = MapConfidenceFaceArea::TRIAREA_VAL_DIFF;
		else if (valueType=="TRIAREA_VAL_MIN_OVER_MAX")
			type = MapConfidenceFaceArea::TRIAREA_VAL_MIN_OVER_MAX;
		else if (valueType=="TRISTRETCH_SV_RATIO")
			type = MapConfidenceFaceArea::TRISTRETCH_SV_RATIO;
		else if (valueType=="TRISTRETCH_SV_L2STRETCH") 
			type = MapConfidenceFaceArea::TRISTRETCH_SV_L2STRETCH;
		else if (valueType=="TRISTRETCH_SV_AREA") 
			type = MapConfidenceFaceArea::TRISTRETCH_SV_AREA;
		else
			assert(false);

		mapEvaluator = new MapConfidenceFaceArea(NULL, type, atSamples, atSamplesRange);		
	}
	else if (type=="MapConfidenceDistanceToGen")
	{
		mapEvaluator = new MapConfidenceDistanceToGen(NULL, "default", "default", atSamples, atSamples);
	}
	else if (type=="MapConfidenceCompareToTruth")
	{
		LoadTruthMap();
		double truthThresh = m_Params.GetDoubleValue(mapEvalDescriptor, "TruthThreshold", valid); 
		assert(valid);
		
		mapEvaluator = new MapConfidenceCompareToTruth(NULL, m_pTruthMap, "default", "default", 
													   atSamples, atSamplesRange, truthThresh, true);
	}
	else if (type=="MapConfidenceCompareToTruthWithSym")
	{
		LoadTruthMap();
		double truthThresh = m_Params.GetDoubleValue(mapEvalDescriptor, "TruthThreshold", valid); 
		assert(valid);
		
		mapEvaluator = new MapConfidenceCompareToTruth(NULL, m_pAlternativeTruth, "default", "default", 
														atSamples, atSamplesRange, truthThresh, true);
		SurfaceMapConfidence::SetClonableSeed(mapEvalDescriptor+"_NotFlipped", mapEvaluator);
		mapEvaluator = new MapConfidenceCompareToTruth(NULL, m_pAlternativeSymmetricTruth, 
													   "default", "default", atSamples, atSamplesRange, 
													   truthThresh, true);
		SurfaceMapConfidence::SetClonableSeed(mapEvalDescriptor+"_Flipped", mapEvaluator);
	}
	else if (type=="MapConfidenceMultiConf")
	{
		VKStringList confidences = m_Params.GetStrValues(mapEvalDescriptor, "Confidences", valid); 
		VKString similarity = m_Params.GetStrValue(mapEvalDescriptor, "Similarity", valid); ;
		if (!valid)
			similarity = "none";
		VKString weightCalcStr = m_Params.GetStrValue(mapEvalDescriptor, "WeightsCalculation", valid);
		assert(valid);
		MapConfidenceMultiConf::WeightCalculation weightCalc;
		weightCalc = MapConfidenceMultiConf::StrToWeightCalculation(weightCalcStr);
		std::cout<<"Creating multi-conf confidence: "<<confidences.count()<<std::endl;
		mapEvaluator = new MapConfidenceMultiConf(NULL, atSamples, atSamplesRange, 
												  confidences, similarity, weightCalc);
	}
	else if (type=="MapConfidenceConstant")
	{
//		mapEvaluator = new MapConfidenceConstant(NULL);
		assert(false);
	}
	else
		assert(false);
	if (mapEvaluator!=NULL)	// might have been added
		SurfaceMapConfidence::SetClonableSeed(mapEvalDescriptor, mapEvaluator);
}

void AnalysisPipeline::AddClonableSimilarity(const VKString & descriptor)
{
//	std::cout<<"Add Clonable Confidence: "<<descriptor.c_str()<<std::endl;	
	if (descriptor=="none")
		return;
	if (SurfaceMapSimilarity::GetClonableSeed(descriptor)!=NULL)
		return;
	
	VKString similarityType = m_Params.GetStrValue(descriptor, "Type", valid);	assert(valid);
	VKString overlapSet = m_Params.GetStrValue(descriptor, "OverlapSet", valid);	assert(valid);
	double simThreshold = m_Params.GetDoubleValue(descriptor, "SimilarityThresh", valid);
	if (!valid)
		simThreshold = 0;
	
	VKString distOnSurface = m_Params.GetStrValue(descriptor, "DistOnSurface", valid);	
	if (!valid)
		distOnSurface = "default";
	
	SurfaceMapSimilarity * similarity=NULL;
	
	if (similarityType=="MapSimilarityOverlap")
	{
		VKString coarsening = m_Params.GetStrValue(descriptor, "CoarseningMethod", valid);	assert(valid);	
		
		similarity = new MapSimilarityOverlap(NULL, NULL, overlapSet, overlapSet, overlapSet, overlapSet);
	}
	else if (similarityType=="MapSimilarityDistance")
	{
		int flags = 0;
		if (m_Params.GetStrValues(descriptor, "Flags", valid).contains("DisallowInconsistent"))
			flags |= MapSimilarityDistance::SIMFLAG_NO_INCONSISTENT_GENS;
		if (m_Params.GetStrValues(descriptor, "Flags", valid).contains("ScaleByConfidence"))
			flags |= MapSimilarityDistance::SIMFLAG_SCALE_BY_CONFIDENCE;
		
		double maxDistance = m_Params.GetDoubleValue(descriptor, "MaxDistance", valid);
		if (!valid)
			maxDistance = -1;
		
		double sigma = m_Params.GetDoubleValue(descriptor, "Sigma", valid);
		assert(valid);
		//	MAPSIM_FORWARD_ONLY			MAPSIM_FAST_DIST_SUBSETS	MAPSIM_FAST_DIST_SUBSETS_BIDIR
		similarity = new MapSimilarityDistance(NULL, NULL, 
											   MapSimilarityDistance::MAPSIM_FAST_DIST_SUBSETS_BIDIR,
											   flags, maxDistance, simThreshold, sigma,
											   overlapSet, overlapSet, distOnSurface, distOnSurface);
		
//		simThreshold = exp(-sqrt(GetSurface(1)->AdjustedRadius(simThreshold)));
	}
	assert(similarity!=NULL);
//	similarity->SetSimThresh(simThreshold);
	
	SurfaceMapSimilarity::SetClonableSeed(descriptor, similarity);
}

MapScoredCollection * AnalysisPipeline::GetCurrCollection()
{
	return m_pCollectionOfMaps;
}

VKString AnalysisPipeline::LocalFilename(const VKString & meshname, const VKString & ext)
{
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);	assert(valid);	
	VKString algOutputName = m_Params.GetStrValue(algorithmName, "OutputName", valid);	assert(valid);
	return LocalFilename(meshname, algOutputName+"_"+meshname, ext);
}

VKString AnalysisPipeline::LocalFilename(const VKString & ext)
{
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);	assert(valid);	
	VKString surfDescr = m_Params.GetStrValue(algorithmName, "Surface", valid);	assert(valid);
	VKStringList meshNames = m_Params.GetStrValues(surfDescr, "MeshName", valid);	assert(valid);
	if (meshNames.count()==1)
		meshNames.push_back(meshNames[0]);
	assert(meshNames.count()==2);	// for now correspondence between two only		

	VKString algOutputName = m_Params.GetStrValue(algorithmName, "OutputName", valid);	assert(valid);
	return LocalFilename(meshNames[0], algOutputName+"_Map_"+meshNames[0]+"_to_"+meshNames[1], ext);
}

VKString AnalysisPipeline::LocalFilename(const VKString & meshname, 
										 const VKString & filename, 
										 const VKString & ext)
{
	VKString subdir="";
	if (m_Params.GetStrValue("Pipeline", "Subdirs", valid)=="DirByMeshname")
		subdir=meshname;
		
	return m_Params.GetStrValue("Pipeline", "WorkFolder", valid) + "/"+subdir+"/"+filename+ext;
}

int AnalysisPipeline::NumSurfaces()
{
	return (int)m_Surfaces.size();
}

SampledSurface * AnalysisPipeline::GetSurface(int i)
{
	assert(i>=0 && i<(int)m_Surfaces.size());
	return m_Surfaces[i];
}

////////////////// RENDERING /////////////////

void AnalysisPipeline::SetWindow(AnalysisWindow * window)
{
	m_pWindow = window;
	
	for (int i=0; i<(int)m_Surfaces.size(); i++)
	{
		m_Surfaces[i]->SetParentWindow(window);
		if (i < (int)m_InitCameraPositions.size() && m_InitCameraPositions[i]!=NULL)
		{
			window->SetMeshCamera(i, 0, 
					m_InitCameraPositions[i][0], m_InitCameraPositions[i][1], m_InitCameraPositions[i][2], 
					m_InitCameraPositions[i][3], m_InitCameraPositions[i][4], m_InitCameraPositions[i][5], 
					m_InitCameraPositions[i][6], m_InitCameraPositions[i][7], m_InitCameraPositions[i][8]);
//			window->ScaleMesh(i, 0, 0.75);
			delete [] m_InitCameraPositions[i];
			m_InitCameraPositions[i] = NULL;
		}
	}
	m_InitCameraPositions.clear();
}

void AnalysisPipeline::Draw(ParamParser & drawParams)
{
	assert(false);	// should be implemented in a child class
}

void AnalysisPipeline::WriteScreenshot(const R2Image & img, const VKString & filename)
{
	VKString algorithmName = m_Params.GetStrValue("Pipeline", "Algorithm", valid);
	assert(valid);

	VKString algNameOut = m_Params.GetStrValue(algorithmName, "OutputName", valid);
	if (!valid)
		algNameOut = algorithmName;
	
	VKString surfDescr = m_Params.GetStrValue(algorithmName, "Surface", valid);
	assert(valid);	
	VKStringList meshNames = m_Params.GetStrValues(surfDescr, "MeshName", valid);
	assert(valid);	
	if (meshNames.count()==1)
		img.WriteJPEG(LocalFilename(meshNames[0], algNameOut+"_"+filename, ".jpg").c_str());
	else if (meshNames.count()==2)
		img.WriteJPEG(LocalFilename(meshNames[0], 
									algNameOut+"_"+meshNames[0]+"_to_"+meshNames[1]+"_"+filename, ".jpg").c_str());
	else 
		assert(false);
}

	


