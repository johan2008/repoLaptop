#include "PipelineGeneral.h"
#include "SampledSurface.h"
#include "MapFlattening.h"
#include "FeatureAGD.h"
#include "FeatureMGD.h"
#include "FeatureCurvature.h"
#include "FeatureIRSA.h"
#include "PlanarTransformMVC.h"
#include "PlanarTransformLSCM.h"
#include "PlanarTransformQuasiConformal.h"
#include "Sortable.h"

PipelineGeneral * PipelineGeneral::s_pGeneralPipeline = NULL;

PipelineGeneral::PipelineGeneral(ParamParser & params)
: AnalysisPipeline(params)
{
	s_pGeneralPipeline = this;
	VKString pipelineType = params.GetStrValue("Pipeline", "Algorithm", valid);
	assert(valid);
	assert(pipelineType=="PipelineGeneral");
	
	VKStringList loadModules = params.GetStrValues("Pipeline", "Load", valid);
	assert(valid);
	VKStringList processModules = params.GetStrValues("Pipeline", "Process", valid);
	assert(valid);
	for (int i=0; i<loadModules.count(); i++)
	{
		std::cout<<"Pipeline Step: "<<loadModules[i].c_str()<<std::endl;
		bool executed = LoadAndProcessModule(loadModules[i], 
											 processModules.contains(loadModules[i]));
		if (!executed)
		{
			std::cout<<"[ERROR] Could not process module: "<<loadModules[i].c_str()<<std::endl;
			assert(false);
		}
	}
	std::cout<<"Pipeline Finished"<<std::endl;
}

VKString PipelineGeneral::GetPipelineType()
{
	return "PipelineGeneral";
}

void PipelineGeneral::AddElementToSurfaces(const VKString & descriptor, bool reprocess, 
										   const VKString & moduleType)
{
	VKStringList toSurfaces = m_Params.GetStrValues(descriptor, "INPUT_Surfaces", valid);
	assert(valid);
	if (toSurfaces.contains("all"))
	{
		for (int i=0; i<(int)m_SurfaceNames.count(); i++)
		{
			if (moduleType=="Feature")
				AddFeaturesToSurface(descriptor, reprocess, m_SurfaceNames[i]);
			else if (moduleType=="Distances")
				AddDistancesToSurface(descriptor, reprocess, m_SurfaceNames[i]);
			else if (moduleType=="SampleSet")
				AddSamplesToSurface(descriptor, reprocess, m_SurfaceNames[i]);
			else
			{
				std::cout<<"[ERROR] Unknown element module type: "<<moduleType.c_str()<<std::endl;
				assert(false);
			}
		}
	}
	else
		assert(false);	//TODO	
}

void PipelineGeneral::AddDistancesToSurface(const VKString & descriptor, bool reprocess, 
											const VKString & surfaceName)
{
	SampledSurface * surface = GetSurfaceByName(surfaceName);
	assert(surface!=NULL);
	SurfaceDistance * distanceMetric=NULL;
	
	VKString metricType = m_Params.GetStrValue(descriptor, "Type", valid);
	VKString calcType = m_Params.GetStrValue(descriptor, "INPUT_Calculation", valid);	
	bool setAsDefault = m_Params.GetStrValue(descriptor, "INPUT_SetAsDefault", valid)=="true";	
	VKString memoryName = m_Params.GetStrValue(descriptor, "MEMORY_Name", valid);		
	
	if (metricType=="DistanceGeodesic")
		distanceMetric = new DistanceGeodesic(surface, DistanceGeodesic::GetTypeFromStr(calcType));
	else if (metricType=="DistanceOnTheFly")
	{
		VKString tryDefault = m_Params.GetStrValue(descriptor, "INPUT_TryDefault", valid);
		int numRowsCached = m_Params.GetIntValue(descriptor, "INPUT_NumCachedVertices", valid);
		
		SurfaceDistance * tryFirst=NULL;
		if (tryDefault!="none")
			tryFirst = surface->GetDistanceMetric(tryDefault);
		distanceMetric = new DistanceOnTheFly(surface, DistanceGeodesic::GetTypeFromStr(calcType),
											  tryFirst, numRowsCached);
	}
	else if (metricType=="DistanceLazySubsets")
	{
		VKStringList setFromNames = m_Params.GetStrValues(descriptor, "INPUT_SetsFrom", valid);
		VKStringList setToNames = m_Params.GetStrValues(descriptor, "INPUT_SetsTo", valid);
		SurfaceSampleSet * setsFrom = NULL;
		SurfaceSampleSet * setsTo = NULL;
		for (int sID=0; sID<2; sID++)
		{
			VKStringList & setNames = (sID==0) ? setFromNames : setToNames;
			SurfaceSampleSet ** setPtr = (sID==0) ? (&setsFrom) : (&setsTo);
			if (setNames.count()==0 || setNames[0]=="all" || setNames[0]=="none")
				*setPtr = new SurfaceSampleSet(surface->GetMesh());
			else
			{
				(*setPtr) = new SurfaceSampleSet();
				for (int i=0; i<setNames.count(); i++)
				{
					SurfaceSampleSet * currSet = surface->GetSampleSet(setNames[i]);
					if (currSet==NULL)
					{
						std::cout<<"[ERROR] Could not find set (PipelineGeneral:AddDistance): "<<setNames[i].c_str()<<std::endl;
						assert(currSet!=NULL);
					}
					(*setPtr)->Union(*currSet);
				}
			}
		}
		assert(setsFrom!=NULL && setsTo!=NULL);
		std::cout<<"Calculating Distances: "<<setsFrom->NumSamples()<<" x "<<setsTo->NumSamples()<<std::endl;
		distanceMetric = new DistanceLazySubsets(surface, setsFrom, setsTo,
												 DistanceGeodesic::GetTypeFromStr(calcType));
	}
	else
	{
		std::cout<<"[ERROR] Unknown Distance Metric: "<<metricType.c_str()<<std::endl;
		assert(false);
	}
	
	VKString filename = m_Params.GetStrValue(descriptor, "OUTPUT_Distances", valid);
	if (!valid)
		filename = "none";
	FilenameAddWrkspace(filename);
	filename.replace("!mshname!", surfaceName);
	assert(distanceMetric!=NULL);
	if (!reprocess && filename!="none" && FileExists(filename))
	{
		bool loaded = distanceMetric->LoadPrecomputed(filename);
		assert(loaded);
	}
	else
	{
		distanceMetric->PrecomputeDistances();
		if (filename!="none")
			distanceMetric->WritePrecomputed(filename);
	}
	if (setAsDefault)
		surface->AddDistanceMetric("default", distanceMetric);
	surface->AddDistanceMetric(memoryName, distanceMetric);
	if (!m_DistanceMetricNames.contains(memoryName))
		m_DistanceMetricNames.push_back(memoryName);
}

void PipelineGeneral::AddFeaturesToSurface(const VKString & descriptor, bool reprocess, 
										   const VKString & surfaceName)
{
	SampledSurface * surface = GetSurfaceByName(surfaceName);
	assert(surface!=NULL);

	SurfaceFeature * feature = NULL;
	
	VKString featureType = m_Params.GetStrValue(descriptor, "Type", valid);
	VKString distanceName = m_Params.GetStrValue(descriptor, "INPUT_Distance", valid);
	if (!valid)
		distanceName = "default";
	VKString memoryFeatureName = m_Params.GetStrValue(descriptor, "MEMORY_Name", valid);
	
	int nInit = m_Params.GetIntValue(descriptor, "INPUT_NInit", valid);
	if (!valid)
		nInit = 256;
	
	if (featureType=="FeatureAGD")
		feature = new FeatureAGD(surface, distanceName, true, memoryFeatureName);
	else if (featureType=="FeatureFastAGD")
		feature = new FeatureFastAGD(surface, nInit, memoryFeatureName);
	else if (featureType=="FeatureCurvature")
		feature = new FeatureCurvature(surface, nInit, memoryFeatureName);
	else if (featureType=="FeatureIRSA")
	{
		VKString mapCollectionName = m_Params.GetStrValue(descriptor, "INPUT_MapCollection", valid);
		assert(valid);
		VKString confidenceName = m_Params.GetStrValue(descriptor, "INPUT_Confidence", valid);
		assert(valid);
		VKString distanceMetric = m_Params.GetStrValue(descriptor, "INPUT_DistanceMetric", valid);
		if (!valid)
			distanceMetric="default";
		feature = new FeatureIRSA(surface, mapCollectionName, confidenceName, 
								  distanceMetric, memoryFeatureName);
	}
	else
	{
		std::cout<<"[ERROR] Unknown Feature Type: "<<featureType.c_str()<<std::endl;
		assert(false);
	}
	
	assert(feature!=NULL);
	VKString filename = m_Params.GetStrValue(descriptor, "OUTPUT_Feature", valid);
	if (!valid)
		filename = "none";
	FilenameAddWrkspace(filename);
	filename.replace("!mshname!", surfaceName);
	
	if (!reprocess && FileExists(filename))
	{
		bool loaded = feature->LoadFeature(filename);
		assert(loaded);
	}
	else
	{
		feature->PrecomputeVertexFeatures();
		if (filename!="none")
			feature->WriteFeature(filename);

		VKString filenameArff = m_Params.GetStrValue(descriptor, "OUTPUT_arff", valid);
		if (!valid)
			filenameArff = "none";
		FilenameAddWrkspace(filenameArff);
		filenameArff.replace("!mshname!", surfaceName);
		if (filenameArff!="none")
			feature->WriteARFFPropertyFile(filenameArff);
	}
	surface->AddFeature(memoryFeatureName, feature);
	if (!m_FeatureNames.contains(memoryFeatureName))
		m_FeatureNames.push_back(memoryFeatureName);

}

void PipelineGeneral::AddSamplesToSurface(const VKString & descriptor, bool reprocess, 
										  const VKString & surfaceName)
{
	SurfaceSampleSet * sampleSet=NULL;
	
	VKString filename = m_Params.GetStrValue(descriptor, "OUTPUT_SampleSet", valid);
	FilenameAddWrkspace(filename);
	filename.replace("!mshname!", surfaceName);
	SampledSurface * surface = GetSurfaceByName(surfaceName);
	assert(surface!=NULL);
	if (!reprocess && FileExists(filename))
	{
		sampleSet = new SurfaceSampleSet();
		if (!sampleSet->LoadSet(filename, surface->GetMesh()))
		{
			delete sampleSet;
			sampleSet = NULL;
		}
	}	
	
	if (sampleSet==NULL)
	{
		VKString type = m_Params.GetStrValue(descriptor, "Type", valid);
		if (type=="AllVertices")
			sampleSet = new SurfaceSampleSet(surface->GetMesh());
		else if (type=="IterativeFurthest")
		{
			int numSamples = m_Params.GetIntValue(descriptor, "INPUT_NumSamples", valid);
			assert (valid);
			VKString sampleSetName = m_Params.GetStrValue(descriptor, "INPUT_InitialSet", valid);
			if (!valid)
				sampleSetName = "none";
			assert(sampleSetName=="none");	// TODO: allow other sets
			std::vector<int> evenVertices;
			FeatureMGD::EvenlySpreadVertices(surface, numSamples, evenVertices);
			sampleSet = new SurfaceSampleSet(surface->GetMesh(), evenVertices);
		}
		else if (type=="Feat_Extrema" || type=="Feat_Max" || type=="Feat_Min")
		{
			VKString featureName = m_Params.GetStrValue(descriptor, "INPUT_Feature", valid);
			assert(valid);
			SurfaceFeature * feature = surface->GetFeature(featureName);
			assert(feature!=NULL);
			
			assert(valid);
			double sizeRing = surface->AdjustedRadius(m_Params.GetDoubleValue(descriptor, "INPUT_ExtremaRingRadius", valid));
			assert(valid);
			int minNumSamples = m_Params.GetIntValue(descriptor, "INPUT_MinSamples", valid);
			assert(valid);
			int maxNumSamples = m_Params.GetIntValue(descriptor, "INPUT_MaxSamples", valid);
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
		else 
		{
			std::cout<<"[ERROR] Unknown sample set type: "<<type.c_str()<<std::endl;
			assert(false);
		}
		
		if (filename!="none")
			sampleSet->WriteSet(filename);
	}
	assert(sampleSet!=NULL);
	VKString sampleSetName = m_Params.GetStrValue(descriptor, "MEMORY_Name", valid);
	surface->AddSampleSet(sampleSetName, sampleSet);
	if (!m_SampleSetNames.contains(sampleSetName))
		m_SampleSetNames.push_back(sampleSetName);
	std::cout<<"\t... "<<sampleSet->NumSamples()<<" samples on surface "<<surfaceName.c_str()<<std::endl;
}

bool PipelineGeneral::LoadAndProcessModule(VKString & moduleName, bool reprocess)
{
	if (moduleName.startsWith("#"))
		return true;
	VKString moduleType = m_Params.GetStrValue(moduleName, "ModuleType", valid);
	assert(valid);
	if (moduleType=="Surface")
		AddSurfaces(moduleName, reprocess);
	else if (moduleType=="SurfaceMap")
		AddMaps(moduleName, reprocess);
	else if (moduleType=="Texture")
		AddTextures(moduleName, reprocess);
	else if (moduleType=="MapConfidence")
		AddConfidences(moduleName, reprocess);
	else if (moduleType=="SampleSet" || moduleType=="Feature" || moduleType=="Distances")
		AddElementToSurfaces(moduleName, reprocess, moduleType);
	else
	{
		std::cout<<"[ERROR] Unknown pipeline module type: "<<moduleType.c_str()<<std::endl;
		assert(false);
	}
	return true;
}

void PipelineGeneral::AddSurfaces(const VKString & descriptor, bool )
{
	VKString instance = m_Params.GetStrValue(descriptor, "Type", valid);
	assert(valid);
	if (instance=="SampledSurface")
	{
		VKStringList inputMeshFiles = m_Params.GetStrValues(descriptor, "INPUT_MeshNames", valid);
		assert(valid);
		VKString memNameTmplt = m_Params.GetStrValue(descriptor, "MEMORY_Name", valid);
		assert(valid);
		for (int i=0; i<inputMeshFiles.count(); i++)
		{
			VKString meshName=GetBasename(inputMeshFiles[i]);			
			VKString currMemTmplt = memNameTmplt;	
			currMemTmplt.replace("!mshname!", meshName);

			// Load mesh
			R3Mesh * mesh = new R3Mesh();
			mesh->ReadFile(inputMeshFiles[i].c_str());

			// Create surface
			SampledSurface * newSurface = new SampledSurface(mesh);
			
			m_mSurfaceNameToID[currMemTmplt] = m_Surfaces.size();
			m_Surfaces.push_back(newSurface);
			m_SurfaceNames.push_back(currMemTmplt);
			
			// Load cameras
			VKStringList cameraFiles = m_Params.GetStrValues(descriptor, "OPT_CameraFile", valid);
			VKStringList cameraPositions = m_Params.GetStrValues(descriptor, "OPT_CameraOrigin", valid);
			VKStringList cameraOrientations = m_Params.GetStrValues(descriptor, "OPT_CameraTowards", valid);
			VKStringList cameraUps = m_Params.GetStrValues(descriptor, "OPT_CameraUp", valid);	
			VKStringList cameraLines = m_Params.GetStrValues(descriptor, "OPT_CameraLine", valid);
			AddSurfaceCamera(cameraFiles, cameraPositions, cameraOrientations, cameraUps,
							 cameraLines, (int)m_Surfaces.size()-1);
		}
	}
	else
	{
		std::cout<<"[ERROR] Unknown surface instance type: "<<instance.c_str()<<std::endl;
		assert(false);
	}
}

void PipelineGeneral::AddTextures(const VKString & descriptor, bool )
{
	VKString instance = m_Params.GetStrValue(descriptor, "Type", valid);
	assert(valid);
	if (instance=="R2Image")
	{
		VKStringList textureFiles = m_Params.GetStrValues(descriptor, "INPUT_Filenames", valid);
		assert(valid);
		for (int i=0; i<textureFiles.count(); i++)
		{
			VKString textName=GetBasename(textureFiles[i]);
			VKString memoryName = m_Params.GetStrValue(descriptor, "MEMORY_Name", valid);
			memoryName.replace("!txname!", textName);
			SurfaceTexture * newTx = new SurfaceTexture(textName, new R2Image(textureFiles[i].c_str()));
			
			m_mTextureNameToID[memoryName] = (int)m_Textures.size();
			m_Textures.push_back(newTx);
			m_TextureNames.push_back(memoryName);
		}
	}	
	else
	{
		std::cout<<"[ERROR] Unknown texture instance type: "<<instance.c_str()<<std::endl;
		assert(false);
	}
}

void PipelineGeneral::AddMaps(const VKString & descriptor, bool reprocess)
{
	VKStringList surfacesFrom = m_Params.GetStrValues(descriptor, "INPUT_SurfacesFrom", valid);
	assert(valid);
	VKStringList surfacesTo = m_Params.GetStrValues(descriptor, "INPUT_SurfacesTo", valid);
	assert(valid);
	bool mapToSelf = surfacesTo.contains("self");

	// Find map between surfaces
	if (surfacesFrom.contains("all"))
	{
		surfacesFrom.clear();
		surfacesFrom = m_SurfaceNames;
	}
	else if (surfacesFrom.contains("none"))
		surfacesFrom.clear();
	
	if (surfacesTo.contains("all"))
	{
		surfacesTo.clear();
		surfacesTo = m_SurfaceNames;
	}
	else if (surfacesFrom.contains("none"))
		surfacesTo.clear();
	
	// Try to load maps if applicable	
	for (int i=0; i<surfacesFrom.count(); i++)
	{
		if (surfacesTo.count()==0)	// load maps that only map one surface (to a simple domain)
		{
			AddMap(descriptor, reprocess, surfacesFrom[i], "none");
		}
		else	// load general maps between two surfaces
		{
			for (int j=0; j<surfacesTo.count(); j++)
			{
				if (surfacesFrom[i]!=surfacesTo[j] || mapToSelf)
					AddMap(descriptor, reprocess, surfacesFrom[i], surfacesTo[j]);
			}
		}
	}
}

void PipelineGeneral::AddMap(const VKString & descriptor, bool reprocess,
							 const VKString & fromMeshName, const VKString & toMeshName)
{
	VKString instance = m_Params.GetStrValue(descriptor, "Type", valid);
	assert(valid);
	VKString outputMapFile = m_Params.GetStrValue(descriptor, "OUTPUT_Map", valid);
	if (!valid)
		outputMapFile = "none";

	if (outputMapFile!="none")
	{
		FilenameAddWrkspace(outputMapFile);
		outputMapFile.replace("!mshfrom!", fromMeshName);
		outputMapFile.replace("!mshto!", toMeshName);
	}
	
	SurfaceMap * loadedMap = NULL;
	if (!reprocess && outputMapFile!="none")	// try to load map
	{
		SampledSurface * fromMesh=GetSurfaceByName(fromMeshName);
		SampledSurface * toMesh=GetSurfaceByName(toMeshName);
		loadedMap = SurfaceMap::CreateMap(outputMapFile, fromMesh, toMesh);
	}

	if (instance=="MapFlatMidEdge")
	{
		assert(toMeshName=="none");
		loadedMap = AddFlatteningMap(descriptor, reprocess, instance, outputMapFile, loadedMap, fromMeshName);
	}
	else if (instance=="MapCollection")
	{
		loadedMap = AddCollectionOfMaps(descriptor, reprocess, instance, outputMapFile, loadedMap, 
										fromMeshName, toMeshName);
	}
	//TODO: we can add other types of maps here
	else
	{
		std::cout<<"[ERROR] Unknown map instance type: "<<instance.c_str()<<std::endl;
		assert(false);
	}
	
	VKString mapName = m_Params.GetStrValue(descriptor, "MEMORY_Name", valid);
	assert(valid);
	assert(loadedMap!=NULL);
	loadedMap->GetSurface(0)->AddMapActingAsFrom(mapName, loadedMap);
	if (toMeshName!="none")
		loadedMap->GetSurface(1)->AddMapActingAsTo(mapName, loadedMap);
	if (!m_SurfaceMapNames.contains(mapName))
		m_SurfaceMapNames.push_back(mapName);
	
	if (outputMapFile!="none")
		loadedMap->SaveMap(outputMapFile);
}

SurfaceMap * PipelineGeneral::AddCollectionOfMaps(const VKString & descriptor, bool reprocess, 
												  const VKString & instance, const VKString & mapFilename,
												  SurfaceMap * loadedMap, 
												  const VKString & fromMeshName, const VKString & toMeshName)
{
	SampledSurface * surface1 = GetSurfaceByName(fromMeshName);
	SampledSurface * surface2 = GetSurfaceByName(toMeshName);
	assert(surface1!=NULL && surface2!=NULL);
	
	VKString surfaceSamples = m_Params.GetStrValue(descriptor, "INPUT_GenerationSet", valid);
	assert(valid);	
	SurfaceSampleSet * genSamples1 = surface1->GetSampleSet(surfaceSamples);
	SurfaceSampleSet * genSamples2 = surface2->GetSampleSet(surfaceSamples);
	assert(genSamples1!=NULL && genSamples2!=NULL);

	VKString stationarySetName = m_Params.GetStrValue(descriptor, "INPUT_GenerationSetStationary", valid);
	if (!valid)
		stationarySetName = "none";
	SurfaceSampleSet * stationarySet1 = surface1->GetSampleSet(stationarySetName);
	SurfaceSampleSet * stationarySet2 = surface2->GetSampleSet(stationarySetName);
	
	VKString unionGenSetName = m_Params.GetStrValue(descriptor, "MEMORY_UnionGenerationSet", valid);
	assert(valid);
	if (!m_SampleSetNames.contains(unionGenSetName))
		m_SampleSetNames.push_back(unionGenSetName);
	assert(surface1->GetSampleSet(unionGenSetName)==NULL);
	assert(surface2->GetSampleSet(unionGenSetName)==NULL);
	SurfaceSampleSet * genSetUnion1 = new SurfaceSampleSet();
	genSetUnion1->Union(*genSamples1);
	if (stationarySet1!=NULL)
		genSetUnion1->Union(*stationarySet1);

	surface1->AddSampleSet(unionGenSetName, genSetUnion1);
	SurfaceSampleSet * genSetUnion2 = NULL;
	if (surface1!=surface2)
	{
		genSetUnion2->Union(*genSamples2);
		if (stationarySet2!=NULL)
			genSetUnion2->Union(*stationarySet2);
		surface2->AddSampleSet(unionGenSetName, genSetUnion2);		
	}
	
	assert(instance=="MapCollection");
	if (loadedMap!=NULL)
	{
		MapScoredCollection * mapCollection = (MapScoredCollection*)loadedMap;
		
		for (int i=0; i<mapCollection->GetNumMaps(); i++)
		{
			VKString flatteningName = m_Params.GetStrValue(descriptor, "INPUT_Flattening", valid);
			SurfaceMap * flatteningMap1 = surface1->GetMapActingAsFrom(flatteningName);
			assert(flatteningMap1!=NULL);
			assert(flatteningMap1->GetSurfaceMapType().startsWith("MapFlat"));
			R3Mesh * flat1 = ((MapFlattening*)flatteningMap1)->GetFlatMesh();
			R2Kdtree<FlatSearchNode*> * tree1 = ((MapFlattening*)flatteningMap1)->GetKdTree();
			
			SurfaceMap * flatteningMap2 = surface2->GetMapActingAsFrom(flatteningName);
			assert(flatteningMap2!=NULL);
			assert(flatteningMap2->GetSurfaceMapType().startsWith("MapFlat"));
			R3Mesh * flat2 = ((MapFlattening*)flatteningMap2)->GetFlatMesh();
			R2Kdtree<FlatSearchNode*> * tree2 = ((MapFlattening*)flatteningMap2)->GetKdTree();
			
			if (mapCollection->GetMapByID(i)->GetSurfaceMapType().startsWith("MapFlat"))
			{
				MapFlattening * flatMap = (MapFlattening*) mapCollection->GetMapByID(i);
				flatMap->SetFlatMesh(flat1, tree1);
			}
			else if (mapCollection->GetMapByID(i)->GetSurfaceMapType()=="MapVia2DPlane")
			{
				MapVia2DPlane * flatMap = (MapVia2DPlane*) mapCollection->GetMapByID(i);
				flatMap->SetFlattenings(flat1, tree1, flat2, tree2);
			}
		}
		
		assert(!reprocess);
		return loadedMap;
	}
	MapScoredCollection * mapCollection = new MapScoredCollection(surface1, surface2);
	loadedMap = mapCollection;
	
	VKString generationType = m_Params.GetStrValue(descriptor, "INPUT_GenerationType", valid);
	assert(valid);
	VKString mapTypes = m_Params.GetStrValue(descriptor, "INPUT_MapTypes", valid);
	assert(valid);
	VKString distanceName = m_Params.GetStrValue(descriptor, "INPUT_DistanceMetric", valid);
	if (!valid)
		distanceName = "default";
	SurfaceDistance * distance1 = surface1->GetDistanceMetric(distanceName);
	SurfaceDistance * distance2 = surface2->GetDistanceMetric(distanceName);	
	
	double dist1Norm = sqrt(surface1->Area());
	double dist2Norm = sqrt(surface2->Area());

	if (generationType=="PairsAndMidpointSymmetry")
	{
		double distThreshold = .1;	// for c1 & c2
		//double distThreshold = 0.;
		
		assert(distance1!=NULL && distance2!=NULL);
		assert(genSamples1==genSamples2);
		assert(distance1==distance2);
		assert(stationarySet1!=NULL && stationarySet2 !=NULL && stationarySet1==stationarySet2);
		
		int N = genSamples1->NumSamples(); 
		std::vector<int> allTriplets;
		for (int c1=0; c1<N; c1++)
		{
			SurfaceSample sc1 = genSamples1->GetSample(c1);
			for (int c2=c1+1; c2<N; c2++)	// iterate through all pairs
			{
				SurfaceSample sc2 = genSamples1->GetSample(c2);
				double d12 = distance1->Distance(sc1, sc2) / dist1Norm;
				if (d12 <= distThreshold)
					continue;
				std::vector<Sortable> allC3Candidates;
				// pick such that abs( (D(c1, c3) - D(c2,c3))  / (D(c1, c3) + D(c2,c3)) ) is minimal
				for (int c3=0; c3 < stationarySet1->NumSamples(); c3++)
				{
					SurfaceSample sc3 = stationarySet1->GetSample(c3);
					double d13 = distance1->Distance(sc1, sc3) / dist1Norm;
					double d23 = distance1->Distance(sc2, sc3) / dist1Norm;
					if (d13 <= distThreshold || d23 <= distThreshold)
						continue;
					double val = vkAbs(d13-d23);
					//double val = vkAbs(((d13-d23) / (d13+d23)));
					//double val = vkAbs(d12-d23) + vkAbs(d12-d13);
					//double val = pow(d12/2.-d23, 2.) + pow(d12/2.-d13, 2.);
					allC3Candidates.push_back(Sortable(val, NULL, c3));
				}
				std::sort(allC3Candidates.begin(), allC3Candidates.end());
				
				// how to select c3
				double takePercentile = .05;
				int maxN = floor(takePercentile * (int)allC3Candidates.size())+1;
				int numC3=100;
				//		int maxN = 1;
				//		int numC3 = 1;
				
				genSetUnion1->ContainsExact(sc1, &c1);
				genSetUnion1->ContainsExact(sc2, &c2);

				for (int i=0; i<numC3; i++)
				{ 
					int locID = rand()%maxN;
					if (locID < 0) locID *= -1.;
					int c3 = allC3Candidates[locID].id;
					
					genSetUnion1->ContainsExact(stationarySet1->GetSample(c3), &c3);
					allTriplets.push_back(c1);
					allTriplets.push_back(c2);
					
					allTriplets.push_back(c2);
					allTriplets.push_back(c1);
					
					allTriplets.push_back(c3);
					allTriplets.push_back(c3);
				}
			}
		}
		
		VKString confidenceName = m_Params.GetStrValue(descriptor, "INPUT_MapConfidenceScore", valid);
		if (!valid)
			confidenceName = "none";
		VKString mapTypes = m_Params.GetStrValue(descriptor, "INPUT_MapTypes", valid);
		std::cout<<"\tCreating Map Collection ("<<allTriplets.size()/6<<"):"<<std::endl;
		for (int i=0; i<(int)allTriplets.size(); i+=6)
		{
			if (i%10==0)
				AnalysisStats::m_GlobalStats.m_Timing.WriteProgress(VKString("MapCollection_")+descriptor, 
																	i/6, allTriplets.size()/6);
			
			SurfaceMap * newMap = NULL;
			if (mapTypes=="Flattenings")
			{
				VKString flatteningName = m_Params.GetStrValue(descriptor, "INPUT_Flattening", valid);
				SurfaceMap * flatteningMap1 = surface1->GetMapActingAsFrom(flatteningName);
				assert(flatteningMap1!=NULL);
				assert(flatteningMap1->GetSurfaceMapType().startsWith("MapFlat"));
				R3Mesh * flat1 = ((MapFlattening*)flatteningMap1)->GetFlatMesh();
				R2Kdtree<FlatSearchNode*> * tree1 = ((MapFlattening*)flatteningMap1)->GetKdTree();
				SurfaceMap * flatteningMap2 = surface2->GetMapActingAsFrom(flatteningName);
				assert(flatteningMap2!=NULL);
				assert(flatteningMap2->GetSurfaceMapType().startsWith("MapFlat"));
				R3Mesh * flat2 = ((MapFlattening*)flatteningMap2)->GetFlatMesh();
				R2Kdtree<FlatSearchNode*> * tree2 = ((MapFlattening*)flatteningMap2)->GetKdTree();
				std::vector<int> corrs;
				for (int j=0; j<6; j++)
					corrs.push_back(allTriplets[i+j]);
				newMap = new MapVia2DPlane(surface1, surface2, unionGenSetName, corrs,
										   flatteningMap1->GetSurfaceMapType(), "MobiusTransformation",
										   true, flat1, tree1, flat2, tree2);
			}
			else
			{
				std::cout<<"[ERROR] Unknown Map Type"<<std::endl;
				assert(false);
			}
			assert(newMap!=NULL);
			
			if (confidenceName=="none")
				mapCollection->AddMap(newMap, 0);
			else
				mapCollection->AddMap(newMap, confidenceName);
		}
		AnalysisStats::m_GlobalStats.m_Timing.WriteProgress(VKString("MapCollection_")+descriptor, allTriplets.size()/6, 
															allTriplets.size()/6);
		std::cout<<std::endl;
		mapCollection->SortMaps();
	}
	else
	{
		std::cout<<"[ERROR] Unknown Map Collection generation type: "<<generationType.c_str()<<std::endl;
		assert(false);
	}
	
	return loadedMap;
}

SurfaceMap * PipelineGeneral::AddFlatteningMap(const VKString & descriptor, bool reprocess, 
											   const VKString & instance, const VKString & mapFilename,
											   SurfaceMap * loadedMap, const VKString & fromMeshName)
{
	R3Mesh * flatMesh = NULL;
	R2Kdtree<FlatSearchNode*> * flatSearch=NULL;
	
	VKString flatMeshFile = m_Params.GetStrValue(descriptor, "OUTPUT_FlatMesh", valid);
	assert(valid);	
	SampledSurface * surface = GetSurfaceByName(fromMeshName);
	assert(surface!=NULL);
	// try to reuse flat mesh file
	VKString reuseFlatMesh = m_Params.GetStrValue(descriptor, "INPUT_ReuseFlat", valid);
	if (valid && reuseFlatMesh!="none")
	{
		SurfaceMap * anotherFlatMesh = surface->GetMapActingAsFrom(reuseFlatMesh);
		assert(anotherFlatMesh!=NULL);
		assert(anotherFlatMesh->GetSurfaceMapType()=="MapFlatMidEdge");
		flatMesh = ((MapFlattening*)anotherFlatMesh)->GetFlatMesh();
		flatSearch = ((MapFlattening*)anotherFlatMesh)->GetKdTree();
	}

	FilenameAddWrkspace(flatMeshFile);
	flatMeshFile.replace("!mshfrom!", fromMeshName);	
	// try to load flat mesh file
	if (!reprocess && flatMesh==NULL && FileExists(flatMeshFile))	
	{
		flatMesh = new R3Mesh();
		flatMesh->ReadFile(flatMeshFile.c_str());
	}
		
	// Create flat map
	if (loadedMap==NULL)	
	{
		if (instance=="MapFlatMidEdge")
		{
			bool reflect = m_Params.GetStrValue(descriptor, "INPUT_Reflect", valid)=="true";
			VKString planarXform = m_Params.GetStrValue(descriptor, "INPUT_PlanarXForm", valid);
			assert(valid);
			VKString cutVertex = m_Params.GetStrValue(descriptor, "INPUT_CutVertex", valid);
			assert(valid);
			
			PlanarTransform * planarTransform = NULL;
			if (planarXform=="MobiusIdentity")
				planarTransform = new MobiusTransformation();
			else if (planarXform=="MVCIdentity")
				planarTransform = new PlanarTransformMVC();
			else if (planarXform=="QuasiConformalIdentity")
				planarTransform = new PlanarTransformQuasiConformal();
			//else if (planarXform=="LSCMIdentity")
			//	planarTransform = new PlanarTransformLSCM();
				
			assert(planarTransform!=NULL);
			if (cutVertex=="any")
				loadedMap = new MapFlatMidEdge(surface, planarTransform, reflect, flatMesh, flatSearch);
			else
			{
				SurfaceFeature * agdFeature = surface->GetFeature(cutVertex);
				assert(agdFeature!=NULL);
				loadedMap = new MapFlatMidEdge(surface, planarTransform, agdFeature, 
											   reflect, flatMesh,  flatSearch);
			}

			assert(loadedMap!=NULL);
			
			// Save Map
			if (mapFilename!="none")
				loadedMap->SaveMap(mapFilename);
			
			// Save flat mesh (if necessary)
			if (flatMesh==NULL)
				((MapFlattening*)loadedMap)->GetFlatMesh()->WriteFile(flatMeshFile.c_str());
		}
		else
		{
			std::cout<<"[ERROR] Unknown flat mesh instance: "<<instance.c_str()<<std::endl;
			assert(false);
		}		
	}
	else	// map loaded (flat mesh should exist)
	{
		assert(loadedMap->GetSurfaceMapType()=="MapFlatMidEdge");
		assert(flatMesh!=NULL);
		((MapFlattening*)loadedMap)->SetFlatMesh(flatMesh, flatSearch);
	}
	
	return loadedMap;
}

void PipelineGeneral::AddConfidences(const VKString & descriptor, bool reprocess)
{
	// add clonable confidence
	AddClonableConfidence(descriptor);
	VKString confidenceName = m_Params.GetStrValue(descriptor, "MEMORY_Name", valid);
	if (!m_ConfidencesNames.contains(confidenceName))
		m_ConfidencesNames.push_back(confidenceName);
	
	VKStringList surfacePairs;
	VKStringList maps;
	
	// TODO fill 3 above
	
	for (int i=0; i<surfacePairs.count(); i+=2)
	{
		VKString surfaceFromStr = surfacePairs[i+0];
		VKString surfaceToStr = surfacePairs[i+1];
		SampledSurface * surfaceFrom = GetSurfaceByName(surfaceFromStr);
		SampledSurface * surfaceTo = GetSurfaceByName(surfaceToStr);
		assert(surfaceFrom!=NULL && surfaceTo!=NULL);
		for (int j=0; j<maps.count(); j++)
		{
			VKString mapName = maps[j];
			
		}
	}
	
	assert(valid);
}

void PipelineGeneral::FilenameAddWrkspace(VKString & name)
{
	VKString wrkspace = m_Params.GetStrValue("Pipeline", "WorkFolder", valid);
	assert(valid);
	name.replace("!wrkdir!", wrkspace);
}

bool PipelineGeneral::FileExists(const VKString & filename)
{
	std::ifstream tempFile(filename.c_str());
	return (tempFile.is_open());
}

VKString PipelineGeneral::GetBasename(const VKString & fullFilename)
{
	int lastSlash = fullFilename.lastIndexOf("/");
	int extension = fullFilename.lastIndexOf(".");
	if (extension==-1)
		extension = fullFilename.length()-1;
	return fullFilename.mid(lastSlash+1, extension-lastSlash-1);
}

VKString PipelineGeneral::GetDir(const VKString & fullFilename)
{
	int lastSlash = fullFilename.lastIndexOf("/");
	if (lastSlash==-1)
		return "";

	return fullFilename.left(lastSlash);
}

SampledSurface * PipelineGeneral::GetSurfaceByName(const VKString & surfaceName)
{
	std::map<VKString, int>::iterator iter = m_mSurfaceNameToID.find(surfaceName);
	if (iter==m_mSurfaceNameToID.end())
		return NULL;
	assert(iter->second < (int)m_Surfaces.size());
	return m_Surfaces[iter->second];
}

//SurfaceMap * PipelineGeneral::GetMapByName(const VKString & mapName)
//{
//	std::map<VKString, int>::iterator iter = m_mMapNameToID.find(mapName);
//	if (iter==m_mMapNameToID.end())
//		return NULL;
//	assert(iter->second < (int)m_SurfaceMaps.size());
//	return m_SurfaceMaps[iter->second];
//}

SurfaceTexture * PipelineGeneral::GetTextureByName(const VKString & textureName)
{
	std::map<VKString, int>::iterator iter = m_mTextureNameToID.find(textureName);
	if (iter==m_mTextureNameToID.end())
		return NULL;
	assert(iter->second < (int)m_Textures.size());
	return m_Textures[iter->second];
}

VKStringList & PipelineGeneral::GetSurfaceNames()
{
	return m_SurfaceNames;
}

VKStringList & PipelineGeneral::GetTextureNames()
{
	return m_TextureNames;
}

VKStringList & PipelineGeneral::GetMapNames()
{
	return m_SurfaceMapNames;
}

VKStringList & PipelineGeneral::GetConfidenceNames()
{
	return m_ConfidencesNames;
}

VKStringList & PipelineGeneral::GetSimilarityNames()
{
	return m_SimilarityNames;
}

VKStringList & PipelineGeneral::GetDistanceMetricNames()
{
	return m_DistanceMetricNames;
}

VKStringList & PipelineGeneral::GetFeatureNames()
{
	return m_FeatureNames;
}

VKStringList & PipelineGeneral::GetSampleSetNames()
{
	return m_SampleSetNames;
}

void PipelineGeneral::Draw(ParamParser & drawParams)
{	
	VKString mapName = drawParams.GetStrValue("RendererDefault", "RenderMap", valid);
	VKString otherMapName = drawParams.GetStrValue("RendererDefault", "RenderOtherMap", valid);
	VKString sampleSetName = drawParams.GetStrValue("RendererDefault", "RenderSampleSet", valid); 
	VKString confidenceName = drawParams.GetStrValue("RendererDefault", "RenderConfidence", valid); 
	VKString similarityName = drawParams.GetStrValue("RendererDefault", "RenderSimilarity", valid);
	
	for (int i=0; i<(int)m_Surfaces.size(); i++)
	{
		SurfaceMap * map1=GetSurface(i)->GetMapActingAsFrom(mapName);
		SurfaceMap * map2=GetSurface(i)->GetMapActingAsFrom(otherMapName);
		if (map1==NULL)
			map1 = GetSurface(i)->GetMapActingAsTo(mapName);
		if (map2==NULL)
			map2 = GetSurface(i)->GetMapActingAsTo(otherMapName);
		
		if (map1!=NULL && map1->GetSurfaceMapType()=="MapScoredCollection")
		{
			MapScoredCollection * mapCollection = (MapScoredCollection*)map1;
			int subMapID = drawParams.GetIntValue("RendererDefault", "RenderSubmap", valid); assert(valid);
			if (subMapID< 0)
			{
				int numVals = vkAbs(subMapID / mapCollection->GetNumMaps()) + 1;
				subMapID += numVals * mapCollection->GetNumMaps();
			}
			assert(subMapID>=0);
			subMapID = subMapID % mapCollection->GetNumMaps();

			mapCollection->SelectMapByID(subMapID);
			//mapCollection->GetMapByID(subMapID)->Draw(m_pWindow, &drawParams);
		}
		if (map1!=NULL)
			map1->Draw(m_pWindow, &drawParams);
		
		GetSurface(i)->SetDrawSurfaceColors(&drawParams, sampleSetName, map1, map2, 
											similarityName, confidenceName);
		GetSurface(i)->Draw(&drawParams);
	}
	DrawSubwindows(drawParams);
	//DrawTextureWindow(drawParams);
}

void PipelineGeneral::DrawSubwindows(ParamParser & drawParams)
{
	glLineWidth(1);
	VKStringList windowTypes;
	windowTypes.push_back("TopRightWindow");
	windowTypes.push_back("TopLeftWindow");
	windowTypes.push_back("BottomRightWindow");
	windowTypes.push_back("BottomLeftWindow");
	
	for (int i=0; i<windowTypes.count(); i++)
	{
		VKString subWindowName = drawParams.GetStrValue("RendererDefault", windowTypes[i], valid);
		if (valid)
		{
			bool top = (windowTypes[i].lastIndexOf("Top")!=-1);
			bool left = (windowTypes[i].lastIndexOf("Left")!=-1);
			m_pWindow->DrawInAdditionalWindowBegin(top, left);
			if (subWindowName=="Textures")
				DrawTextureWindow(drawParams);
			else if (subWindowName=="FlatMapView")
				DrawFlatMapWindow(drawParams);
			m_pWindow->DrawInAdditionalWindowEnd();
		}
	}
}

void PipelineGeneral::DrawTextureWindow(ParamParser & drawParams)
{
	VKString textureName = drawParams.GetStrValue("RendererDefault", "RenderTexture", valid);
	SurfaceTexture * texture = GetTextureByName(textureName);
	if (valid && textureName!="none")
	{
		if (texture==NULL)
			std::cout<<"[WARNING] Texture "<<textureName.c_str()<<" does not exist!"<<std::endl;
		else
			texture->Draw(drawParams);
	}
}

void PipelineGeneral::DrawFlatMapWindow(ParamParser & drawParams)
{
	VKString surfaceName = drawParams.GetStrValue("RendererDefault", "RenderSurface", valid);
	if (valid && surfaceName!="none")
	{
		SampledSurface * surface = GetSurfaceByName(surfaceName);
		if (surface!=NULL)
		{
			VKString mapName = drawParams.GetStrValue("RendererDefault", "RenderMap", valid);
			if (valid && mapName!="none")
			{
				SurfaceMap * surfaceMap = surface->GetMapActingAsFrom(mapName);
				if (surfaceMap!=NULL)
				{
					if (surfaceMap->GetSurfaceMapType().startsWith("MapFlat"))
						((MapFlattening*)surfaceMap)->DrawFlatMesh(drawParams);
				}
				else
					std::cout<<"[WARNING] Cannot find map (DrawFlatMap): "<<mapName.c_str()<<std::endl;
			}
		}
		else
			std::cout<<"[WARNING] Cannot find surface (DrawFlatMap): "<<surfaceName.c_str()<<std::endl;
	}
}

