#ifndef __ANALYSIS_PIPELINE_H
#define __ANALYSIS_PIPELINE_H

#include <fstream>
#include <sstream>

#include "gaps.h"

#include "VKString.h"
#include "ParamParser.h"

#include "SampledSurface.h"
#include "SurfaceMap.h"
#include "SurfaceMidEdgeConf.h"
#include "DistanceGeodesic.h"
#include "FeatureAGD.h"
#include "FeatureMGD.h"
#include "MapCoarse.h"
#include "MapConformal.h"
#include "MapEuclidean.h"
#include "MapGeoFeature.h"
#include "SurfaceMap.h"
#include "AnalysisWindow.h"
#include "SampleSetMapGenerator.h"
#include "MapScoredCollection.h"
#include "AnalysisStats.h"
#include "SurfaceMapSimilarity.h"
#include "DistanceLazySubsets.h"
#include "DistanceOnTheFly.h"
#include "SurfaceMapConfidence.h"
#include "MapConfidenceMutuallyClosest.h"
#include "MapConfidenceDiscreteArea.h"
#include "MapConfidenceFaceArea.h"
#include "MapConfidenceCompareToTruth.h"
#include "MapConfidenceMultiConf.h"
#include "MapConfidenceDistanceToGen.h"
#include "SurfaceMapSimilarity.h"
#include "MapSimilarityOverlap.h"
#include "MapSimilarityDistance.h"

using namespace std;

class AnalysisPipeline 
{
	public:	
		static AnalysisPipeline * CreatePipeline(ParamParser & params);
		virtual ~AnalysisPipeline();

		// Some static members
		static bool Verb(const VKString & verbousness );
	
	// Below is used by algorithms
		virtual void AddBenchmark(const VKString & confidenceName, SurfaceMap * predictionMap, 
								  ParamParser * benchParams, int statID, double score, 
								  bool includesSymFlip);	
		virtual void LoadSamplesFromCoarseAlgorithm(const VKString & otherMapPrefix);
		virtual void LoadTruthMap();
		virtual void LoadAndProcessSurfaces(const VKString & algorithmName);
		virtual void LoadAndProcessDistances(SampledSurface * surface, const VKString & surfaceName);
		virtual void LoadAndProcessFeatures(SampledSurface * surface, const VKString & surfaceName,
											SurfaceMidEdgeConf ** needFlattening, 
											std::vector<SampledSurface*> &symmetricSurfaces,
											const VKString & meshName);	
		virtual void LoadAndProcessSampleSets(SampledSurface * surface, const VKString & surfaceName);
		virtual void LoadSetSamplers(const VKString & algorithmName);
		virtual void LoadAndProcessBenchmark();

		virtual int ClustersToAnalyze(); 
		virtual void TryLoadCoarseMap(const VKString & filename, MapCoarse ** saveTo);
		
		virtual void FillSymmetricSurfacesConf(VKStringList & typeSymmetries, 
									   R3Mesh * mesh, R3Mesh * flatMesh,
									   std::vector<SampledSurface*> & fillMe);
		virtual SampledSurface * CreateSurface(const VKString & surfaceDescriptor, 
									   std::vector<SampledSurface*> &symmetricSurfaces, int meshID);
	
		virtual void AddSurfaceCamera(const VKStringList & cameraFiles, 
									  const VKStringList & cameraPositions,
									  const VKStringList & cameraOrientations,
									  const VKStringList & cameraUps,
									  const VKStringList & cameraLines,
									  int meshID);

		virtual SurfaceSampleSet * CreateSampleSet(const VKString & setName, 
												   SampledSurface * surface);
		virtual SurfaceDistance * CreateDistanceMetric(const VKString & metricDescriptor, 
											   SampledSurface * surface);
		virtual SurfaceFeature * CreateFeature(const VKString & featureDescriptor, 
									   SampledSurface * surface);
		virtual SampleSetMapGenerator * CreateSetSampler(const VKString & samplerDescriptor, 
												 const VKString & sampleSetName,
												 std::vector<SampledSurface *>& surfaces);
		virtual void AddClonableConfidence(const VKString & mapEvalDescriptor);
		virtual void AddClonableSimilarity(const VKString & descriptor);
	
		virtual void WriteLog();
	
		virtual int NumSurfaces();
		virtual SampledSurface * GetSurface(int i);
	
		virtual VKString GetPipelineType();
		
	// Rendering
		virtual void SetWindow(AnalysisWindow * window);
		virtual void Draw(ParamParser & drawParams);
		virtual void WriteScreenshot(const R2Image & img, const VKString & filename);
		
		MapScoredCollection * GetCurrCollection();
		
	protected:
		AnalysisPipeline(ParamParser & params);
		R3Mesh * CopyMesh(R3Mesh * mesh);
	
	// collection of meshes
		std::vector<SampledSurface*> m_Surfaces;
	
	// algorithm (some algorithms might use only a subset of below structures): 
		// used by load / create set samplers
		std::vector<SampleSetMapGenerator*> m_CandidateCorrGenerators;
		VKString m_GeneratorSetName1;
		VKString m_GeneratorSetName2;	
	
		// used in load map accumulator / write log
		MapScoredCollection * m_pCollectionOfMaps;
		
		// the final map
		SurfaceMap * m_pTheFinalMap;
	
		
		// used for benchmarking results
		MapCoarse * m_pBenchmarkQuery;	
		MapCoarse * m_pTruthMap;
	// to verify whether there is a symmetry flip
		MapCoarse * m_pAlternativeTruth;
		MapCoarse * m_pAlternativeSymmetricTruth;
	
	// Param parsing:
		virtual VKString LocalFilename(const VKString & meshname, const VKString & ext);
		virtual VKString LocalFilename(const VKString & ext);
		virtual VKString LocalFilename(const VKString & meshname, 
									   const VKString & filename, const VKString & ext);
		ParamParser & m_Params;
		bool valid;
		static VKString s_Verbous;
	
	// Rendering and GUI
		AnalysisWindow * m_pWindow;	
		std::vector<double *> m_InitCameraPositions;
		
		ParamParser * m_pExternalAlgorithmTime;
};

#endif




