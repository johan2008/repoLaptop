#ifndef __PIPELINE_GENERAL_H
#define __PIPELINE_GENERAL_H

#include "AnalysisPipeline.h"
#include "SurfaceTexture.h"

class PipelineGeneral : public AnalysisPipeline
{
	public:
		static PipelineGeneral * s_pGeneralPipeline;
	
		PipelineGeneral(ParamParser & params);
		virtual VKString GetPipelineType();
	
	// Pipeline modules
		virtual bool LoadAndProcessModule(VKString & moduleName, bool reprocess);
		virtual void AddSurfaces(const VKString & descriptor, bool reprocess);
		virtual void AddElementToSurfaces(const VKString & descriptor, bool reprocess, 
										  const VKString & moduleType);
		virtual void AddFeaturesToSurface(const VKString & descriptor, bool reprocess, 
										  const VKString & surfaceName);
		virtual void AddDistancesToSurface(const VKString & descriptor, bool reprocess, 
										   const VKString & surfaceName);
		virtual void AddSamplesToSurface(const VKString & descriptor, bool reprocess, 
										 const VKString & surfaceName);
		virtual void AddTextures(const VKString & descriptor, bool reprocess);
		virtual void AddMaps(const VKString & descriptor, bool reprocess);
		virtual void AddMap(const VKString & descriptor, bool reprocess,
							const VKString & fromMesh, const VKString & toMesh);
		virtual SurfaceMap * AddCollectionOfMaps(const VKString & descriptor, bool reprocess, 
												 const VKString & instance, const VKString & mapFilename,
												 SurfaceMap * loadedMap, const VKString & fromMeshName, 
												 const VKString & toMeshName);
		
		virtual SurfaceMap * AddFlatteningMap(const VKString & descriptor, bool reprocess, 
										  const VKString & instance, const VKString & mapFilename,
										  SurfaceMap * loadedMap, const VKString & fromMesh);
		virtual void AddConfidences(const VKString & descriptor, bool reprocess);
	
	//	getters
		VKStringList & GetSurfaceNames();
		VKStringList & GetTextureNames();
		VKStringList & GetMapNames();
		VKStringList & GetConfidenceNames();
		VKStringList & GetSimilarityNames();
		VKStringList & GetSampleSetNames();
		VKStringList & GetDistanceMetricNames();
		VKStringList & GetFeatureNames();
		SampledSurface * GetSurfaceByName(const VKString & surfaceName);
		SurfaceTexture * GetTextureByName(const VKString & textureName);
	
	// helper functions
		virtual void FilenameAddWrkspace(VKString & name);
		static bool FileExists(const VKString & filename); 
		static VKString GetBasename(const VKString & fullFilename);
		static VKString GetDir(const VKString & fullFilename);
	
	// drawing
		virtual void Draw(ParamParser & drawParams);
		virtual void DrawSubwindows(ParamParser & drawParams);
		virtual void DrawTextureWindow(ParamParser & drawParams);
		virtual void DrawFlatMapWindow(ParamParser & drawParams);
	
	protected:
		std::vector<SurfaceTexture*> m_Textures;
		VKStringList m_TextureNames;
		std::map<VKString, int> m_mTextureNameToID;
	
	// note AnalysiPipeline::m_Surfaces contains surfaces
		VKStringList m_SurfaceNames;
		std::map<VKString, int> m_mSurfaceNameToID;
	
		VKStringList m_SurfaceMapNames;	
		VKStringList m_ConfidencesNames;
		VKStringList m_SimilarityNames;
		VKStringList m_SampleSetNames;
		VKStringList m_DistanceMetricNames;
		VKStringList m_FeatureNames;
};

#endif

