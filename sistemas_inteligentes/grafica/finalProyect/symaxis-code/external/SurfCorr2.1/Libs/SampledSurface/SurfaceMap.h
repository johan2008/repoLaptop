#ifndef __SURFACE_MAP_H
#define __SURFACE_MAP_H

#include <fstream>
#include <sstream>
#include "gaps.h"
#include "VKString.h"
#include "SampledSurface.h"


using namespace std;

class MapCoarse;
class SurfaceMapSimilarity;
class SurfaceMapConfidence;

/**
 * This is the parent class for all maps between surfaces. 
 * Subclasses need to implement two essential functions: ForwardMap and InverseMap
 * In general maps support: saving/loading, similarity (consistency between maps) 
 * and confidence (opposite of distortion), drawing and caching (vertex-to-vertex)
 */
class SurfaceMap 
{
	public:	
		SurfaceMap(SampledSurface * M1, SampledSurface * M2);
		virtual ~SurfaceMap();
		virtual VKString GetSurfaceMapType() = 0;
	
		virtual void InteractiveUpdateFromCorrs(const VKString & sampleSet,
												const VKString & textureName);
	
	// map related: map, inverse, domain and range functions
		virtual SurfaceSample ForwardMap(const SurfaceSample & s) = 0;
		virtual SurfaceSample InverseMap(const SurfaceSample & s) = 0;

		virtual bool ValidForward(const SurfaceSample & s) const;
		virtual bool ValidInverse(const SurfaceSample & s) const;
		virtual const SurfaceSampleSet * GetValidDomain() const;	// return NULL if smooth (i.e. no discrete set)
		virtual const SurfaceSampleSet * GetValidRange() const;		// return NULL if smooth (i.e. no discrete set)

		SampledSurface * GetSurface(int id);
		SampledSurface * GetSurface(const SurfaceSample & samp);
		SampledSurface * GetOtherSurface(const SurfaceSample & samp);
	
	// acccess additional data
		virtual void SetAdditionalData(const VKString & additionalData);
		virtual void AppendAdditionalData(const VKString & appendText);
		virtual VKString GetAdditionalData();
	
	// confidence and similarity
		virtual void ClearAllConfidences();
		virtual void ClearAllSimilarities();
		virtual void ClearConfidenceCalculator(const VKString & mapValueName);
		virtual void ClearMapSimilarity(const VKString & mapValueName, SurfaceMap * otherMap);
		virtual void ClearAllSimilarities(const VKString & mapValueName);
		
		virtual SurfaceMapConfidence * GetMapConfidenceCalculator(const VKString & mapValueName);
		virtual SurfaceMapSimilarity * GetSurfaceMapSimilarity(const VKString & mapValueName, 
															   SurfaceMap * otherMap);
	
	// Drawing routines
		virtual void Draw(AnalysisWindow * window, ParamParser * params,
						  const VKString & renderingParams="RendererDefault", 
						  const VKString & surfaceName="none");
	
		
		virtual void DrawCorrespondenceForSelected(ParamParser * params,
												   const VKString & renderingParams="RendererDefault", 
												   const VKString & surfaceName="none",
												   bool bidir = false);
		virtual void DrawCorrespondenceForSelected(ParamParser * params, const VKString & renderingParams,
												   const VKString & surfaceName,
												   SampledSurface * surfFrom, SampledSurface * surfTo,
												   bool bidir=false);
		virtual SurfaceSample DrawCorrespondenceForSample(ParamParser * params, 
														  const VKString & renderingParams,
														  const VKString & surfaceName,
														  const SurfaceSample & drawCorr);
	
	
	// Loading / saving maps
		virtual void SaveMap(const VKString & outFilename);
		virtual bool LoadMap(const VKString & outFilename);
		static SurfaceMap * CreateMap(const VKString & outFilename, 
									  SampledSurface * surf1, SampledSurface * surf2);
	
		virtual void SaveMap(std::ofstream & textStream) = 0;
		virtual void LoadMap(std::ifstream & textStream) = 0;
		static SurfaceMap * CreateMap(std::ifstream & textStream,
									  SampledSurface * surf1, SampledSurface * surf2);
		
	// Cache a corresponding coarse map (optionally used for ad-hoc speedup)
	//		e.g. map score calculation involves finding a coarse map, that can be used in the future
		void SetCorrespondingCoarseMap(MapCoarse * coarseMap,
									   const VKString & set1, const VKString & set2,
									   const VKString & dist1, const VKString & dist2,
									   int method);
		MapCoarse * GetCurrentCoarseMap();
		MapCoarse * GetCurrentCoarseMap(const VKString & set1, const VKString & set2,
										const VKString & dist1, const VKString & dist2,
										int method, bool & correctCache);
	
		MapCoarse * GetOrCreateCoarseMap(const VKString & set1, const VKString & set2,
										 const VKString & dist1, const VKString & dist2,
										 int method);
	
	// Cache correspondences
		void CacheVertexToVertexMap(bool forwardOnly = true);
		bool IsCachedVertexToVertex();
		void InitializeCache(int maxCacheSize=-1);		
		void EnableCache();
		void DisableCache();	
		bool IsCacheEnabled();
	
		void LockCacheForSave();
		void UnlockCacheForSave();
		bool IsCacheLockedForSave();
			
	protected:
		SampledSurface * m_pM1;
		SampledSurface * m_pM2;
		VKString m_AdditionalData;
	
		std::map<VKString, SurfaceMapConfidence *> m_mConfidences;
		std::map<VKString, std::map<SurfaceMap *, SurfaceMapSimilarity *> > m_mSimilarities;
	
	
	
	// Cached correspondences - can optionally be used by subclasses
		void AddSampleToForwMap(const SurfaceSample & sFrom, const SurfaceSample & sTo);
		void AddSampleToBackMap(const SurfaceSample & sFrom, const SurfaceSample & sTo);
	
		SurfaceSample GetCachedForw(const SurfaceSample & samp, bool & exists);
		SurfaceSample GetCachedBack(const SurfaceSample & samp, bool & exists);	
		
		void WriteCache(std::ofstream & textStream);
		void ReadCache(std::ifstream & textStream);	
	
		int GetCacheSize();
		void ClearCache();
	
	private:
		std::map<int, SurfaceSample> m_ForwCache;
		std::map<int, SurfaceSample> m_BackCache;
		std::vector<int> m_CacheOrder;
		std::vector<bool> m_CacheOrderForwBack;
		int m_iMaxCacheSize;
		bool m_bCacheEnabled;
		bool m_bCacheLockedForSave;

	// Cached corresponding coarse map
		MapCoarse * m_pCorrespondingCoarseMap;
		VKString m_CachedSetName1, m_CachedSetName2;
		VKString m_CachedDistanceName1, m_CachedDistanceName2;	
		int m_CachedMethod;
	
};

#endif
