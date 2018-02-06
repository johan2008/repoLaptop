#ifndef __SURFACE_MAP_SIMILARITY_H
#define __SURFACE_MAP_SIMILARITY_H

#include "gaps.h"
#include "VKString.h"
#include "SurfaceMap.h"
#include "MapCoarse.h"
#include "MapGeoFeature.h"
#include "MapConformal.h"
#include "MapScoredCollection.h"
#include "SurfaceDistance.h"
#include "AnalysisStats.h"
#include "DistanceGeodesic.h"
#include "DistanceOnTheFly.h"
#include "SurfacePerVertexValue.h"

/**
 * Compares a pair of maps (higher value means maps are more similar)
 */
class SurfaceMapSimilarity
{
	public:
		SurfaceMapSimilarity(SurfaceMap * M1, SurfaceMap * M2,
							 bool smoothSimilarity, bool smoothDissimilarity,
							 double smoothingSigma,
							 const VKString & sampleSet1, const VKString & sampleSet2,
							 const VKString & distance1, const VKString & distance2,
							 double simThreshold, bool cacheSimilarity, bool cacheDissimilarity);
		virtual ~SurfaceMapSimilarity();
	
		virtual VKString GetSurfaceMapSimilarityType() = 0;
		
		virtual double Similarity();
		virtual double Dissimilarity();
	
		virtual double SimilarityAtSample(const SurfaceSample & samp);
		virtual double DissimilarityAtSample(const SurfaceSample & samp);
		
		virtual double GetSimilarityThreshold(SampledSurface * surface1, SampledSurface * surface2);
	
		virtual void GetPerVertexSimilarity(SampledSurface * surf, 
											std::vector<double> & vals, bool normalize);
		virtual void GetPerVertexDissimilarity(SampledSurface * surf, 
											   std::vector<double> & vals, bool normalize);
	
		virtual SurfaceMapSimilarity * Clone(SurfaceMap * M1, SurfaceMap * M2) = 0;
	
	
		static SurfaceMapSimilarity * GetClonableSeed(const VKString & seedName);
		static void SetClonableSeed(const VKString & seedName, SurfaceMapSimilarity * seed);

		virtual void PrecomputeSmoothSimilarities();
		virtual void PrecomputeSmoothDissimilarities();	
	
	protected:
		static std::map<VKString, SurfaceMapSimilarity*> m_mNameToSeed;	
	
		virtual double CalculateTotalSimilarity();
		virtual double CalculateTotalDissimilarity();
		
		virtual double CalculateSimilarityAtSample(const SurfaceSample & samp) = 0;
		virtual double CalculateDissimilarityAtSample(const SurfaceSample & samp);
	

		virtual void PrecomputeSmoothSimilarities(SampledSurface * surface);	
	
		virtual void PrecomputeSmoothDissimilarities(SampledSurface * surface);	

	
		double m_fSimilarityThreshold;
		SurfaceMap * m_pSurfaceMap1;
		SurfaceMap * m_pSurfaceMap2;	
	
		VKString m_CoarseSet1Name;
		VKString m_CoarseSet2Name;
		VKString m_Distance1Name;
		VKString m_Distance2Name;
	
		double m_fSimilarityValue;
		double m_fDissimlarityValue;
	
		bool m_bSmoothSimilarity;
		bool m_bSmoothDissimilarity;
		double m_fSmoothSigma;
	
		bool m_bCacheSimilarityValues;
		bool m_bCacheDissimilarityValues;
		std::map<SampledSurface *, SurfacePerVertexValue * > m_mSimilarityCache;
		std::map<SampledSurface *, SurfacePerVertexValue * > m_mDissimilarityCache;
};


#endif


