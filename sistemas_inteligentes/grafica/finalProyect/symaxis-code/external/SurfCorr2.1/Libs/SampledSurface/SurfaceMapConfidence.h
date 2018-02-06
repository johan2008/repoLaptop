#ifndef __SURFACE_MAP_CONFIDENCE_H
#define __SURFACE_MAP_CONFIDENCE_H


#include "VKString.h"
#include "SurfaceMap.h"
#include "MapCoarse.h"
#include "MapGeoFeature.h"
#include "MapConformal.h"
#include "SurfaceSampleSet.h"
#include "SurfaceDistance.h"
#include "AnalysisStats.h"
#include "SurfacePerVertexValue.h"

#include "gaps.h"
#include <vector>

class MapScoredCollection;

/**
 * Class contains map confidence (estimated according to some measure of quality of a map).
 * Confidence is calculated over the whole surfaces (Surface(0), Surface(1)) or at some point 
 *	on either of the surfaces. Confidence value must be in a range [0..1], but this not checked.
 * By default class supports Error (1-Confidence), and Smooth Confidence (evaluated at some samples
 *	with gaussian weighting. 
 */
class SurfaceMapConfidence
{
	public:
		SurfaceMapConfidence(SurfaceMap * correspondingMap, bool smoothConfidence, bool smoothError,
							 double smoothingSigma,
							 const VKString & testSamples1, const VKString & testSamples2,
							 bool cacheConfidence, bool cacheError);
		
		virtual ~SurfaceMapConfidence();
	
		virtual MapCoarse * GetBestSparseCorrs(SurfaceSampleSet * sampleSetFrom, 
											   int numCorrs, SurfaceSampleSet ** bestSparse);
	
		virtual bool DisplayProgressBar();
	
		virtual VKString GetMapConfidenceType() = 0;
		virtual void ClearConfidenceValue();
		virtual double Confidence();
		virtual double Error();
		virtual double ConfidenceAtSampleNoSmoothing(const SurfaceSample & sample);
	
		virtual double ConfidenceAtSample(const SurfaceSample & sample);
		virtual double ErrorAtSample(const SurfaceSample & sample);
	
		virtual void GetPerVertexConfidence(SampledSurface * surf, 
											  std::vector<double> & vals);
		virtual bool GetPerVertexErrors(SampledSurface * surf, 
										  std::vector<double> & vals, bool normalize);

		virtual void GetPerSampleConfidence(SampledSurface * surf, 
											std::vector<double> & vals);
		virtual void GetPerSampleErrors(SampledSurface * surf, 
										std::vector<double> & vals, bool normalize);
		
		virtual SurfaceMapConfidence * Clone(SurfaceMap * currMap) = 0;
		
		static SurfaceMapConfidence * GetClonableSeed(const VKString & seedName);
		static void SetClonableSeed(const VKString & seedName, SurfaceMapConfidence * seed);
	
		void SaveConfidenceValues(std::ofstream & textStream);
		void LoadConfidenceValues(std::ifstream & textStream);
		
//	protected:
		static std::map<VKString, SurfaceMapConfidence*> m_mNameToSeed;
	
		virtual double ErrorThreshold(SampledSurface * surf);
	
		virtual void PrecomputeSmoothConfidence();
		virtual void PrecomputeSmoothConfidence(SampledSurface * surface);

		virtual void PrecomputeSmoothError();
		virtual void PrecomputeSmoothError(SampledSurface * surface);
	
		virtual double CalculateTotalConfidence();
		virtual double CalculateTotalError();	
	
		virtual double CalculateConfidenceAtSample(const SurfaceSample & sample) = 0;
		virtual double CalculateErrorAtSample(const SurfaceSample & sample);
	
		double m_TotalConfidence;
		double m_TotalError;
	
		bool m_bCacheConfidence;
		bool m_bCacheError;

		std::map<SampledSurface *, SurfacePerVertexValue *> m_CachedConfidenceOnSurf;
		std::map<SampledSurface *, SurfacePerVertexValue *> m_CachedErrorOnSurf;
		
		SurfaceMap * m_pSurfaceMap;
	
		VKString m_SampleSet1;
		VKString m_SampleSet2;	
		
		bool m_bSmoothConfidence;
		bool m_bSmoothError;
		double m_SmoothSigma;
};

#endif
