#ifndef __MAP_CONFIDENCE_COMPARE_TO_TRUTH_H
#define __MAP_CONFIDENCE_COMPARE_TO_TRUTH_H

#include "SurfaceMapConfidence.h"

/**
 * Compare based on truth
 */
class MapConfidenceCompareToTruth : public SurfaceMapConfidence
{
	public:
		MapConfidenceCompareToTruth(SurfaceMap * mapToAnalyze, SurfaceMap * trueMap, 
									const VKString & metric1Name, const VKString & metric2Name, 
									const VKString & set1ToVerify, const VKString & set2ToVerify, 
									double unadjustedThreshold, bool bidir=true);	// arbitrary set
	
		virtual VKString GetMapConfidenceType();
		virtual bool DisplayProgressBar();
		virtual double ErrorThreshold();	
	protected:
		virtual SurfaceMapConfidence * Clone(SurfaceMap * currMap);

		virtual double CalculateConfidenceAtSample(const SurfaceSample & sample);
		virtual double CalculateErrorAtSample(const SurfaceSample & sample);
	
		std::map<SampledSurface *, double *> m_CachedError;
	
		VKString m_sDistName1;
		VKString m_sDistName2;
		SurfaceMap * m_pTrueMap;
	
		double m_fUnadjustedThreshold;
};

#endif


