#ifndef __MAP_CONFIDENCE_DISTANCE_TO_GEN_H
#define __MAP_CONFIDENCE_DISTANCE_TO_GEN_H

#include "SurfaceMapConfidence.h"

/**
 * Compare based on truth
 */
class MapConfidenceDistanceToGen : public SurfaceMapConfidence
{
	public:
		MapConfidenceDistanceToGen(SurfaceMap * mapToAnalyze, 
								   const VKString & metric1Name, const VKString & metric2Name, 
								   const VKString & set1ToVerify, const VKString & set2ToVerify);	// arbitrary set
		virtual ~MapConfidenceDistanceToGen();
	
		virtual VKString GetMapConfidenceType();

		virtual SurfaceMapConfidence * Clone(SurfaceMap * currMap);
	
	protected:
		virtual double CalculateConfidenceAtSample(const SurfaceSample & sample);
	
		VKString m_sDistName1;
		VKString m_sDistName2;
		
		double m_fEpsilon;
	
		std::map<SampledSurface *, double> m_mMaxConfidence;
	
		std::map<SampledSurface *, std::vector<SurfaceSample> > m_mGenerators;
};

#endif


