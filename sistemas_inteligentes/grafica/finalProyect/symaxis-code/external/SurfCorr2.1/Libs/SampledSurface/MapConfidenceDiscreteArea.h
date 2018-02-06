#ifndef __MAP_CONFIDENCE_DISCRETE_AREA_H
#define __MAP_CONFIDENCE_DISCRETE_AREA_H

#include "SurfaceMapConfidence.h"

/**
 * A better(?) approximation to estimate area preservation
 * For each point neighborhood, number of even samples is compared
 * (mapped and samples on surface) must be approximately same
 */
class MapConfidenceDiscreteArea : public SurfaceMapConfidence
{
	public:	
		MapConfidenceDiscreteArea(SurfaceMap * corrMap, double neighborhoodRadius, 
								  const VKString & sampleSet1, const VKString & sampleSet2, 
								  const VKString & dist1="default", const VKString & dist2="default");
		virtual ~MapConfidenceDiscreteArea();
	
		virtual VKString GetMapValueType();
			
	protected:
		virtual SurfaceMapConfidence * Clone(SurfaceMap * currMap);
	
		double CalculateConfidenceAtSample(const SurfaceSample & sample);	
		double m_fNeighborhoodRadius;
		VKString m_Distance1Name;
		VKString m_Distance2Name;	
};


#endif


