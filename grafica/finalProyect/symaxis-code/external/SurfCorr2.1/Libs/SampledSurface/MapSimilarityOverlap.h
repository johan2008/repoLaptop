#ifndef __MAP_SIMILARITY_OVERLAP_H
#define __MAP_SIMILARITY_OVERLAP_H

#include "SurfaceMapSimilarity.h"

class DistanceOnTheFly;

/**
 * Fraction of points overlapping
 */
class MapSimilarityOverlap : public SurfaceMapSimilarity
{
	public:
		// to compare any maps
		MapSimilarityOverlap(SurfaceMap * M1, SurfaceMap * M2,
							 const VKString & sampleSet1, const VKString & sampleSet2,
							 const VKString & distance1, const VKString & distance2);
		virtual ~MapSimilarityOverlap();
		virtual VKString GetSurfaceMapSimilarityType();
		SurfaceMapSimilarity * Clone(SurfaceMap * M1, SurfaceMap * M2);
	
	protected:
		virtual double CalculateSimilarityAtSample(const SurfaceSample & samp);
};



#endif

