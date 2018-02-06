#ifndef __FEATURE_IRSA_H
#define __FEATURE_IRSA_H

#include "SurfaceFeature.h"

/**
 * This feature values represent a symmetry axis on a surface.
 * it is similar to partial intrinsic symmetry method by Xu et al, but 
 * uses Moebius Transformations to cast better votes.
 */
class FeatureIRSA : public SurfaceFeature
{
	public:
	/** NOTE: mapCollection should reference a MapScoredCollection with a set of 
			MapVia2DPlane maps, each scored by 'symmetry' (e.g. area preservation), 
			with generators such that c1 <-> c2 are flipped under symmetry, and 
			c3 -> c3 is a stationary point
	 */
		FeatureIRSA(SampledSurface * surface, const VKString & mapCollection,
					const VKString & confidenceName, const VKString & distanceMetric,
					const VKString & featureName);
		
		virtual ~FeatureIRSA();
		
		virtual VKString GetFeatureName();
		virtual void PrecomputeVertexFeatures();
		
	protected:	
		virtual double CalculateValue(const SurfaceSample & sample) const;
		VKString m_MapCollectionName;
		VKString m_ConfidenceName;
		VKString m_DistanceMetric;
};

#endif