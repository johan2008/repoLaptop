#include "SurfaceFeature.h"
#include "VKString.h"

#ifndef __FEATURE_CURVATURE_H
#define __FEATURE_CURVATURE_H

class SampledSurface;
class SurfaceDistance;

class FeatureCurvature : public SurfaceFeature
{
	public:
		FeatureCurvature(SampledSurface * surface, int numInitialSamples,
						 const VKString & featurename = "FeatureCurvature");
		virtual ~FeatureCurvature();
		virtual VKString GetFeatureName();
		virtual void PrecomputeVertexFeatures();
	
	protected:
		double m_fNorm;
		virtual double CalculateValue(const SurfaceSample & sample) const;

		int m_iNumInitialSamples;
};

#endif


