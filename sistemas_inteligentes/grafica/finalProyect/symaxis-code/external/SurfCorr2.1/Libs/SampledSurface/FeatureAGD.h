#include "SurfaceFeature.h"
#include "VKString.h"

#ifndef __FEATURE_AGD_H
#define __FEATURE_AGD_H

class SampledSurface;
class SurfaceDistance;


class FeatureAGD : public SurfaceFeature
{
	public:
		FeatureAGD(SampledSurface * surface, const VKString & distanceName, 
					bool vertexOnly = true, const VKString & featurename = "AGDFeature");
		virtual VKString GetFeatureName();
	
	protected:
		double m_fNorm;
		virtual double CalculateValue(const SurfaceSample & sample) const;

		SurfaceDistance * m_pDistance;
	
		double * m_aVerticesDistances;
		void * m_aVertexData;	
};

class FeatureFastAGD : public SurfaceFeature
{
	public:
		FeatureFastAGD(SampledSurface * surface, int numInitialSamples,
					   const VKString & featurename = "FastAGDFeature");
		virtual ~FeatureFastAGD();
		
		virtual VKString GetFeatureName();
		virtual void PrecomputeVertexFeatures();

	protected:
		virtual double CalculateValue(const SurfaceSample & sample) const;
		int m_iNumInitialSamples;
};

#endif


