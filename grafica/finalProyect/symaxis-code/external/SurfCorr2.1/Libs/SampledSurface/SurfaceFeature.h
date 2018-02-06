#include "SurfaceSampleSet.h"
#include "SurfacePerVertexValue.h"

class SampledSurface; 

#ifndef __SURFACE_FEATURE_H
#define __SURFACE_FEATURE_H

class SurfaceFeature
{
	public:
		SurfaceFeature(SampledSurface * surface, bool perVertexFeature, VKString name);
		virtual ~SurfaceFeature(){}
	
		virtual double Value(const SurfaceSample & sample);
	
		virtual void Extremas(SurfaceSampleSet * minSet, SurfaceSampleSet * maxSet, double minDistance);
	
		virtual double GetMaximalValue();
		virtual double GetMinimalValue();	
		virtual void PrecomputeVertexFeatures();
	
		virtual void WriteARFFPropertyFile(const VKString & filename);
		virtual void WriteFeature(const VKString & featureName);
		virtual bool LoadFeature(const VKString & featureName);
	
		virtual VKString GetFeatureName() = 0;
	
	protected:
		virtual double CalculateValue(const SurfaceSample & sample) const = 0;
	
		SurfacePerVertexValue * m_pPerVertexValues;	
	
		SampledSurface * m_pSurface;

		VKString m_sFeatureName;
};

#endif


