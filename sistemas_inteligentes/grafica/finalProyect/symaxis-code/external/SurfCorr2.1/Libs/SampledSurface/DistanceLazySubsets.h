#ifndef __DistanceLazySubsets_H
#define __DistanceLazySubsets_H

#include "DistanceGeodesic.h"
#include "SurfaceSampleSet.h"

class DistanceLazySubsets : public DistanceGeodesic
{
	public:
		DistanceLazySubsets(SampledSurface * surface, 
							SurfaceSampleSet * fromSet, SurfaceSampleSet * toSet,
							DistanceGeodesic::CalculationType method);
		virtual ~DistanceLazySubsets();
	
		virtual bool Valid(const SurfaceSample & s1, const SurfaceSample & s2) const;
	
	protected:
		virtual void PrecomputeFunkhouser();	
		virtual void PrecomputeKirsanov();	
		virtual double GetPrecomp(int i, int j);
		virtual void SetPrecomp(int i, int j, double val);

		int * m_MapVertexFrom;	
		int * m_MapVertexTo;
	
	private:
		bool DistanceExists(int i, int j) const;
};


#endif

