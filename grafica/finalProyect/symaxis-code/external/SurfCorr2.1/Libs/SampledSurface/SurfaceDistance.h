#include <fstream>
#include <sstream>

#include "gaps.h"

#include "VKString.h"

#include "SurfaceSample.h"

using namespace std;

class SampledSurface;

#ifndef __SURFACE_DISTANCE_H
#define __SURFACE_DISTANCE_H

class SurfaceDistance 
{
	public:
		SurfaceDistance(SampledSurface * surface, bool vertexToVertexOnly);
		virtual ~SurfaceDistance();

		virtual double Distance(const SurfaceSample &s1, const SurfaceSample &s2);

		virtual void PrecomputeDistances();	
	
		virtual void WritePrecomputed(VKString outFile);
		virtual bool LoadPrecomputed(VKString inFile);
	
		virtual bool Valid(const SurfaceSample &s1, const SurfaceSample &s2) const;

	protected:
		virtual double CalculateDistance(const SurfaceSample & s1, const SurfaceSample & s2) const;
		virtual double GetPrecomp(int i, int j);
		virtual void SetPrecomp(int i, int j, double val);

		virtual double DistanceToVertex(const SurfaceSample & s, int v); 

		double * m_aPrecomputedDistances;
		int m_iPrecompWidth;
		int m_iPrecompHeight;

		bool m_bVertexToVertex;
		SampledSurface * m_pSurface;
};

#endif
