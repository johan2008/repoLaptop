#ifndef __DistanceOnTheFly_H
#define __DistanceOnTheFly_H

#include "DistanceGeodesic.h"
#include "SurfaceSampleSet.h"

class DistanceOnTheFly : public SurfaceDistance
{
	public:
		DistanceOnTheFly(SampledSurface * surface, DistanceGeodesic::CalculationType method,
						 SurfaceDistance * tryMeFirst, int maxNumCached);
		DistanceOnTheFly(SampledSurface * surface, DistanceGeodesic::CalculationType method,
						 SurfaceDistance * tryMeFirst, int maxNumCached, double maxDistance);

		virtual ~DistanceOnTheFly();

		virtual void WritePrecomputed(VKString outFile);
		virtual bool LoadPrecomputed(VKString inFile);
	
		virtual int GetCacheSize();
	
		virtual void PrecomputeRow(int vertexFrom);
	
		virtual std::map<int, double> & GetNeighborhood(int vertexID);
		
//	protected:
		virtual void PrecomputeDistances();
		virtual double GetPrecomp(int i, int j);
		virtual void SetPrecomp(int i, int j, double val);
	
		DistanceGeodesic::CalculationType m_Method;
		int m_iMaxNumCached;
		SurfaceDistance * m_pTryMeFirst;
		
		double ** m_PrecomputedRow;
		std::vector<int> m_RowAddedInOrder;

		DistanceGeodesic::GeoVertexData * m_pInitializedFunkComp;

		std::map<int, double> * m_PrecomputedSmallNeighborhoood;
		double m_fDistanceThreshold;
	
};


#endif

