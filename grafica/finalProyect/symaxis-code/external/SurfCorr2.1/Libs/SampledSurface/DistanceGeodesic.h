#include <fstream>
#include <sstream>

#include "gaps.h"

#include "VKString.h"

#include <vector>
#include <map>
#include <iostream>
#include <assert.h>
#include <fstream>
#include <sstream>


#include "SurfaceDistance.h"
#include "SurfaceSample.h"

using namespace std;

class SampledSurface;
class SurfaceSampleSet;

#ifndef __DISTANCE_GEODESIC_H
#define __DISTANCE_GEODESIC_H

namespace geodesic
{
	class GeodesicAlgorithmBase;
	class Mesh;
}


class DistanceGeodesic : public SurfaceDistance 
{	
	public:
		enum CalculationType
		{
			GEODIST_DIJKSTRA_FUNKHOUSER, 
			GEODIST_DIJKSTRA_KIRSANOV1,
			GEODIST_DIJKSTRA_KIRSANOV2,
			GEODIST_EXACT_KIRSANOV,
			GEODIST_SUBDIVISION_KIRSANOV
		};
	
		struct GeoVertexData 
		{
			R3MeshVertex *vertex;
			double distance;
			RNBoolean selected;
			GeoVertexData **heappointer;
		};
	
	
		DistanceGeodesic(SampledSurface * surface, CalculationType geodesicType=GEODIST_DIJKSTRA_FUNKHOUSER);
		virtual ~DistanceGeodesic(){}
		virtual void PrecomputeDistances();
	
	
		static CalculationType GetTypeFromStr(const VKString & str);

	/** NOTE if exists value passed - no new GeoVertexData is created*/ 
		static GeoVertexData * InitializeFunkPrecomputation(R3Mesh * mesh, bool verb, GeoVertexData * exists=NULL);	
		static void PrecomputeFillFunkhouser(GeoVertexData * tempVertexData, R3Mesh * mesh, int fromVertex, 
											 double * toDistanceArray, 
											 const std::vector<int> & toVerticesRemember=std::vector<int>());
		static void FillThresholdedDistances(GeoVertexData * tempVertexData, R3Mesh * mesh, int fromVertex,
											 std::map<int, double> & verticesWithValues, double maxDistance);
		static double GetThresholdedDistances(int toVertex, R3Mesh * mesh, 
											  std::map<int, double> & verticesWithValues, double maxDistance);
	
	protected:	
		virtual void PrecomputeFunkhouser();
	
		virtual void InitializeKirsanovPrecomputation(geodesic::Mesh ** mesh, 
													  geodesic::GeodesicAlgorithmBase ** algorithm);
		virtual void PrecomputeKirsanov();
		virtual void PrecomputeFillRowKirsanov(int fromVertex, 
											   geodesic::GeodesicAlgorithmBase * algorithm, 
											   geodesic::Mesh * mesh, double * toDistanceArray, 
											   const std::vector<int> & toVerticesRemember=std::vector<int>());
		virtual void FreeKirsanovPrecomputation(geodesic::Mesh * mesh, 
												geodesic::GeodesicAlgorithmBase * algorithm);

		CalculationType m_GeodesicAlgorithm;
		bool m_bVerbous;
};

#endif
