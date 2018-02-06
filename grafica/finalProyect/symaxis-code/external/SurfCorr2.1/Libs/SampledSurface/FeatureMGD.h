#include "SurfaceFeature.h"
#include "VKString.h"

#ifndef __FEATURE_MGD_H
#define __FEATURE_MGD_H

class SampledSurface;
class SurfaceDistance;
class SurfaceSampleSet;

class FeatureMGD : public SurfaceFeature
{
	public:
		FeatureMGD(SampledSurface * surface, const VKString & distanceName, 
				   const VKString & initSetName,
				   bool vertexOnly = true, const VKString & featurename = "MGDFeature");
		FeatureMGD(SampledSurface * surface, const VKString & distanceName, 
				   SurfaceSampleSet * initSampleSet,
				   bool vertexOnly = true, const VKString & featurename = "MGDFeature");

		virtual VKString GetFeatureName();
	
		static SurfaceSampleSet * ConstructSymmetryInvariantSet(SampledSurface * surface, int maxIter, 
																double tau, int maxSamples, 
																const VKString & initSetName,
																const VKString & distanceName = "default",
																std::vector<SurfaceSampleSet *> * historySets = NULL,
																std::vector<FeatureMGD *> * historyFeatures = NULL);	
		static void EvenlySpreadVertices(SampledSurface * surface, int numPnts, std::vector<int> & evenVertices);
		static void EvenlySpreadVertices(SampledSurface * surface, int numPnts, std::vector<int> & evenVertices,
										 std::vector<int> & seeds);
		
	protected:
		virtual double CalculateValue(const SurfaceSample & sample) const;

		SurfaceDistance * m_pDistance;
		SurfaceSampleSet * m_pSampleSet;
	
		struct VertexData {
			R3MeshVertex *vertex;
			RNBoolean selected;
			double distance;
			VertexData **heappointer;
		};	

		static R3MeshVertex * FindFurthestVertex(R3Mesh *mesh, const std::vector<R3MeshVertex *>& seeds, 
												 VertexData *vertex_data);
		static std::map<SampledSurface*, std::map<int, std::vector<int> > > m_EvenVertices;
};

#endif


