#ifndef __MAP_GEO_FEATURE_H
#define __MAP_GEO_FEATURE_H

#include <fstream>
#include <sstream>

#include "gaps.h"

#include "VKString.h"
#include "VkFunctions.h"

#include "SampledSurface.h"
#include "SurfaceMap.h"

using namespace std;

class MapCoarse;

class MapGeoFeature : public SurfaceMap
{
	public:
		MapGeoFeature(SampledSurface * M1, SampledSurface * M2);
		MapGeoFeature(SampledSurface * M1, SampledSurface * M2, double interpolationEpsilon,
					  MapCoarse * coarseMap, const VKString & distNameM1, const VKString & distNameM2, 
					  SurfaceSampleSet * validDomainSet=NULL, SurfaceSampleSet * validRangeSet=NULL);
		virtual ~MapGeoFeature();
	
		virtual void PrecomputeDistanceFeatures();
		virtual void DeallocatePrecomputed();
		virtual SurfaceSample ForwardMap(const SurfaceSample & s);
		virtual SurfaceSample InverseMap(const SurfaceSample & s);

		virtual bool ValidForward(const SurfaceSample & s) const;
		virtual bool ValidInverse(const SurfaceSample & s) const;
	
		virtual void SaveMap(std::ofstream & textStream);
		virtual void LoadMap(std::ifstream & textStream);
	
		virtual void Draw(AnalysisWindow * window, ParamParser * params,
						  const VKString & renderingParams="RendererDefault", 
						  const VKString & surfaceName="none");
	
	
		virtual VKString GetSurfaceMapType();

	protected:
		void InitializeMap(SampledSurface * M1, SampledSurface * M2, 
						   SurfaceSampleSet * setNameM1, SurfaceSampleSet * setNameM2,
						   const VKString & distNameM1, const VKString & distNameM2, 
						   SurfaceSampleSet * validDomainSet, SurfaceSampleSet * validRangeSet);
	
		static SurfaceSample AssignVertex(SampledSurface * fromSurf, SampledSurface * toSurf, 
					SurfaceSampleSet * knownFromSet, SurfaceSampleSet * knownToSet,
					SurfaceDistance * fromDist, SurfaceDistance * toDist,
					const SurfaceSample & sFrom, SurfaceSampleSet * validSet);

		static SurfaceSample AssignVertexPrecomputed(SampledSurface * fromSurf, SampledSurface * toSurf, 
					SurfaceSampleSet * knownFromSet, SurfaceSampleSet * knownToSet,
					double ** featuresFrom, double ** featuresTo, int vertexFrom, 
					SurfaceSampleSet * validToSet, double interpolationEpsilon);

		double m_fInterpolationEpsilon;
		virtual void FreeDistances(double ** distanceFeatures, int rows);
		virtual double ** AllocateDistances(int rows, int cols);

		SurfaceSampleSet * m_pKnownCorrM1;
		SurfaceSampleSet * m_pKnownCorrM2;
		SurfaceDistance * m_pDistanceM1;
		SurfaceDistance * m_pDistanceM2;
	
		VKString m_DistanceMetricName1;
		VKString m_DistanceMetricName2;
	
		SurfaceSampleSet * m_pValidDomainSet;	
		SurfaceSampleSet * m_pValidRangeSet;

		int m_iM1NVertices;
		int m_iM2NVertices;
		int m_iNKnownSamples;
	// initialized with PrecomputeDistanceFeatures
	// distance: m_M1*[i][j] = distance from vertex i to known sample j on mesh 1
		double ** m_M1DistanceFeatures;
		double ** m_M2DistanceFeatures;
	
		bool m_bVerbous;
		MapCoarse * m_pInputCoarse;

};

#endif
