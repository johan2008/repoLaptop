#ifndef __MapMultiConformal_H
#define __MapMultiConformal_H

#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "gaps.h"

#include "VKString.h"
#include "LinAlgMatrixReal.h"

#include "SampledSurface.h"
#include "SurfaceMap.h"
#include "MapScoredCollection.h"
#include "SurfaceDistance.h"
#include "DistanceGeodesic.h"
#include "DistanceOnTheFly.h"
#include "MapConformal.h"

using namespace std;

class MapConfidenceMultiConf;

/**
 * Given a set of (consistent) conformal maps - generate a smooth map over the surface
 */
class MapMultiConformal : public SurfaceMap
{
	public:	
		friend class MapConfidenceMultiConf;
	
		enum InterpolationMethod
		{
			MULTICONF_GEODESIC_CENTROID,
			MULTICONF_AVERAGE_MOBIUS,
			MULTICONF_TAKE_BEST
		};
	
		MapMultiConformal(SampledSurface * surf1, SampledSurface * surf2);
		
		MapMultiConformal(std::vector<MapConformal*> & maps,
						  std::vector< std::map<int, double> > & weights,
						  const VKString & multiConfConfidence,
						  const VKString & interpMethodStr);
	
		MapMultiConformal(MapConformal * confMap, const VKString & multiConfConfidence);
		MapMultiConformal(MapScoredCollection * scoredCollection, 
						  const std::vector<int> & cluster, 
						  const VKString & multiConfConfidence,
						  const VKString & interpMethodStr);	// create multi-map

		MapMultiConformal(MapCoarse * goodCorrs, const VKString & setName,
						  const VKString & multiConfConfidence,
						  const VKString & interpMethodStr, 
						  int MAX_CONF_MAPS);	// interpolate corrs

		MapMultiConformal(MapCoarse * goodCorrs, const VKString & setName,
						  const VKString & multiConfConfidence,
						  const VKString & interpMethodStr, 
						  std::vector<int> & triplets, std::vector<double> & weights);	
				// interpolate with prescribed triplets
	
		void SetGeodesicCentroidOutlierThreshold(double threshold);
		std::vector<int> GetValidCorrs(MapCoarse * goodCorrs, const VKString & setName);
	
		virtual ~MapMultiConformal();
	
		virtual void LocallyRefineCorrs(const VKString & confidenceName, 
										const VKString & searchInSet);
		virtual void PruneConfMaps(const VKString & confidenceName, int maxMaps);	
	
		virtual void SetInterpolationDetails(InterpolationMethod interpolationMethod);
		virtual void SetInterpolationDetails(const VKString & interpolationMethodStr);
		static InterpolationMethod StrToInterMethod(const VKString & interpolationMethod);
	
		MapCoarse * GetCoarseMapForSamples();
	
		virtual SurfaceSample ForwardMap(const SurfaceSample & s);
		virtual SurfaceSample InverseMap(const SurfaceSample & s);
		
		virtual VKString GetSurfaceMapType();
	
		virtual void SaveMap(std::ofstream & outFilename);
		virtual void LoadMap(std::ifstream & outFilename);
	
		virtual void Draw(AnalysisWindow * window, ParamParser * params,
						  const VKString & renderingParams="RendererDefault", 
						  const VKString & surfaceName="none");

		virtual void DrawMultiCorrespondenceForSelected(AnalysisWindow * window, ParamParser * params,
														const VKString & renderingParams,
														const VKString & surfaceName);
	
		virtual void DrawMultiCorrespondenceForSelected(AnalysisWindow * window, ParamParser * params,
														const VKString & renderingParams, 
														const VKString & surfaceName,
														SampledSurface * surfMapFrom, 
														SampledSurface * surfMapTo);
	
		virtual void FillPerVertexWeights(int confMapID, const VKString & weightFlags, 
										  SampledSurface * surface,
										  std::vector<double> & perVertexWeights);

		virtual void FindOptimalWeights();	
		virtual void SmoothOptimalWeights();
		virtual int GetNumConformalMaps();
	
	protected:
		void Initialize(std::vector<MapConformal *> & confMaps, 
						const VKString & multiConfConfidence,
						InterpolationMethod interpMethodStr);
	
		SurfaceSample MultiConfMap(const SurfaceSample & s, bool forward);
		SurfaceSample MultiConfMapGeodesicCentroid(const SurfaceSample & s, bool forward);
		SurfaceSample MultiConfMapAverageMobius(const SurfaceSample & s, bool forward);
		
		double Weight(int mapID, const SurfaceSample & samp);
		double GetSampleConsistency(int confID, const SurfaceSample & samp);
		double GetSampleConfidence(int conformalID, int confidenceID, const SurfaceSample & samp);
	
		// at the end we will probably only need a set of conformal maps
		std::vector<MapConformal *> m_ConformalMaps;
		std::vector<double> m_MapWeight;
	
		InterpolationMethod m_InterpolationMethod;
		VKString m_MultiConformalConfidenceName;

		MapCoarse * m_pCoarseForGenerators;
		
//		bool m_bPreparedWeights;
	
		double UpdateAndGetValueOfConfMaps(std::vector<int> & updateConfMaps, 
										   SurfaceSampleSet * generatingSet2,
										   int s2, const SurfaceSample & newS2,
										   const VKString &confidenceName);
		void GetNeighborhoodInSet(SurfaceSampleSet * searchInSet, 
								  const SurfaceSample & nearMe,
								  std::vector<int> & searchedIDs);
		
	
		bool m_bModifiedGenCorrSet;
		VKString m_sModifiedGenCorrSetName;
		static int s_iGlobalGeneratingCorrSetUniqueIndex;	
		double m_fGeodesicCentroidOutlierThreshold;
};

#endif




