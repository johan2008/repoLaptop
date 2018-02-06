#ifndef __MAP_SCORED_COLLECTION_H
#define __MAP_SCORED_COLLECTION_H

#include "SurfaceMap.h"
#include "SurfaceSample.h"
#include "SurfaceSampleSet.h"

class MapScoredCollection : public SurfaceMap
{
	public:
		MapScoredCollection(SampledSurface * surface1, SampledSurface * surface2);
		virtual ~MapScoredCollection(){}
	
		int AddMap(SurfaceMap * map, double score);
		int AddMap(SurfaceMap * map, const VKString & confidenceName);	
	
		void ClearAndDeleteMap(int id);

		SurfaceMap * GetBestMap(int * id = NULL);

		MapCoarse * GetBestCoarseMap(int * id = NULL);
	
		SurfaceMap * GetMapByID(int id, double * score = NULL);
		
		int GetNumMaps();
		void SortMaps(std::map<int, int> * permutation = NULL);

		double GetMapValue(int id);
		void SelectBestMap();
		void SelectMapByID(int id);
		int GetSelectedMapID();
		void IncreaseSelectedMapID();
		void DecreaseSelectedMapID();
	
		void RecalculateScores(const VKString & confidenceName);
	
		virtual SurfaceSample ForwardMap(const SurfaceSample & s);
		virtual SurfaceSample InverseMap(const SurfaceSample & s);
		
		virtual bool ValidForward(const SurfaceSample & s) const;
		virtual bool ValidInverse(const SurfaceSample & s) const;
		virtual const SurfaceSampleSet * GetValidDomain() const;	// return NULL if smooth (i.e. no discrete set)
		virtual const SurfaceSampleSet * GetValidRange() const;		// return NULL if smooth (i.e. no discrete set)
	
		virtual VKString GetSurfaceMapType();	
	
		virtual void Draw(AnalysisWindow * window,
						  ParamParser * params,
						  const VKString & renderingParams="RendererDefault", 
						  const VKString & surfaceName="none");
	
		virtual void SaveMap(std::ofstream & textStream);
		virtual void LoadMap(std::ifstream & textStream);
	
		static MapScoredCollection * CreateConfMapsReflSymmFromPairs(const VKString & genSampleSetName,
																	 SampledSurface * surface, 
																	 int MAX_MAPS);	
		static MapScoredCollection * CreateConfMapsReflSymmFromTriplets(const VKString & genSampleSetName,
																		SampledSurface * surface,
																		int MAX_MAPS);
		static MapScoredCollection * CreateConfMapsReflFromTriplets(const VKString & genSampleSetName,
																	SampledSurface * surface1,
																	SampledSurface * surface2,
																	int MAX_MAPS);
		static MapScoredCollection * Union(MapScoredCollection * map1, MapScoredCollection * map2);
	
	protected:
		int m_iSelectedMap;	
		int m_iLastEvaluatedMap;
		std::vector<double> m_Scores;
		std::vector<SurfaceMap * > m_Maps;
};

#endif

