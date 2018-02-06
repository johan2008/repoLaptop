#include "MapScoredCollection.h"
#include "SurfaceMapConfidence.h"
#include "Sortable.h"
#include <algorithm>

MapScoredCollection::MapScoredCollection(SampledSurface * surface1, SampledSurface * surface2)
: SurfaceMap(surface1, surface2)
{
	m_iSelectedMap = 0;
	m_iLastEvaluatedMap = 0;
}

int MapScoredCollection::AddMap(SurfaceMap * map, double score)
{
	m_Maps.push_back(map);
	m_Scores.push_back(score);
	return (int)m_Scores.size()-1;
}

int MapScoredCollection::AddMap(SurfaceMap * currMap, const VKString & confidenceName)
{	
	SurfaceMapConfidence * confidence = currMap->GetMapConfidenceCalculator(confidenceName);
	return AddMap(currMap, confidence->Confidence());
}

void MapScoredCollection::ClearAndDeleteMap(int id)
{
	assert(id < (int)m_Scores.size() && id < (int)m_Maps.size());
	delete m_Maps[id];
	m_Scores.erase(m_Scores.begin()+id);	
	m_Maps.erase(m_Maps.begin()+id);
}

void MapScoredCollection::RecalculateScores(const VKString & confidenceName)
{
	for (int i=0; i<(int)m_Maps.size(); i++)
	{
//		m_Maps[i]->ClearConfidenceCalculator(confidenceName);
		SurfaceMapConfidence * confidence = m_Maps[i]->GetMapConfidenceCalculator(confidenceName);		
		confidence->ClearConfidenceValue();
		m_Scores[i] = confidence->Confidence();
	}
}

int MapScoredCollection::GetNumMaps()
{
	return (int)m_Maps.size();
}

void MapScoredCollection::SortMaps(std::map<int, int> * permutation)
{
	std::vector<Sortable> allMaps;
	for (int i=0; i<GetNumMaps(); i++)
		allMaps.push_back(Sortable(m_Scores[i], m_Maps[i], i));
	std::sort(allMaps.begin(), allMaps.end());
	
	m_Maps.clear();
	m_Scores.clear();
	for (int i=(int)allMaps.size()-1; i>=0; i--)
	{
		if (permutation!=NULL)
			(*permutation)[allMaps[i].id] = m_Maps.size();
		m_Maps.push_back((SurfaceMap*)allMaps[i].ptr);
		m_Scores.push_back(allMaps[i].value);
	}
}

SurfaceMap * MapScoredCollection::GetBestMap(int * id)
{
	int bestID = -1;
	for (int i=0; i<(int)m_Scores.size(); i++)
	{
		if (bestID==-1 || m_Scores[i] > m_Scores[bestID])
			bestID = i;
	}
	if (id!=NULL)
		*id = bestID;
	
	if (bestID!=-1)
		return m_Maps[bestID];
	else
		return NULL;
}

MapCoarse * MapScoredCollection::GetBestCoarseMap(int * id)
{
	assert(GetNumMaps()>0);
	SurfaceMap * surfMap = GetBestMap(id);
	if (surfMap->GetSurfaceMapType()=="MapCoarse")
		return (MapCoarse*)surfMap;
	else
		return surfMap->GetCurrentCoarseMap();
}

SurfaceMap * MapScoredCollection::GetMapByID(int id1, double * score1)
{
	if (!(id1>=0 && id1<(int)m_Maps.size() && id1<(int)m_Scores.size()))
	{
		std::cout<<"Getting invalid map id in scored collection. id="<<id1<<std::endl;
		assert(false);
	}
	if (score1!=NULL)
		*score1 = m_Scores[id1];
	return m_Maps[id1];
}

void MapScoredCollection::SelectBestMap()
{
	GetBestMap(&m_iSelectedMap);
}

int MapScoredCollection::GetSelectedMapID()
{
	return m_iSelectedMap;
}

void MapScoredCollection::IncreaseSelectedMapID()
{
	if (++m_iSelectedMap == (int)m_Maps.size())
		m_iSelectedMap = 0;
}

void MapScoredCollection::DecreaseSelectedMapID()
{
	if (--m_iSelectedMap < 0)
		m_iSelectedMap = m_Maps.size()-1;
}

void MapScoredCollection::SelectMapByID(int id)
{
	m_iSelectedMap = id;
	assert(m_iSelectedMap>=0 && m_iSelectedMap<(int)m_Maps.size());	
}

double MapScoredCollection::GetMapValue(int id)
{
	if(id<0 || id >= (int)m_Scores.size())
		return -1;
	else
		return m_Scores[id];
}

VKString MapScoredCollection::GetSurfaceMapType()
{
	return "MapScoredCollection";
}

////////// SEND TO INDIVIDUAL MAPS /////////////
SurfaceSample MapScoredCollection::ForwardMap(const SurfaceSample & s)
{
	assert(m_iSelectedMap>=0 && m_iSelectedMap<(int)m_Maps.size());
	return m_Maps[m_iSelectedMap]->ForwardMap(s);
}

SurfaceSample MapScoredCollection::InverseMap(const SurfaceSample & s)
{
	assert(m_iSelectedMap>=0 && m_iSelectedMap<(int)m_Maps.size());
	return m_Maps[m_iSelectedMap]->InverseMap(s);
}

bool MapScoredCollection::ValidForward(const SurfaceSample & s) const
{
	assert(m_iSelectedMap>=0 && m_iSelectedMap<(int)m_Maps.size());
	return m_Maps[m_iSelectedMap]->ValidForward(s);	
}

bool MapScoredCollection::ValidInverse(const SurfaceSample & s) const
{
	assert(m_iSelectedMap>=0 && m_iSelectedMap<(int)m_Maps.size());
	return m_Maps[m_iSelectedMap]->ValidInverse(s);		
}

const SurfaceSampleSet * MapScoredCollection::GetValidDomain() const
{
	assert(m_iSelectedMap>=0 && m_iSelectedMap<(int)m_Maps.size());
	return m_Maps[m_iSelectedMap]->GetValidDomain();		
}

const SurfaceSampleSet * MapScoredCollection::GetValidRange() const
{
	assert(m_iSelectedMap>=0 && m_iSelectedMap<(int)m_Maps.size());
	return m_Maps[m_iSelectedMap]->GetValidRange();			
}

void MapScoredCollection::Draw(AnalysisWindow * window,
							   ParamParser * params,
							   const VKString & renderingParams,
							   const VKString & surfaceName)
{
	if (m_Maps.size()>0)
	{
		int mapID = m_Maps.size()-1; 
		if (m_iSelectedMap>=0) 
			mapID = m_iSelectedMap;
//		if (m_Maps[mapID]->GetSurfaceMapType()=="MapConformal")
//		{
//			assert(m_pScoreCalculator!=NULL);
//			if (params->GetStrValue(renderingParams, "ConfSurf", valid)=="RenderOriginal"
//				&& m_Maps[mapID]->GetCurrentCoarseMap()==NULL)	// TODO: can do colors as well
//					m_pScoreCalculator->CalculateValue(m_Maps[mapID]);
//		}
		m_Maps[mapID]->Draw(window, params, renderingParams, surfaceName);
	}
}

void MapScoredCollection::SaveMap(std::ofstream & textStream)
{
	textStream<<"Format MapScoredCollection\n";
	textStream<<m_iSelectedMap<<"\n";	
	textStream<<m_iLastEvaluatedMap<<"\n";
	textStream<<m_Scores.size()<<"\n";
	for (int i=0; i<(int)m_Scores.size(); i++)
	{
		textStream<<m_Scores[i]<<"\n";
		((SurfaceMap*)m_Maps[i])->SaveMap(textStream);
	}
}

void MapScoredCollection::LoadMap(std::ifstream & textStream)
{
//	std::cout<<"Loading MapScoredCollection "<<std::endl;
	std::string tempStr;
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="Format");
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="MapScoredCollection");

	textStream>>m_iSelectedMap;		m_iSelectedMap = 0;
	textStream>>m_iLastEvaluatedMap;
	int numScores;
	textStream>>numScores;
	for (int i=0; i<numScores; i++)
	{
//		std::cout<<"\tCreating Map "<<i<<" / "<<numScores<<"\n";
		double val;
		textStream>>val;
		m_Scores.push_back(val);
		m_Maps.push_back(SurfaceMap::CreateMap(textStream, GetSurface(0), GetSurface(1)));
		assert(m_Maps[m_Maps.size()-1]!=NULL);
	}
//	std::cout<<std::endl;
}

MapScoredCollection * MapScoredCollection::CreateConfMapsReflSymmFromPairs(const VKString & genSampleSetName,
																		   SampledSurface * surface, 
																		   int MAX_MAPS)
{
	assert(false);
	return NULL;
}

MapScoredCollection * MapScoredCollection::CreateConfMapsReflSymmFromTriplets(const VKString & genSampleSetName,
																			  SampledSurface * surface,
																			  int MAX_MAPS)
{
	assert(false);
	return NULL;
}

MapScoredCollection * MapScoredCollection::CreateConfMapsReflFromTriplets(const VKString & genSampleSetName,
																		  SampledSurface * surface1,
																		  SampledSurface * surface2,
																		  int MAX_MAPS)
{
	assert(false);
	return NULL;
}

MapScoredCollection * MapScoredCollection::Union(MapScoredCollection * map1, 
												 MapScoredCollection * map2)
{
	assert(false);
	return NULL;
}


