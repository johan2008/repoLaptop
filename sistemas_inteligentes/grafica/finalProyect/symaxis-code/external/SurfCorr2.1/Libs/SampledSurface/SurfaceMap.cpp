#include "SurfaceMap.h"
#include "MapMultiConformal.h"
#include "MapScoredCollection.h"
#include "MapCoarse.h"
#include "MapConformal.h"
#include "MapEuclidean.h"
#include "MapGeoFeature.h"
#include "SurfaceMapConfidence.h"
#include "SurfaceMapSimilarity.h"
#include "MapFlattening.h"

SurfaceMap::SurfaceMap(SampledSurface * M1, SampledSurface * M2)
{
	m_pCorrespondingCoarseMap = NULL;		
	m_iMaxCacheSize = -1;
	m_bCacheEnabled = false;
	m_bCacheLockedForSave = false;	
	m_AdditionalData = "";
	m_pM1 = M1;
	m_pM2 = M2;
}

SurfaceMap::~SurfaceMap()
{
	ClearAllConfidences();
	ClearAllSimilarities();
}

void SurfaceMap::InteractiveUpdateFromCorrs(const VKString & interactiveSet,
											const VKString & textureName)
{
	std::cout<<"[WARNING] Interactive Update is not implemented"<<std::endl;
}

bool SurfaceMap::ValidForward(const SurfaceSample & s) const
{
	return true;
}

bool SurfaceMap::ValidInverse(const SurfaceSample & s) const
{
	return true;
}

const SurfaceSampleSet * SurfaceMap::GetValidDomain() const
{
	return NULL;	
}

const SurfaceSampleSet * SurfaceMap::GetValidRange() const
{
	return NULL;
}

SampledSurface * SurfaceMap::GetSurface(const SurfaceSample & samp)
{
	if (samp.CheckMeshIsSame(m_pM1->GetMesh()))
		return m_pM1;
	else if (samp.CheckMeshIsSame(m_pM2->GetMesh()))
		return m_pM2;
	else
		assert(false);
	return NULL;
}

SampledSurface * SurfaceMap::GetOtherSurface(const SurfaceSample & samp)
{
	if (samp.CheckMeshIsSame(m_pM1->GetMesh()))
		return m_pM2;
	else if (samp.CheckMeshIsSame(m_pM2->GetMesh()))
		return m_pM1;
	else
		assert(false);
	return NULL;
}

SampledSurface * SurfaceMap::GetSurface(int id)
{
	if (id==0)
		return m_pM1;
	else if (id==1)
		return m_pM2;
	else 
		assert(false);
	return NULL;
}

void SurfaceMap::SetAdditionalData(const VKString & additionalData)
{
	m_AdditionalData = additionalData;
}

void SurfaceMap::AppendAdditionalData(const VKString & appendText)
{
	m_AdditionalData += appendText;
}

VKString SurfaceMap::GetAdditionalData()
{
	return m_AdditionalData;
}

void SurfaceMap::ClearConfidenceCalculator(const VKString & mapValueName)
{
	std::map<VKString, SurfaceMapConfidence *>::iterator iter = m_mConfidences.find(mapValueName);
	if (iter!=m_mConfidences.end())
	{
		delete iter->second;
		m_mConfidences.erase(iter);
	}
}

void SurfaceMap::ClearMapSimilarity(const VKString & mapSimilarityName, SurfaceMap * otherMap)
{
	if (m_mSimilarities.find(mapSimilarityName)!=m_mSimilarities.end()
		|| m_mSimilarities[mapSimilarityName].find(otherMap)!=m_mSimilarities[mapSimilarityName].end())
	{
		delete m_mSimilarities[mapSimilarityName][otherMap];
		m_mSimilarities[mapSimilarityName].erase(otherMap);
		otherMap->m_mSimilarities[mapSimilarityName].erase(this);
	}
}

void SurfaceMap::ClearAllConfidences()
{
	VKStringList confToClear;
	std::map<VKString, SurfaceMapConfidence *>::iterator iter;
	for (iter = m_mConfidences.begin(); iter!=m_mConfidences.end(); iter++)
		confToClear.push_back(iter->first);
	for (int i=0; i<confToClear.count(); i++)
		ClearConfidenceCalculator(confToClear[i]);
}

void SurfaceMap::ClearAllSimilarities()
{
	VKStringList simToClear;
	std::map<VKString, std::map<SurfaceMap *, SurfaceMapSimilarity *> >::iterator iter;
	for (iter = m_mSimilarities.begin(); iter != m_mSimilarities.end(); iter++)
		simToClear.push_back(iter->first);
	for (int i=0; i<simToClear.count(); i++)
		ClearAllSimilarities(simToClear[i]);
}

void SurfaceMap::ClearAllSimilarities(const VKString & mapSimName)
{
	if (m_mSimilarities.find(mapSimName)!=m_mSimilarities.end())
	{
		std::vector<SurfaceMap * > otherMaps;
		
		std::map<SurfaceMap *, SurfaceMapSimilarity *>::iterator iter = m_mSimilarities[mapSimName].begin();
		for (; iter!=m_mSimilarities[mapSimName].end(); iter++)
			otherMaps.push_back(iter->first);
		
		for (int i=0; i<(int)otherMaps.size(); i++)
			ClearMapSimilarity(mapSimName, otherMaps[i]);
		
		m_mSimilarities.erase(mapSimName);
	}
}

SurfaceMapConfidence * SurfaceMap::GetMapConfidenceCalculator(const VKString & mapValueName)
{	
	VKString newMapValueName = mapValueName;
	if (GetSurfaceMapType()=="MapScoredCollection")
		newMapValueName=VKString(((MapScoredCollection*)this)->GetSelectedMapID())+"_"+mapValueName;
	
	if (m_mConfidences.find(newMapValueName)==m_mConfidences.end())
	{
		SurfaceMapConfidence * seed = SurfaceMapConfidence::GetClonableSeed(mapValueName);
		if (seed==NULL)
			std::cout<<"[ERROR] Cannot find seed for confidence: "<<mapValueName.c_str()<<std::endl;
		assert(seed!=NULL);
		if (GetSurfaceMapType()=="MapScoredCollection")
		{
			int selected = ((MapScoredCollection*)this)->GetSelectedMapID();
			if (selected==-1)
				selected = 0;
			m_mConfidences[newMapValueName] = seed->Clone(((MapScoredCollection*)this)->GetMapByID(selected));
		}
		else
			m_mConfidences[newMapValueName] = seed->Clone(this);
	}
	return m_mConfidences[newMapValueName];
}

SurfaceMapSimilarity * SurfaceMap::GetSurfaceMapSimilarity(const VKString & mapSimilarityName, 
														   SurfaceMap * otherMap)
{
	VKString newMapSimilarityName = mapSimilarityName;
	if (GetSurfaceMapType()=="MapScoredCollection")
		newMapSimilarityName=VKString(((MapScoredCollection*)this)->GetSelectedMapID())+"_"+mapSimilarityName;
	
	if (m_mSimilarities.find(newMapSimilarityName)==m_mSimilarities.end()
		|| m_mSimilarities[newMapSimilarityName].find(otherMap)==m_mSimilarities[newMapSimilarityName].end())
	{
		SurfaceMapSimilarity * seed = SurfaceMapSimilarity::GetClonableSeed(mapSimilarityName);
		if (seed==NULL)
			std::cout<<"[ERROR] Cannot find seed for similarity: "<<mapSimilarityName.c_str()<<std::endl;
		assert(seed!=NULL);
		
		SurfaceMap * thisPtr = NULL;
		SurfaceMap * otherPtr = NULL;
		if (GetSurfaceMapType()=="MapScoredCollection")
		{
			int selected = ((MapScoredCollection*)this)->GetSelectedMapID();
			if (selected==-1)
				selected = 0;
			thisPtr = ((MapScoredCollection*)this)->GetMapByID(selected);
		}
		else
		{
			thisPtr = this;
		}		
		if (otherMap->GetSurfaceMapType()=="MapScoredCollection")
		{
			int selected = ((MapScoredCollection*)otherMap)->GetSelectedMapID();
			if (selected==-1)
				selected = 0;
			otherPtr = ((MapScoredCollection*)otherMap)->GetMapByID(selected);
		}
		else
		{
			otherPtr = otherMap;
		}		
		m_mSimilarities[newMapSimilarityName][otherPtr] = seed->Clone(thisPtr, otherPtr);
		otherMap->m_mSimilarities[newMapSimilarityName][thisPtr] = m_mSimilarities[newMapSimilarityName][otherPtr];

	}
	
	return m_mSimilarities[newMapSimilarityName][otherMap];
}

void SurfaceMap::Draw(AnalysisWindow * window, ParamParser * params,
					  const VKString & renderingParams,
					  const VKString & surfaceName)
{
	std::cout<<"[WARNING] SurfaceMap draw function invoked - not implemented."<<std::endl;
}

SurfaceSample SurfaceMap::DrawCorrespondenceForSample(ParamParser * params, 
													  const VKString & renderingParams,
													  const VKString & surfaceName,
													  const SurfaceSample & samp1)
{
	SampledSurface * surfFrom = GetSurface(samp1);
	SampledSurface * surfTo = GetOtherSurface(samp1);	
	SurfaceSample samp2;
	if (surfFrom==GetSurface(0))
	{
		samp2 = ForwardMap(samp1);
//		std::cout<<"\ts1=["<<samp1.VertID(0)<<", "<<samp1.VertID(1)<<", "<<samp1.VertID(2)<<"]"<<std::flush;
//		std::cout<<"\tweights=["<<samp1.B(0)<<", "<<samp1.B(1)<<", "<<samp1.B(2)<<"]"<<std::endl;		
//		std::cout<<"\tMap(s1)=["<<samp2.VertID(0)<<", "<<samp2.VertID(1)<<", "<<samp2.VertID(2)<<"]"<<std::endl;
//		std::cout<<"\tweights=["<<samp2.B(0)<<", "<<samp2.B(1)<<", "<<samp2.B(2)<<"]"<<std::endl;				
	}
	else if (surfFrom==GetSurface(1))
	{
		samp2 = InverseMap(samp1);
//		std::cout<<"\ts1=["<<samp1.VertID(0)<<", "<<samp1.VertID(1)<<", "<<samp1.VertID(2)<<"]"<<std::flush;
//		std::cout<<"\tweights=["<<samp1.B(0)<<", "<<samp1.B(1)<<", "<<samp1.B(2)<<"]"<<std::endl;		
//		std::cout<<"\tInvMap(s1)=["<<samp2.VertID(0)<<", "<<samp2.VertID(1)<<", "<<samp2.VertID(2)<<"]"<<std::endl;
//		std::cout<<"\tweights=["<<samp2.B(0)<<", "<<samp2.B(1)<<", "<<samp2.B(2)<<"]"<<std::endl;						
	}
	else
		assert(false);
	
	if (samp2.Invalid())
		return samp2;
	
	R3Point pnt = surfTo->GetDrawablePosition(samp2, params, renderingParams, surfaceName);
	double rad = surfTo->GetStandardRadius(params);
	R3Point org = surfFrom->GetDrawablePosition(samp1, params, renderingParams, surfaceName);
	
	bool valid;
	glLineWidth(1.);
	VKString color = params->GetStrValue(renderingParams, "SingleCorrColor", valid);
	if (valid && color=="none")
		glColor3d(.3, .3, .3);
	
	glDisable(GL_LIGHTING);
	R3Sphere(pnt, rad).Draw();
	glBegin(GL_LINES);
	glVertex3d(org.X(), org.Y(), org.Z());
	glVertex3d(pnt.X(), pnt.Y(), pnt.Z());
	glEnd();
	return samp2;
}

void SurfaceMap::DrawCorrespondenceForSelected(ParamParser * params, const VKString & renderingParams,
											   const VKString & surfaceName,
											   SampledSurface * surfFrom, SampledSurface * ,
											   bool bidir)
{
	int selectedVertex = surfFrom->GetSelectedVertex(params);

	if (selectedVertex==-1)
		return;
	
//	std::cout<<"Mapping: "<<std::endl;	
	SurfaceSample samp1=SurfaceSample(selectedVertex, surfFrom->GetMesh());
	
	const SurfaceSample & samp2=DrawCorrespondenceForSample(params, renderingParams, surfaceName, samp1);
	
	if (bidir && !samp2.Invalid())
	{
//		std::cout<<std::endl;
		DrawCorrespondenceForSample(params, renderingParams, surfaceName, samp2);
	}
}

void SurfaceMap::DrawCorrespondenceForSelected(ParamParser * params, 
											   const VKString & renderingParams,
											   const VKString & surfaceName,
											   bool bidir)
{
	DrawCorrespondenceForSelected(params, renderingParams, surfaceName, 
								  GetSurface(0), GetSurface(1), bidir);
	DrawCorrespondenceForSelected(params, renderingParams, surfaceName, 
								  GetSurface(1), GetSurface(0), bidir);	
}


void SurfaceMap::SaveMap(const VKString & outFilename)
{
	std::ofstream textStream(outFilename.c_str(), std::ios::out);
	if (!textStream.is_open())
		std::cout<<"[ERROR] Failed to open: "<<outFilename.c_str()<<std::endl;
	assert(textStream.is_open());
	SaveMap(textStream);
}

bool SurfaceMap::LoadMap(const VKString & outFilename)
{
	std::ifstream textStream(outFilename.c_str(), std::ios::in);
	if (!textStream.is_open())
		return false;
	LoadMap(textStream);
	return true;
}

SurfaceMap * SurfaceMap::CreateMap(const VKString & outFilename, SampledSurface * surf1, SampledSurface * surf2)
{
	std::ifstream textStream(outFilename.c_str(), std::ios::in);
	SurfaceMap * retVal=NULL;
	if (textStream.is_open())
	{
		retVal = CreateMap(textStream, surf1, surf2);
		textStream.close();
	}
	return retVal;
}

SurfaceMap * SurfaceMap::CreateMap(std::ifstream & textStream, 
								   SampledSurface * surf1, SampledSurface * surf2)
{
	std::streampos pos=textStream.tellg();
	
	std::string tempStr;
	textStream>>tempStr;	
	if (VKString(tempStr.c_str())!="Format")
	{
		std::cout<<"[ERROR] MAP file is supposed to specify format. Invalid string: "<<tempStr.c_str()<<std::endl;
		assert(VKString(tempStr.c_str())=="Format");
	}
	textStream>>tempStr;	VKString mapName(tempStr.c_str());
	textStream.seekg(pos);
	
	SurfaceMap * newMap = NULL;
	
	if (mapName=="MapMultiConformal")
		newMap = new MapMultiConformal(surf1, surf2);
	else if (mapName=="MapScoredCollection")
		newMap = new MapScoredCollection(surf1, surf2);
	else if (mapName=="MapCoarse")
		newMap = new MapCoarse(surf1, surf2);
	else if (mapName=="MapConformal")
	{
		MobiusTransformation m;		
		assert(surf1->GetSurfaceType()=="SurfaceMidEdgeConf");
		assert(surf2->GetSurfaceType()=="SurfaceMidEdgeConf");		
		newMap = new MapConformal((SurfaceMidEdgeConf*)surf1, (SurfaceMidEdgeConf*)surf2, m, m);
	}
	else if (mapName=="MapEuclidean")
		newMap = new MapEuclidean(surf1, surf2);
	else if (mapName=="MapGeoFeature")
		newMap = new MapGeoFeature(surf1, surf2);
	else if (mapName=="MapFlatMidEdge")
		newMap = new MapFlatMidEdge(surf1);
	else if (mapName=="MapVia2DPlane")
		newMap = new MapVia2DPlane(surf1, surf2);
	
	newMap->LoadMap(textStream);	
	
	return newMap;
}


////////////////// CACHING ///////////////////
bool SurfaceMap::IsCachedVertexToVertex()
{
	return (int)m_ForwCache.size()>=GetSurface(0)->GetMesh()->NVertices();
}

void SurfaceMap::CacheVertexToVertexMap(bool forwardOnly)
{
	ClearCache();
	InitializeCache();
	EnableCache();
	
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("FinalVertexToVertex");
	std::cout<<"Caching Vertex To Vertex Map: "<<std::endl;
	int numSurf = (forwardOnly ? 1 : 2);
	int MAX_MAPS = 0;
	for (int surfID=0; surfID < numSurf; surfID++)
		MAX_MAPS+=GetSurface(surfID)->GetMesh()->NVertices();
	
	for (int surfID=0; surfID < numSurf; surfID++)
	{
		R3Mesh * mesh = GetSurface(surfID)->GetMesh();
		std::cout<<"\tSurf"<<surfID<<" ("<<mesh->NVertices()<<"): "<<std::flush;		
		for (int i=0; i<mesh->NVertices(); i++)
		{
//			if (i%500==0)
//				std::cout<<i<<" "<<std::flush;
			if (i%100==0)
				AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("FinalVertexToVertex", 
																	(surfID*mesh->NVertices() + i)/50,
																	(MAX_MAPS/50));
			SurfaceSample samp(i, mesh);
			if (surfID==0)
				ForwardMap(samp);
			else
				InverseMap(samp);
		}
//		std::cout<<std::endl;
	}
	std::cout<<std::endl;
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("FinalVertexToVertex");	
}

void SurfaceMap::InitializeCache(int maxCacheSize)
{
	m_bCacheEnabled = true;
	m_iMaxCacheSize = maxCacheSize;
}

void SurfaceMap::AddSampleToForwMap(const SurfaceSample & sFrom, const SurfaceSample & sTo)
{
	assert(m_bCacheEnabled && !m_bCacheLockedForSave);
	assert(sFrom.CheckMeshIsSame(GetSurface(0)->GetMesh()));
	assert(sTo.CheckMeshIsSame(GetSurface(1)->GetMesh()) || sTo.Invalid());
	
	bool valid;
	int vertexID = sFrom.NearestVertex(&valid);
	if (!valid)
		return;
	
	m_CacheOrder.push_back(vertexID);
	m_CacheOrderForwBack.push_back(true);
	m_ForwCache[vertexID] = sTo;
	
	while (m_iMaxCacheSize>=0 && (int)m_CacheOrder.size() > m_iMaxCacheSize)
	{
		int deleteVertex = m_CacheOrder[0];
		bool forw = m_CacheOrderForwBack[0];
		if (forw)
			m_ForwCache.erase(deleteVertex);
		else
			m_BackCache.erase(deleteVertex);
		m_CacheOrder.erase(m_CacheOrder.begin());
		m_CacheOrderForwBack.erase(m_CacheOrderForwBack.begin());
	}
}

void SurfaceMap::AddSampleToBackMap(const SurfaceSample & sFrom, const SurfaceSample & sTo)
{
	assert(m_bCacheEnabled && !m_bCacheLockedForSave);
	assert(sFrom.CheckMeshIsSame(GetSurface(1)->GetMesh()));
	assert(sTo.CheckMeshIsSame(GetSurface(0)->GetMesh()) || sTo.Invalid());
	
	bool valid;
	int vertexID = sFrom.NearestVertex(&valid);
	if (!valid)
		return;
	
	m_CacheOrder.push_back(vertexID);
	m_CacheOrderForwBack.push_back(false);
	m_BackCache[vertexID] = sTo;
	
	while (m_iMaxCacheSize>=0 && (int)m_CacheOrder.size() > m_iMaxCacheSize)
	{
		int deleteVertex = m_CacheOrder[0];
		bool forw = m_CacheOrderForwBack[0];
		if (forw)
			m_ForwCache.erase(deleteVertex);
		else
			m_BackCache.erase(deleteVertex);
		m_CacheOrder.erase(m_CacheOrder.begin());
		m_CacheOrderForwBack.erase(m_CacheOrderForwBack.begin());
	}	
}

SurfaceSample SurfaceMap::GetCachedForw(const SurfaceSample & sFrom, bool & exists)
{
	exists = false;
	if (sFrom.Invalid())
		return SurfaceSample();
	assert(sFrom.CheckMeshIsSame(GetSurface(0)->GetMesh()));
	int vertexID = sFrom.NearestVertex(&exists);
	if (!exists)
		return SurfaceSample();
	
	std::map<int, SurfaceSample>::iterator iter = m_ForwCache.find(vertexID);
	exists = (iter!=m_ForwCache.end());
	if (exists)
		return iter->second;
	else
		return SurfaceSample();
}

SurfaceSample SurfaceMap::GetCachedBack(const SurfaceSample & sFrom, bool & exists)
{
	exists = false;
	if (sFrom.Invalid())
		return SurfaceSample();

	assert(sFrom.CheckMeshIsSame(GetSurface(1)->GetMesh()));
	
	int vertexID = sFrom.NearestVertex(&exists);
	if (!exists)
		return SurfaceSample();
	
	std::map<int, SurfaceSample>::iterator iter = m_BackCache.find(vertexID);
	exists = (iter!=m_BackCache.end());
	if (exists)
		return iter->second;
	else
		return SurfaceSample();
}

void SurfaceMap::WriteCache(std::ofstream & textStream)
{
	textStream<<"Cache "<<m_CacheOrder.size()<<" "<<m_iMaxCacheSize<<"\n";
	for (int i=0; i<(int)m_CacheOrder.size(); i++)
	{
		int vertexID = m_CacheOrder[i];
		bool forw = m_CacheOrderForwBack[i];
		SurfaceSample samp;
		if (forw)
			samp = m_ForwCache[vertexID];
		else
			samp = m_BackCache[vertexID];
		
		textStream<<vertexID<<" "<<forw<<" "<<samp.TriID()<<" "<<samp.B(0)<<" "<<samp.B(1)<<" "<<samp.B(2)<<"\n";
	}
	
	// save corresponding coarse map (if any)
	if (m_pCorrespondingCoarseMap!=NULL)
	{
		textStream<<"1\n";
		textStream<<m_CachedSetName1.c_str()<<" "<<m_CachedSetName2.c_str()<<"\n";
		textStream<<m_CachedDistanceName1.c_str()<<" "<<m_CachedDistanceName2.c_str()<<"\n";
		textStream<<m_CachedMethod<<"\n";
		m_pCorrespondingCoarseMap->SaveMap(textStream);
	}
	else
		textStream<<"0\n";	
}

void SurfaceMap::ReadCache(std::ifstream & textStream)
{	
	std::string tempStr;
	textStream>>tempStr;		assert(tempStr=="Cache");
	int cacheSize;
	textStream>>cacheSize>>m_iMaxCacheSize;
	
	for (int i=0; i<cacheSize; i++)	
	{
		int vertexID;
		bool forw;
		int triID;
		double b1, b2, b3;
		
		textStream>>vertexID>>forw>>triID>>b1>>b2>>b3;
		
		m_CacheOrder.push_back(vertexID);
		m_CacheOrderForwBack.push_back(forw);
		
		if (forw)
			m_ForwCache[vertexID] = SurfaceSample(triID, b1, b2, b3, GetSurface(1)->GetMesh());
		else
			m_BackCache[vertexID] = SurfaceSample(triID, b1, b2, b3, GetSurface(0)->GetMesh());		
	}
	
	// save corresponding coarse map (if any)
	bool loadCoarse;
	textStream>>loadCoarse;
	if (loadCoarse)
	{
		textStream>>tempStr;
		m_CachedSetName1 = VKString(tempStr.c_str());
		textStream>>tempStr;
		m_CachedSetName2 = VKString(tempStr.c_str());
		textStream>>tempStr;
		m_CachedDistanceName1 = VKString(tempStr.c_str());
		textStream>>tempStr;
		m_CachedDistanceName2 = VKString(tempStr.c_str());
		int cachedInt;
		textStream>>cachedInt;
		m_CachedMethod = cachedInt;
		
		m_pCorrespondingCoarseMap = (MapCoarse*)SurfaceMap::CreateMap(textStream, GetSurface(0), GetSurface(1));
		assert(m_pCorrespondingCoarseMap!=NULL);
	}	
}

int SurfaceMap::GetCacheSize()
{
	return m_iMaxCacheSize;
}

void SurfaceMap::ClearCache()
{
	m_ForwCache.clear();
	m_BackCache.clear();
	m_CacheOrder.clear();
	m_CacheOrderForwBack.clear();
}

void SurfaceMap::EnableCache()
{
	m_bCacheEnabled = true;
}

void SurfaceMap::DisableCache()
{
	m_bCacheEnabled = false;
}

bool SurfaceMap::IsCacheEnabled()
{
	return m_bCacheEnabled;
}

bool SurfaceMap::IsCacheLockedForSave()
{
	return m_bCacheLockedForSave;
}

void SurfaceMap::LockCacheForSave()
{
	m_bCacheLockedForSave = true;
}

void SurfaceMap::UnlockCacheForSave()
{
	m_bCacheLockedForSave = false;
}

MapCoarse * SurfaceMap::GetOrCreateCoarseMap(const VKString & set1Name, const VKString & set2Name,
								 const VKString & dist1, const VKString & dist2,
								 int method)
{		
	bool exists;
	MapCoarse * coarseMap = GetCurrentCoarseMap(set1Name, set2Name, dist1, dist2, 
												(MapCoarse::FineToCoarseGeneration)method, exists);
	if (!exists || coarseMap==NULL)
	{
		SurfaceSampleSet * set1 = GetSurface(0)->GetSampleSet(set1Name);
		SurfaceSampleSet * set2 = GetSurface(1)->GetSampleSet(set2Name);		
		coarseMap = new MapCoarse(set1, set2, this, (MapCoarse::FineToCoarseGeneration)method, 
								  dist1, dist2);
		SetCorrespondingCoarseMap(coarseMap, set1Name, set2Name, dist1, dist2, 
								  (MapCoarse::FineToCoarseGeneration)method);
	}
	
	return coarseMap;
}

void SurfaceMap::SetCorrespondingCoarseMap(MapCoarse * coarseMap,
										   const VKString & set1, const VKString & set2,
										   const VKString & dist1, const VKString & dist2,
										   int method)
{
	// this is to avoid bugs - for now only support one corresponding coarse map
	// need to set corresponding coarse map to NULL before using it.
	assert(m_pCorrespondingCoarseMap==NULL ||coarseMap==NULL);	
	
	m_pCorrespondingCoarseMap = coarseMap;
	m_CachedSetName1 = set1;
	m_CachedSetName2 = set2;
	m_CachedDistanceName1 = dist1;
	m_CachedDistanceName2 = dist2;	
	m_CachedMethod = method;
}
MapCoarse * SurfaceMap::GetCurrentCoarseMap(const VKString & set1, const VKString & set2,
											const VKString & dist1, const VKString & dist2,
											int method, bool & correctCache)
{
	correctCache = (set1==m_CachedSetName1) && (set2==m_CachedSetName2) && (method==m_CachedMethod);
	correctCache = correctCache && (dist1==m_CachedDistanceName1) && (dist2==m_CachedDistanceName2);
	
	return m_pCorrespondingCoarseMap;
}

MapCoarse * SurfaceMap::GetCurrentCoarseMap()
{
	return m_pCorrespondingCoarseMap;
}

