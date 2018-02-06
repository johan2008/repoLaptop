#include "MapGeoFeature.h"
#include "SurfaceDistance.h"
#include "DistanceOnTheFly.h"
#include "MapCoarse.h"
#include <assert.h>
#include "AnalysisStats.h"

MapGeoFeature::MapGeoFeature(SampledSurface * M1, SampledSurface * M2)
: SurfaceMap(M1, M2)
{
	m_pInputCoarse = NULL;
	m_fInterpolationEpsilon = 0;
	m_DistanceMetricName1 = "";
	m_DistanceMetricName2 = "";
	
	m_pKnownCorrM1 = NULL;
	m_pKnownCorrM2 = NULL;
	m_pDistanceM1 = NULL;
	m_pDistanceM2 = NULL;
	
	m_M1DistanceFeatures = NULL;
	m_M2DistanceFeatures = NULL;
	
	m_pValidDomainSet = NULL;
	m_pValidRangeSet = NULL;
}

MapGeoFeature::MapGeoFeature(SampledSurface * M1, SampledSurface * M2, double interpolationEpsilon,
							 MapCoarse * coarseMap, 
							 const VKString & distNameM1, const VKString & distNameM2, 
							 SurfaceSampleSet * validDomainSet, SurfaceSampleSet * validRangeSet)
: SurfaceMap(M1, M2)
{
	// for now only support vertex-to-vertex
	assert(validDomainSet==NULL || validDomainSet->NumSamples()==M1->GetMesh()->NVertices());	
	assert(validRangeSet==NULL || validRangeSet->NumSamples()==M2->GetMesh()->NVertices());		
	m_fInterpolationEpsilon = interpolationEpsilon;
	SurfaceSampleSet * setM1;
	SurfaceSampleSet * setM2;	
	coarseMap->FillBestCorr(&setM1, &setM2, -1);
	InitializeMap(M1, M2, setM1, setM2, distNameM1, distNameM2, validDomainSet, validRangeSet);
}


void MapGeoFeature::InitializeMap(SampledSurface * M1, SampledSurface * M2, 
								  SurfaceSampleSet * setM1, SurfaceSampleSet * setM2,
								  const VKString & distNameM1, const VKString & distNameM2, 
								  SurfaceSampleSet * validDomainSet, SurfaceSampleSet * validRangeSet)
{
	
	m_pInputCoarse = new MapCoarse(M1, M2, setM1, setM2);	
	std::vector<int> corrs;
	for (int i=0; i<setM1->NumSamples(); i++)
		corrs.push_back(i);
	m_pInputCoarse->SetFinalCorrMap(corrs);
	
	EnableCache();
	m_DistanceMetricName1 = distNameM1;
	m_DistanceMetricName2 = distNameM2;	
	
	m_bVerbous = true;
	m_pKnownCorrM1 = setM1;
	m_pKnownCorrM2 = setM2;
	m_pDistanceM1 = M1->GetDistanceMetric(distNameM1);
	m_pDistanceM2 = M2->GetDistanceMetric(distNameM2);
	
	m_M1DistanceFeatures = NULL;
	m_M2DistanceFeatures = NULL;
	
	if (validDomainSet==NULL)
		m_pValidDomainSet = new SurfaceSampleSet(M1->GetMesh());
	else
		m_pValidDomainSet = validDomainSet;
	
	if (validRangeSet==NULL)
		m_pValidRangeSet = new SurfaceSampleSet(M2->GetMesh());
	else
		m_pValidRangeSet = validRangeSet;
	
	assert(m_pKnownCorrM1->NumSamples()==m_pKnownCorrM2->NumSamples());	
	PrecomputeDistanceFeatures();
}


MapGeoFeature::~MapGeoFeature()
{
	DeallocatePrecomputed();
}

void MapGeoFeature::DeallocatePrecomputed()
{
	FreeDistances(m_M1DistanceFeatures, m_iM1NVertices);
	FreeDistances(m_M2DistanceFeatures, m_iM2NVertices);
	m_M1DistanceFeatures = NULL;
	m_M2DistanceFeatures = NULL;
}

void MapGeoFeature::FreeDistances(double ** distanceFeatures, int rows)
{
	if (distanceFeatures!=NULL)
	{
		for (int i=0; i<rows; i++)
			delete [] distanceFeatures[i];
		delete [] distanceFeatures;
	}
}

double ** MapGeoFeature::AllocateDistances(int rows, int cols)
{
	double ** distanceFeatures = new double * [rows];
//	std::cout<<"distance features: "<<distanceFeatures<<" : "<<rows<<" x "<<cols<<std::endl;
	for (int i=0; i<rows; i++)
		distanceFeatures[i] = new double[cols];
	return distanceFeatures;
}

void MapGeoFeature::PrecomputeDistanceFeatures()
{
	std::cout<<"Preparing GeoFeatureMap (GMDS)"<<std::endl;
	TimeProfiler profiler;
	R3Mesh * mesh1 = m_pM1->GetMesh();
	R3Mesh * mesh2 = m_pM2->GetMesh();	
	//re-allocate if needed
	if (m_M1DistanceFeatures!=NULL 
		&& (m_iM1NVertices!=mesh1->NVertices() || m_iNKnownSamples!=m_pKnownCorrM2->NumSamples()))	
	{
		FreeDistances(m_M1DistanceFeatures, m_iM1NVertices);
		m_M1DistanceFeatures = NULL;
	}

	if (m_M2DistanceFeatures!=NULL 
		&& (m_iM2NVertices!=mesh2->NVertices() || m_iNKnownSamples!=m_pKnownCorrM1->NumSamples()))	
	{
		FreeDistances(m_M2DistanceFeatures, m_iM2NVertices);
		m_M2DistanceFeatures = NULL;
	}

	m_iM1NVertices = m_pValidDomainSet->NumSamples();
	m_iM2NVertices = m_pValidRangeSet->NumSamples();
	m_iNKnownSamples = m_pKnownCorrM1->NumSamples();
	assert(m_pKnownCorrM1->NumSamples()==m_pKnownCorrM2->NumSamples());

	if (m_M1DistanceFeatures==NULL)
		m_M1DistanceFeatures = AllocateDistances(m_iM1NVertices, m_iNKnownSamples);

	if (m_M2DistanceFeatures==NULL)
		m_M2DistanceFeatures = AllocateDistances(m_iM2NVertices, m_iNKnownSamples);
	
	int TOTAL_TIME = 4 * m_iNKnownSamples;
	int currTime=0;
	DistanceOnTheFly * distances1 = GetSurface(0)->GetOnTheFlyDistanceMetric(-1, "default", -1);
	bool atVertex;
	for (int j=0; j<m_iNKnownSamples; j++)
	{
		if (currTime%20==0)
			profiler.WriteProgress("GeoFeatureDistanceCalc", currTime, TOTAL_TIME);
		currTime++;
		distances1->PrecomputeRow(m_pKnownCorrM1->GetSample(j).NearestVertex(&atVertex));
		assert(atVertex);	// if does not hold - do for each vertex of triangle
	}	
	for (int j=0; j<m_iNKnownSamples; j++)
	{
		if (currTime%20==0)
			profiler.WriteProgress("GeoFeatureDistanceCalc", currTime, TOTAL_TIME);
		currTime++;				
		for (int i=0; i<m_iM1NVertices; i++)
		{		
			m_M1DistanceFeatures[i][j] = distances1->Distance(m_pValidDomainSet->GetSample(i), 
															  m_pKnownCorrM1->GetSample(j));
		}
	}
	
	DistanceOnTheFly * distances2 = GetSurface(1)->GetOnTheFlyDistanceMetric(-1, "default", -1);	
	for (int j=0; j<m_iNKnownSamples; j++)
	{
		if (currTime%20==0)
			profiler.WriteProgress("GeoFeatureDistanceCalc", currTime, TOTAL_TIME);
		currTime++;			
		distances2->PrecomputeRow(m_pKnownCorrM2->GetSample(j).NearestVertex(&atVertex));
		assert(atVertex);	// if does not hold - do for each vertex of triangle
	}	

	for (int j=0; j<m_iNKnownSamples; j++)
	{
		if (currTime%20==0)
			profiler.WriteProgress("GeoFeatureDistanceCalc", currTime, TOTAL_TIME);
		currTime++;			
		for (int i=0; i<m_iM2NVertices; i++)		
		{
			m_M2DistanceFeatures[i][j] = distances2->Distance(m_pValidRangeSet->GetSample(i), 
															  m_pKnownCorrM2->GetSample(j));
		}
	}
	profiler.WriteProgress("GeoFeatureDistanceCalc", TOTAL_TIME, TOTAL_TIME);
	std::cout<<std::endl;
}

	
	// s on M1, return s2 on M2
SurfaceSample MapGeoFeature::ForwardMap(const SurfaceSample & s1)
{
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("MapGeoFeatureCache");
	SurfaceSample retVal;	
	if (IsCacheEnabled())
	{
		bool valid;
		retVal = GetCachedForw(s1, valid);
		if (valid)
			return retVal;
	}
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("MapGeoFeatureCache");	
	
//	std::cout<<"Precomputing"<<std::endl;
//	if (m_M1DistanceFeatures==NULL || m_M2DistanceFeatures==NULL)
//		PrecomputeDistanceFeatures();
//		retVal = AssignVertex(m_pM1, m_pM2, m_pKnownCorrM1, m_pKnownCorrM2, 
//							  m_pDistanceM1, m_pDistanceM2, s1, m_pValidRangeSet);
//	else	//TODO: note this is redundant (caching is same as precomputing here)
//
//	std::cout<<"Assign Precomputed"<<std::endl;	
	assert(m_M1DistanceFeatures!=NULL && m_M2DistanceFeatures!=NULL);
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("AssignPrecomputed");
	retVal = AssignVertexPrecomputed(m_pM1, m_pM2, m_pKnownCorrM1, m_pKnownCorrM2,
									 m_M1DistanceFeatures, m_M2DistanceFeatures, 
									 s1.NearestVertex(), m_pValidRangeSet, m_fInterpolationEpsilon);
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("AssignPrecomputed");	
//	std::cout<<"Mapped"<<std::endl;
	
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("MapGeoFeatureCache");	
	if (IsCacheEnabled() && !IsCacheLockedForSave() && !retVal.Invalid())
		AddSampleToForwMap(s1, retVal);
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("MapGeoFeatureCache");		
	return retVal;
}

SurfaceSample MapGeoFeature::InverseMap(const SurfaceSample & s2)
{
	SurfaceSample retVal;
	if (IsCacheEnabled())
	{
		bool valid;
		retVal = GetCachedBack(s2, valid);
		if (valid)
			return retVal;
	}
	
//	if (m_M1DistanceFeatures==NULL || m_M2DistanceFeatures==NULL)
//		retVal = AssignVertex(m_pM2, m_pM1, m_pKnownCorrM2, m_pKnownCorrM1, 
//								  m_pDistanceM2, m_pDistanceM1, s2, m_pValidDomainSet);
//	else	//TODO: note precomputing vertices is redundant - parent map supports caching
	assert(m_M1DistanceFeatures!=NULL && m_M2DistanceFeatures!=NULL);
	retVal = AssignVertexPrecomputed(m_pM2, m_pM1, m_pKnownCorrM2, m_pKnownCorrM2,
										 m_M2DistanceFeatures, m_M1DistanceFeatures, 
										 s2.NearestVertex(), m_pValidDomainSet, m_fInterpolationEpsilon);
	
	if (IsCacheEnabled() && !IsCacheLockedForSave() && !retVal.Invalid())
		AddSampleToBackMap(s2, retVal);
	
	return retVal;	
}

SurfaceSample MapGeoFeature::AssignVertex(SampledSurface * fromSurf, SampledSurface * toSurf, 
										  SurfaceSampleSet * knownFromSet, SurfaceSampleSet * knownToSet,
										  SurfaceDistance * fromDist, SurfaceDistance * toDist,
										  const SurfaceSample & sFrom, SurfaceSampleSet * validToSet)
{
	assert(false);	// only use precomputed versions
	// check if sFrom is same as some known corr
	for (int i=0; i<knownFromSet->NumSamples(); i++)
	{
		if (sFrom==knownFromSet->GetSample(i))
			return knownToSet->GetSample(i);
	}
	R3Mesh * mesh1 = fromSurf->GetMesh();
	R3Mesh * mesh2 = toSurf->GetMesh();

	int retVToID=-1;

	double scale1 = 1.0 / mesh1->AverageRadius();
	double scale2 = 1.0 / mesh2->AverageRadius();

	double minDelta = -1;
	for (int i=0; i<validToSet->NumSamples(); i++)
	{
		double delta=0;
		for (int j=0; j<knownFromSet->NumSamples(); j++)
		{	
				// distance to sample
			double d1 = fromDist->Distance(sFrom, knownFromSet->GetSample(j)) * scale1;	
				// distance to potential correspondence
			double d2 = toDist->Distance(validToSet->GetSample(i), knownToSet->GetSample(j)) * scale2;	
			double d = vkAbs(d1 - d2);
			if (d1>0)	// if not at sample - add delta
				delta+=d/d1;	
			else if (d>0)	// if at sample - add infinity
				delta += RN_INFINITY;
		}
		if (minDelta<0 || delta < minDelta)
		{
			minDelta = delta;
			retVToID = i;
		}		
	}
//	if (error!=NULL)
//		*error = minDelta;
	return validToSet->GetSample(retVToID);
}


SurfaceSample MapGeoFeature::AssignVertexPrecomputed(SampledSurface * fromSurf, SampledSurface * toSurf, 
				SurfaceSampleSet * knownFromSet, SurfaceSampleSet * knownToSet,
				double ** featuresFrom, double ** featuresTo,
				int vertexFrom, SurfaceSampleSet * validToSet, double interpolationEpsilon)
{
//	R3Mesh * mesh1 = fromSurf->GetMesh();
//	R3Mesh * mesh2 = toSurf->GetMesh();
	
	// NOTE: do not use AverageRadius() call - it's slow!
//	double scale1 = 1.0 / mesh1->AverageRadius();
//	double scale2 = 1.0 / mesh2->AverageRadius();
//	std::cout<<"AssignVertexPrecomputed - 1"<<std::endl;
	double scale1 = 1. / sqrt(fromSurf->Area());
	double scale2 = 1. / sqrt(toSurf->Area());	
	int retVToID=-1;
	
	double minDelta = FLT_MAX;
	for (int i=0; i < validToSet->NumSamples(); i++)
	{
		double delta=0;
		for (int j=0; j<knownFromSet->NumSamples(); j++)
		{
			double d0 = featuresFrom[vertexFrom][j] * scale1;
			double d1 = featuresTo[i][j] * scale2;
			double d = vkAbs((d0 - d1));
			if (interpolationEpsilon==0)	// original tom's method
			{
				if (d0 > 0 ) delta += d / d0;
				else if (d>0) delta += RN_INFINITY;
			}
			else
				assert(false);
				//delta += d / (interpolationEpsilon + d0);	// non-interpolatory method
		}
		if (delta < minDelta)
		{
			minDelta = delta;
			retVToID = i;
		}
	}

//	if (error!=NULL)
//		*error = minDelta;

	//return validToSet->GetSample(retVToID);
	return SurfaceSample(retVToID, toSurf->GetMesh());
}


bool MapGeoFeature::ValidForward(const SurfaceSample & s) const
{
	return m_pValidDomainSet->ContainsExact(s);
}

bool MapGeoFeature::ValidInverse(const SurfaceSample & s) const
{
	return m_pValidRangeSet->ContainsExact(s);
}

VKString MapGeoFeature::GetSurfaceMapType()
{
	return "MapGeoFeature";
}

void MapGeoFeature::Draw(AnalysisWindow * window, ParamParser * params,
						 const VKString & renderingParams,
						 const VKString & surfaceName)
{
	DrawCorrespondenceForSelected(params, renderingParams, surfaceName, 
								  GetSurface(0), GetSurface(1));
	if (m_pInputCoarse!=NULL)
		m_pInputCoarse->Draw(window, params, renderingParams, surfaceName);
}

void MapGeoFeature::SaveMap(std::ofstream & textStream)
{
	textStream<<"Format MapGeoFeature\n";
	m_pKnownCorrM1->WriteSet(textStream);
	m_pKnownCorrM2->WriteSet(textStream);
	m_pValidDomainSet->WriteSet(textStream);
	m_pValidRangeSet->WriteSet(textStream);
	textStream<<m_DistanceMetricName1.c_str()<<"\n";
	textStream<<m_DistanceMetricName2.c_str()<<"\n";		
	
	textStream<<"CacheExists "<<IsCacheEnabled()<<"\n";
	if (IsCacheEnabled())
		WriteCache(textStream);	
}

void MapGeoFeature::LoadMap(std::ifstream & textStream)
{
	std::string tempStr;
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="Format");
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="MapGeoFeature");	
	
	m_pKnownCorrM1 = new SurfaceSampleSet();
	m_pKnownCorrM1->LoadSet(textStream, GetSurface(0)->GetMesh());

	m_pKnownCorrM2 = new SurfaceSampleSet();
	m_pKnownCorrM2->LoadSet(textStream, GetSurface(1)->GetMesh());
	
	m_pValidDomainSet = new SurfaceSampleSet();
	m_pValidDomainSet->LoadSet(textStream, GetSurface(0)->GetMesh());

	m_pValidRangeSet = new SurfaceSampleSet();
	m_pValidRangeSet->LoadSet(textStream, GetSurface(1)->GetMesh());
	
	textStream>>tempStr;		m_DistanceMetricName1 = VKString(tempStr.c_str());
	textStream>>tempStr;		m_DistanceMetricName2 = VKString(tempStr.c_str());	

	
	InitializeMap(GetSurface(0), GetSurface(1), m_pKnownCorrM1, m_pKnownCorrM2, 
				  m_DistanceMetricName1, m_DistanceMetricName2, m_pValidDomainSet, m_pValidRangeSet);	
	
	textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "CacheExists")==0);	
	bool cacheSaved;		textStream>>cacheSaved;
	
	if (cacheSaved)
		ReadCache(textStream);	
}






