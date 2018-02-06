#include "DistanceOnTheFly.h"
#include "SampledSurface.h"

DistanceOnTheFly::DistanceOnTheFly(SampledSurface * surface, DistanceGeodesic::CalculationType method,
								   SurfaceDistance * tryMeFirst, int maxNumCached)
: SurfaceDistance(surface, true)
{
	m_pTryMeFirst = tryMeFirst;
	m_iMaxNumCached = maxNumCached;
	assert(method==DistanceGeodesic::GEODIST_DIJKSTRA_FUNKHOUSER);
	m_Method = method;
	m_aPrecomputedDistances = NULL;
	m_iPrecompWidth = surface->GetMesh()->NVertices();
	m_iPrecompHeight = 0;
	m_pInitializedFunkComp = NULL;
	m_fDistanceThreshold = -1;
	
	m_PrecomputedRow = new double*[m_iPrecompWidth];
	for (int i=0; i<m_iPrecompWidth; i++)
		m_PrecomputedRow[i] = NULL;
}

DistanceOnTheFly::DistanceOnTheFly(SampledSurface * surface, DistanceGeodesic::CalculationType method,
								   SurfaceDistance * tryMeFirst, int maxNumCached, double maxDistance)
: SurfaceDistance(surface, true)
{
	assert(maxNumCached==-1);
	m_pTryMeFirst = tryMeFirst;
	m_iMaxNumCached = maxNumCached;
	assert(method==DistanceGeodesic::GEODIST_DIJKSTRA_FUNKHOUSER);
	m_Method = method;
	m_aPrecomputedDistances = NULL;
	m_iPrecompWidth = surface->GetMesh()->NVertices();
	m_iPrecompHeight = 0;
	m_pInitializedFunkComp = NULL;	
	m_fDistanceThreshold = maxDistance;
	
	m_PrecomputedRow = NULL;
	m_PrecomputedSmallNeighborhoood = new std::map<int, double>[m_iPrecompWidth];
}

DistanceOnTheFly::~DistanceOnTheFly()
{
	if (m_PrecomputedRow!=NULL)
	{
		for (int i=0; i<m_iPrecompWidth; i++)
			if (m_PrecomputedRow[i]!=NULL)
				delete [] m_PrecomputedRow[i];
		delete [] m_PrecomputedRow;
	}
	else
	{
		delete [] m_PrecomputedSmallNeighborhoood;
	}
}

void DistanceOnTheFly::WritePrecomputed(VKString outFile)
{
	return;
}

bool DistanceOnTheFly::LoadPrecomputed(VKString inFile)
{
	return false;
}

void DistanceOnTheFly::PrecomputeDistances()
{
	// don't do anything
}

void DistanceOnTheFly::PrecomputeRow(int vertexFrom)
{
	if (m_fDistanceThreshold==-1 && m_PrecomputedRow[vertexFrom] == NULL)
	{
		m_RowAddedInOrder.push_back(vertexFrom);
		double * rowToFill=NULL;
		
		if (m_iMaxNumCached!=-1 && (int)m_RowAddedInOrder.size()>m_iMaxNumCached)	// cache filled
		{
			rowToFill = m_PrecomputedRow[vertexFrom];
			m_PrecomputedRow[vertexFrom] = NULL;
			m_RowAddedInOrder.erase(m_RowAddedInOrder.begin());
		}	
		else		// cache is not filled
			rowToFill = new double[m_pSurface->GetMesh()->NVertices()];
		
		m_PrecomputedRow[vertexFrom] = rowToFill;	
		
		m_pInitializedFunkComp = DistanceGeodesic::InitializeFunkPrecomputation(m_pSurface->GetMesh(), 
																				false, m_pInitializedFunkComp);	
		
		DistanceGeodesic::PrecomputeFillFunkhouser(m_pInitializedFunkComp, m_pSurface->GetMesh(), vertexFrom, rowToFill);
	}
	else if (m_fDistanceThreshold>=0 && m_PrecomputedSmallNeighborhoood[vertexFrom].size()==0)
	{
		m_pInitializedFunkComp = DistanceGeodesic::InitializeFunkPrecomputation(m_pSurface->GetMesh(), false,
																				m_pInitializedFunkComp);	
		DistanceGeodesic::FillThresholdedDistances(m_pInitializedFunkComp, m_pSurface->GetMesh(), vertexFrom, 
												   m_PrecomputedSmallNeighborhoood[vertexFrom], m_fDistanceThreshold);		
	}
}

std::map<int, double> & DistanceOnTheFly::GetNeighborhood(int vertexID)
{
	assert(m_fDistanceThreshold>=0);
	if (m_PrecomputedSmallNeighborhoood[vertexID].size()==0)
		PrecomputeRow(vertexID);
	return m_PrecomputedSmallNeighborhoood[vertexID];
}

double DistanceOnTheFly::GetPrecomp(int i, int j)
{
	double retVal = -1.;
	// check if existing distance metric has this value cached
	if (m_pTryMeFirst!=NULL && m_pTryMeFirst->Valid(SurfaceSample(i, m_pSurface->GetMesh()), 
													SurfaceSample(j, m_pSurface->GetMesh())))
	{
		retVal = m_pTryMeFirst->Distance(SurfaceSample(i, m_pSurface->GetMesh()), 
										 SurfaceSample(j, m_pSurface->GetMesh()));
		assert(retVal>=0);
		return retVal;
	}
	
	if (m_fDistanceThreshold==-1)
	{
		// check if I have this value cached
		if (m_PrecomputedRow[i]!=NULL)
		{
			retVal = m_PrecomputedRow[i][j];
			assert(retVal>=0);
			return retVal;
		}
		
		if (m_PrecomputedRow[j]!=NULL)
		{
			retVal = m_PrecomputedRow[j][i];
			assert(retVal>=0);
			return retVal;
		}
		
		// value is not cached - create new row (potentially erase something
		PrecomputeRow(i);
		assert(m_PrecomputedRow[i]!=NULL);
		retVal = m_PrecomputedRow[i][j];;
		assert(retVal>=0);
		return retVal;
	}
	else
	{
		if ((int)m_PrecomputedSmallNeighborhoood[i].size()>0)
		{
			retVal = DistanceGeodesic::GetThresholdedDistances(j, m_pSurface->GetMesh(), 
															   m_PrecomputedSmallNeighborhoood[i], 
															   m_fDistanceThreshold);
			assert(retVal>=0);
			return retVal;
		}
		
		if (m_PrecomputedSmallNeighborhoood[j].size()>0)
		{
			retVal = DistanceGeodesic::GetThresholdedDistances(i, m_pSurface->GetMesh(), 
															   m_PrecomputedSmallNeighborhoood[j], 
															   m_fDistanceThreshold);
			assert(retVal>=0);
			return retVal;
		}
		PrecomputeRow(i);
		retVal = DistanceGeodesic::GetThresholdedDistances(i, m_pSurface->GetMesh(), 
														   m_PrecomputedSmallNeighborhoood[j], 
														   m_fDistanceThreshold);
		assert(retVal>=0);
		return retVal;
	}
	assert(false);
	return -1;
}

void DistanceOnTheFly::SetPrecomp(int i, int j, double val)
{
	assert(false);
}

int DistanceOnTheFly::GetCacheSize()
{
	return (int)m_RowAddedInOrder.size();
}




