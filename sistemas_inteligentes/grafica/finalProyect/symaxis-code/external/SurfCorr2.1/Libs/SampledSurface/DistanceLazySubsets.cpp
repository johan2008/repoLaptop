#include "DistanceLazySubsets.h"
#include "SampledSurface.h"
#include "VkFunctions.h"
#include "AnalysisStats.h"

DistanceLazySubsets::DistanceLazySubsets(SampledSurface * surface, 
										 SurfaceSampleSet * fromSet, SurfaceSampleSet * toSet,
										 DistanceGeodesic::CalculationType method)
: DistanceGeodesic(surface, method)
{
	assert(fromSet!=NULL);
	m_iPrecompHeight = fromSet->NumSamples();
	if (toSet==NULL)
		m_iPrecompWidth = surface->GetMesh()->NVertices();
	else
		m_iPrecompWidth = toSet->NumSamples();
	m_aPrecomputedDistances = NULL;

	m_MapVertexFrom = new int[surface->GetMesh()->NVertices()];
	m_MapVertexTo = new int[surface->GetMesh()->NVertices()];
	
	for (int i=0; i<surface->GetMesh()->NVertices(); i++)
	{
		m_MapVertexFrom[i] = -1;
		if (toSet!=NULL)
			m_MapVertexTo[i] = -1;
		else
			m_MapVertexTo[i] = i;
	}
	
	for (int i=0; i<fromSet->NumSamples(); i++)
	{
		bool ok;
		int id = fromSet->GetSample(i).NearestVertex(&ok);
		assert(ok);
		m_MapVertexFrom[id] = i;
	}
	
	if (toSet!=NULL)
	{
		for (int i=0; i<toSet->NumSamples(); i++)
		{
			bool ok;
			int id = toSet->GetSample(i).NearestVertex(&ok);
			assert(ok);
			m_MapVertexTo[id] = i;
		}
	}
}

DistanceLazySubsets::~DistanceLazySubsets()
{
	delete [] m_MapVertexFrom;
	delete [] m_MapVertexTo;
}

void DistanceLazySubsets::PrecomputeFunkhouser()
{
	// allocate vertex data
	if (m_bVerbous)	
		std::cout<<"\tLazy Subset ("<<m_iPrecompWidth<<" x "<<m_iPrecompHeight<<"): "<<flush;

	GeoVertexData * vertex_data = DistanceGeodesic::InitializeFunkPrecomputation(m_pSurface->GetMesh(), m_bVerbous);
		
	
	std::vector<int> verticesToRemember;
	for (int i=0; i<m_iPrecompWidth; i++)
		verticesToRemember.push_back(-1);
	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
		if (m_MapVertexTo[i]!=-1)
			verticesToRemember[m_MapVertexTo[i]]=i;
	
	std::vector<int> verticesFrom;
	for (int i=0; i<m_iPrecompHeight; i++)
		verticesFrom.push_back(-1);
	for (int i=0; i<m_pSurface->GetMesh()->NVertices(); i++)
		if (m_MapVertexFrom[i]!=-1)
			verticesFrom[m_MapVertexFrom[i]]=i;
	
	for (int i=0; i<m_iPrecompWidth * m_iPrecompHeight; i++)
		m_aPrecomputedDistances[i] = -1;
	
	//std::cout<<"Vertices to remember = "<<verticesToRemember.size()<<std::endl;
	
	// Run Dijkstra for each from vertex
	for (int i=0; i<m_iPrecompHeight; i++)
	{
//		if (m_bVerbous)
//		{	
//			if (i%10000==0 && i!=0)
//				std::cout<<"\t"<<i<<" / "<<m_iPrecompHeight<<"\n\t\t"<<std::flush;
//			else if (i%500==0 && i!=0)
//				std::cout<<"."<<std::flush;
//		}
		if (i%10==0)
			AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("CreatingSurface_Distances",
																i, m_iPrecompHeight);

		assert(verticesFrom[i]!=-1);
		DistanceGeodesic::PrecomputeFillFunkhouser(vertex_data, m_pSurface->GetMesh(), verticesFrom[i], 
												   &(m_aPrecomputedDistances[i * m_iPrecompWidth]), 
												   verticesToRemember);
	}

	delete [] vertex_data;
//	if (m_bVerbous)
//		std::cout<<" - Done!"<<std::endl;	
	std::cout<<std::endl;
}

void DistanceLazySubsets::PrecomputeKirsanov()
{
	assert(false);	//todo see as above
	R3Mesh * mesh = m_pSurface->GetMesh();
	geodesic::Mesh * kirsanovMesh;
	geodesic::GeodesicAlgorithmBase * algorithm;
	InitializeKirsanovPrecomputation(&kirsanovMesh, &algorithm);
	
	std::vector<int> verticesToRemember;
	for (int i=0; i<m_iPrecompHeight; i++)
		verticesToRemember.push_back(-1);
	for (int i=0; i<mesh->NVertices(); i++)
		if (m_MapVertexTo[i]!=-1)
			verticesToRemember[m_MapVertexTo[i]]=i;	
	
	assert(m_aPrecomputedDistances!=NULL);
	// run Kirsanov code for each vertex
	for (int i=0; i<m_iPrecompWidth; i++)
	{
		if (m_bVerbous)
		{
			if (i%10000==0 && i!=0)
				std::cout<<"\t"<<i<<" / "<<mesh->NVertices()<<"\n\t\t"<<std::flush;
			else if (i%50==0 && i!=0)
				std::cout<<"."<<std::flush;
		}
		
		PrecomputeFillRowKirsanov(m_MapVertexFrom[i], algorithm, kirsanovMesh, 
								  &(m_aPrecomputedDistances[i * m_iPrecompWidth]), verticesToRemember);
	}	
	FreeKirsanovPrecomputation(kirsanovMesh, algorithm);
}

double DistanceLazySubsets::GetPrecomp(int i, int j)
{
	if (m_MapVertexFrom[i]==-1 || m_MapVertexTo[j]==-1)
	{
		if (m_MapVertexFrom[j]==-1 || m_MapVertexTo[i]==-1)
		{
			std::cout<<"[ERROR] DistanceLazySubsets Querying: ["<<i<<", "<<j<<"] maps to ";
			std::cout<<"["<<m_MapVertexFrom[j]<<", "<<m_MapVertexTo[i]<<"]"<<std::flush;
			std::cout<<" or ["<<m_MapVertexFrom[i]<<", "<<m_MapVertexTo[j]<<"]"<<std::endl;
			assert(false);
		}
		double retVal = m_aPrecomputedDistances[m_MapVertexFrom[j]*m_iPrecompWidth+m_MapVertexTo[i]];
		assert(retVal>=0);
		return retVal;
	}
	else
	{
		assert(m_MapVertexTo[j]!=-1);
		double retVal = m_aPrecomputedDistances[m_MapVertexFrom[i]*m_iPrecompWidth+m_MapVertexTo[j]];
		assert(retVal>=0);
		return retVal;
	}
}

void DistanceLazySubsets::SetPrecomp(int i, int j, double val)
{
	assert(false);	// no need
}

bool DistanceLazySubsets::Valid(const SurfaceSample & s1, const SurfaceSample & s2) const
{
	int v11, v12=-1, v13=-1;
	int v21, v22=-1, v23=-1;
	bool nearest;
	v11 = s1.NearestVertex(&nearest);
	if (!nearest)
	{
		v11 = s1.VertID(0);		v12 = s1.VertID(1);		v13 = s1.VertID(2);		
	}
	v21 = s2.NearestVertex(&nearest);
	if (!nearest)
	{
		v21 = s2.VertID(0);		v22 = s2.VertID(1);		v23 = s2.VertID(2);		
	}

	
	return DistanceExists(v11, v21) && DistanceExists(v11, v22) && DistanceExists(v11, v23) 
			&& DistanceExists(v12, v21) && DistanceExists(v12, v22) && DistanceExists(v12, v23)
			&& DistanceExists(v13, v21) && DistanceExists(v13, v22) && DistanceExists(v13, v23);	
}

bool DistanceLazySubsets::DistanceExists(int i, int j) const
{
	if (i==-1 || j==-1)
		return true;	// special exception see functio above
	return ((m_MapVertexFrom[i]!=-1 && m_MapVertexTo[j]!=-1)
			|| (m_MapVertexFrom[j]!=-1 && m_MapVertexTo[i]!=-1));
}


