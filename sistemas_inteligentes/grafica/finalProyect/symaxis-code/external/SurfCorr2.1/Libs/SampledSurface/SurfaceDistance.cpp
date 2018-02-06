#include "SurfaceDistance.h"
#include "SampledSurface.h"

SurfaceDistance::SurfaceDistance(SampledSurface * surface, bool vertexToVertexOnly)
{
	m_pSurface = surface;
	m_bVertexToVertex = vertexToVertexOnly;
	m_aPrecomputedDistances = NULL;
	m_iPrecompWidth=-1;
	m_iPrecompHeight=-1;
	
}

SurfaceDistance::~SurfaceDistance()
{
	delete [] m_aPrecomputedDistances;
}

double SurfaceDistance::CalculateDistance(const SurfaceSample & s1, const SurfaceSample & s2) const
{
	assert(false);
	return -1;
}

double SurfaceDistance::GetPrecomp(int i, int j)
{
	assert(m_aPrecomputedDistances!=NULL);
	assert(j>=0 && i>=0 && i<m_iPrecompWidth && j<m_iPrecompHeight);
	assert(m_iPrecompWidth>=0 && m_iPrecompHeight>=0);
	return m_aPrecomputedDistances[j * m_iPrecompWidth + i];
}

void SurfaceDistance::SetPrecomp(int i, int j, double val)
{
	m_aPrecomputedDistances[j * m_iPrecompWidth + i] = val;
}

bool SurfaceDistance::Valid(const SurfaceSample &, const SurfaceSample &) const
{
	return true;
}

double SurfaceDistance::Distance(const SurfaceSample &s1, const SurfaceSample &s2)
{	
	assert(!s1.Invalid() && !s2.Invalid());

//	assert(s1.CheckMeshIsSame(m_pSurface->GetMesh()));
//	assert(s2.CheckMeshIsSame(m_pSurface->GetMesh()));
	// NOTE: this checked was removed because symmetric surfaces share distance
	
	if (!m_bVertexToVertex)
	{
		double retVal = CalculateDistance(s1, s2);;
		assert(retVal>=0);
		return retVal;
	}
	else
	{
		if (m_aPrecomputedDistances==NULL)
			PrecomputeDistances();
		
		bool atVertex2;
		int v2 = s2.NearestVertex(&atVertex2);
		if (atVertex2)	// one of samples is on vertex
		{
			bool atVertex1; 
			int v1 = s1.NearestVertex(&atVertex1);
			if (atVertex1)	// another sample is on vertex
			{
				double retVal = GetPrecomp(v1, v2);
				assert(retVal>=0);
//				std::cout<<"\tDistance("<<s1.NearestVertex();
//				std::cout<<", "<<s2.NearestVertex()<<") = "<<retVal<<std::endl;				
				return retVal;
			}
			else	// NOTE: only s2 is on vertex
			{
				double retVal = DistanceToVertex(s1, v2);
				assert(retVal>=0);
				return retVal;
			}
		}	
		else
		{
			bool atVertex1; 
			int v1 = s1.NearestVertex(&atVertex1);
			if (atVertex1)
			{
				double retVal = DistanceToVertex(s2, v1);
				assert(retVal>=0);
				return retVal;
			}
			else
			{
				double d1 = DistanceToVertex(s1, s2.VertID(0));
				double d2 = DistanceToVertex(s1, s2.VertID(1));
				double d3 = DistanceToVertex(s1, s2.VertID(2));
				double retVal = s2.Interpolate(d1, d2, d3);		
				assert(retVal>=0);
				return retVal;
			}
		}
	}
	assert(false);
	return -1.;
}

double SurfaceDistance::DistanceToVertex(const SurfaceSample & s1, int v)
{
	double d1 = GetPrecomp(v, s1.VertID(0));
	assert(d1>=0);
	double d2 = GetPrecomp(v, s1.VertID(1));
	assert(d2>=0);	
	double d3 = GetPrecomp(v, s1.VertID(2));	
	assert(d3>=0);	
	if (s1.B(0)<0 || s1.B(1)<0 || s1.B(2)<0)
		std::cout<<"s.b=["<<s1.B(0)<<", "<<s1.B(1)<<", "<<s1.B(2)<<"]"<<std::endl;
	return s1.Interpolate(d1, d2, d3);
}

void SurfaceDistance::PrecomputeDistances()
{
	R3Mesh * mesh = m_pSurface->GetMesh();
	m_aPrecomputedDistances = new double[mesh->NVertices() * mesh->NVertices()];
	m_iPrecompWidth = mesh->NVertices(); 
	m_iPrecompHeight = mesh->NVertices();
	for (int i=0; i<mesh->NVertices(); i++)
	for (int j=0; j<mesh->NVertices(); j++)
		SetPrecomp(i, j, CalculateDistance(SurfaceSample(i,mesh), SurfaceSample(j,mesh)));
}

void SurfaceDistance::WritePrecomputed(VKString outFile)
{
	assert(m_aPrecomputedDistances!=NULL);
	
	FILE * dataStream = fopen(outFile.c_str(), "wb");
	assert(dataStream!=NULL);
	int nVert =  m_pSurface->GetMesh()->NVertices();
	int format=0;
	// num vertices - compare to mesh's num vertices	
	fwrite( &nVert, sizeof(int), 1, dataStream);	
	// dist type - 0 by default		
	fwrite(&format, sizeof(int), 1, dataStream);
	
	// grid dimensions
	fwrite( &m_iPrecompWidth, sizeof(int), 1, dataStream);
	fwrite( &m_iPrecompHeight, sizeof(int), 1, dataStream);	

	// grd values		
	fwrite(m_aPrecomputedDistances, sizeof(double), m_iPrecompWidth*m_iPrecompHeight, dataStream);
	
	fclose(dataStream);
}

bool SurfaceDistance::LoadPrecomputed(VKString inFile)
{
	
	FILE * dataStream = fopen(inFile.c_str(), "rb");
	if(dataStream==NULL)
		return false;
	
	if (m_aPrecomputedDistances!=NULL)
		delete [] m_aPrecomputedDistances;
	
	int nVert =  0;
	int format=0;
	//	// num vertices - compare to mesh's num vertices	
	assert(fread( &nVert, sizeof(int), 1, dataStream )!=0);
	//	// dist type - 0 by default		
	assert(fread(&format, sizeof(int), 1, dataStream)!=0);
	
	assert(nVert==m_pSurface->GetMesh()->NVertices());
	
	// grid dimensions
	assert(fread(&m_iPrecompWidth, sizeof(int), 1, dataStream)!=0);
	assert(fread(&m_iPrecompHeight, sizeof(int), 1, dataStream)!=0);
	
	//	// grd values	
	m_aPrecomputedDistances = new double[m_iPrecompHeight*m_iPrecompWidth]; 
	assert(fread(m_aPrecomputedDistances, sizeof(double), m_iPrecompHeight*m_iPrecompWidth, dataStream)!=0);

	fclose(dataStream);
	
	return true;
}




