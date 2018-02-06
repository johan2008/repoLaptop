#include "FeatureCurvature.h"
#include "FeatureMGD.h"
#include "AnalysisStats.h"

#include "SampledSurface.h"
#include "SurfaceDistance.h"
#include "DistanceGeodesic.h"
#include "DistanceOnTheFly.h"
#include "VkFunctions.h"

///////////////////// Feature Calculation - Fast Calculation ///////////////
FeatureCurvature::FeatureCurvature(SampledSurface * surface, int numInitialSamples, 
								   const VKString & featurename)
: SurfaceFeature(surface, true, featurename)
{
	m_iNumInitialSamples = numInitialSamples;
}

FeatureCurvature::~FeatureCurvature()
{
}

VKString FeatureCurvature::GetFeatureName()
{
	return "FeatureCurvature";
}

void FeatureCurvature::PrecomputeVertexFeatures()
{
	assert(m_pPerVertexValues!=NULL);
	
	// spread N points
	R3Mesh * mesh = m_pSurface->GetMesh();	
	std::vector<int> evenVertices;
//	FeatureMGD::EvenlySpreadVertices(m_pSurface, m_iNumInitialSamples, evenVertices);
	for (int i=0; i<mesh->NVertices(); i++)
		evenVertices.push_back(i);
	std::cout<<"Calculating Curvature (NOT subsampled)"<<std::endl;	
	
	int MAX_ITER = evenVertices.size();
	// calculate Curvature value at evenly spread vertices
	for (int i=0; i<(int)evenVertices.size(); i++)
	{
		if (i%10==0)
			AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("CreatingSurface_Curvature", i, MAX_ITER);
		R3MeshVertex * vertex = mesh->Vertex(evenVertices[i]);
		double curvature = pow(vkAbs(mesh->VertexCurvature(vertex)), .2);
		//double curvature = vkAbs(mesh->VertexCurvature(vertex));
		m_pPerVertexValues->SetVertexValue(evenVertices[i], curvature);
	}

	//m_pPerVertexValues->SmoothGaussianFromKnown(.005);
	//m_pPerVertexValues->SmoothGaussianFromKnown(.02);
	//m_pPerVertexValues->SmoothGaussianFromKnown(-1);
	//m_pPerVertexValues->SmoothGaussianFromKnownBySplatting(, .01, .01);
	m_pPerVertexValues->SmoothSimpleLaplacian(10);
	std::cout<<std::endl;	
}

double FeatureCurvature::CalculateValue(const SurfaceSample & sample) const
{
	assert(false);
}




