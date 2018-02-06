#ifndef __MAP_EUCLIDEAN_H
#define __MAP_EUCLIDEAN_H

#include <fstream>
#include <sstream>

#include "gaps.h"

#include "VKString.h"

#include "SampledSurface.h"
#include "SurfaceMap.h"

using namespace std;


class MapEuclidean : public SurfaceMap
{
	public:
		MapEuclidean(SampledSurface * M1, SampledSurface * M2);

		SurfaceSample ForwardMap(const SurfaceSample & s);
		SurfaceSample InverseMap(const SurfaceSample & s);

		virtual VKString GetSurfaceMapType();

		virtual void SaveMap(std::ofstream & outFilename);
		virtual void LoadMap(std::ifstream & outFilename);

	protected:
		static SurfaceSample MapEuclideanDist(const SurfaceSample & s, 
					R3MeshSearchTree * nearestOnMesh, R3Mesh * corrMesh);

		R3MeshSearchTree * m_pSearchTreeOnM1;
		R3MeshSearchTree * m_pSearchTreeOnM2;
};

#endif




