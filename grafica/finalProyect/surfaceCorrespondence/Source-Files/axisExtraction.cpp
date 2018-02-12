#include "Angel-yjc.h"
#include <string>
#include <vector>
#include <fstream>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>
#include <cmath>
#include "ppm.h"
#include "texture.h"
#include "stb_image.h"
#include "model.h"
#include "axisExtraction.h"
using namespace std;

/*
MinimumEuclideanDistances(R3Mesh* mesh, vector<R3Point>* points)
{
	RNLength* dists = new RNLength[mesh->NVertices()];
	for (int i=0; i<mesh->NVertices(); i++)
	{
		R3MeshVertex* vertex = mesh->Vertex(i);
		R3Point pos = mesh->VertexPosition(vertex);
		RNLength min_dist = FLT_MAX;
		for (int j=0; j<int(points->size()); j++)
		{
			R3Point _p = (*points)[j];
			RNLength d = R3Distance(pos, _p);
			if (d < min_dist)
				min_dist = d;
		}
		dists[i] = min_dist;
	}
	return dists;
}
*/


MinimumEuclideanDistances(model2* g_model, vector<point3>* point){



	for (int i = 0; i < mesh->; ++i)
	{
		/* code */
	}



}

