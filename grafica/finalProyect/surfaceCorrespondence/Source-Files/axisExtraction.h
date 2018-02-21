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
using namespace std;


class axisExtraction{
	public:
		//axisExtraction();
		double MinimumEuclideanDistances(model* g_model, vector<point3>* point);
		void printData(	model* g_model);
		point3 create_plane(model* g_model);
		point3 getCentroide(model* g_model);
		double averageDistance(model* g_model);

		

};
