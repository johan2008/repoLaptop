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

class axisAlignment{
	public:
		//axisExtraction();
		void printData(	model* g_model);
		vector<curveStruct> find_curve(vector<curveStruct> curves1, vector<curveStruct> curves2);


};
