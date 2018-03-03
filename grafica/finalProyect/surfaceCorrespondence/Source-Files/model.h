#pragma once
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
using namespace std;


typedef Angel::vec3  point3;
typedef Angel::vec3  color3;



struct facade
{
	GLuint v1;
	GLuint v2;
	GLuint v3;
};

struct curveStruct{
    vector<point3> curve;
    vector<color3> color;
    double q;
};


class model{
	public:
			point3 *positions;
			facade *facades;
			color3 *colors;
			int nv, nf;
			model(string a,  int b , int scale);

};
