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



model::model(string fileName, int numImg, int escale){
	string readLine;
	int delPos1, delPos2, delPos3, delPos4;
	ifstream in(fileName.c_str());
	getline(in,readLine);
	if(readLine!="OFF"){
		cout<<"not off format"<<endl;
	}

	getline(in,readLine);
	delPos1 = readLine.find(" ",0);
	nv = atoi(readLine.substr(0,delPos1+1).c_str());

	delPos2 = readLine.find(" ",delPos1);
	nf = atoi(readLine.substr(delPos1,delPos2+1).c_str());

	positions = new point3[nv];
	
	for(int n=0; n<nv; n++){
		getline(in,readLine);
		delPos1 = readLine.find(" ",0);
		positions[n].x = atof(readLine.substr(0,delPos1).c_str())/escale +3*numImg;
		//cout<<"positions x:   "<<positions[n].x<<endl;
		delPos2 = readLine.find(" ", delPos1+1);
		positions[n].y = atof(readLine.substr(delPos1,delPos2 ).c_str())/escale;
		//cout<<"positions y:   "<<positions[n].y<<endl;
		delPos3 = readLine.find(" ", delPos2+1);
		positions[n].z = atof(readLine.substr(delPos2,delPos3 ).c_str())/escale;
		//cout<<"positions z:   "<<positions[n].z<<endl;
	}		 
	facades = new facade[nf];
	for (int n = 0; n < nf; ++n)
	{
		getline(in,readLine);
		delPos1 = readLine.find(" ",0);
		delPos2 = readLine.find(" ",delPos1+1);
		facades[n].v1 = atoi(readLine.substr(delPos1,delPos2 ).c_str());
		delPos3 = readLine.find(" ",delPos2+1);
		facades[n].v2 = atoi(readLine.substr(delPos2,delPos3 ).c_str());
		delPos4 = readLine.find(" ",delPos3+1);
		facades[n].v3 = atoi(readLine.substr(delPos3,delPos4 ).c_str());
	}
	cout<<"nv: "<<nv<<endl;;
}

