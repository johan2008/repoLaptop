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






double R3Distance(const point3 span1, const point3 span2)
{

	double dx = span1.x - span2.x;
	double dy = span1.y - span2.y;
	double dz = span1.z - span2.z;


	return sqrt(dx*dx + dy*dy + dz*dz);

}



double MinimumEuclideanDistances(model* g_model, vector<point3>* point){

	double min_dist;

	for (int i = 0; i < g_model->nv; ++i)
	{
		point3 vertex = g_model->positions[i]; 
		min_dist = 100000;
		for (int j = 0; j < g_model->nv; ++i)
		{
			double aux_dist = R3Distance(g_model->positions[i], g_model->positions[j]) ;
			

		}
	}

}


void axisExtraction::printData(model* g_model){

}


double max3(double a , double b, double c){
	double max = a;
	if(max<b) max = b;
	if(max<c) max = c;

	return max;
}


point3 axisExtraction::create_plane(model* g_model){
	double x = 0, y = 0, z=0;
	point3 sum = 0;
	for (int i = 0; i < g_model->nv; ++i)
	{
		/* code */
		sum = sum + g_model->positions[i];
	}
	point3 centroide = sum*(1.0/g_model->nv);

	double xx = 0.0;
	double xy = 0.0;
	double xz = 0.0;
	double yy = 0.0;
	double yz = 0.0;	
	double zz = 0.0;

	point3 r = 0;
	for (int i = 0; i < g_model->nv; ++i)
	{
		r = g_model->positions[i] - centroide;
		xx += r.x * r.x;
		xy += r.x * r.y;
		xz += r.x * r.z;
		yy += r.y * r.y;
		yz += r.y * r.z;
		zz += r.z * r.z;
	}
	double det_x = yy*zz - yz*yz;
	double det_y = xx*zz - xz*xz;
	double det_z = xx*yy - xy*xy;
	double det_max = max3(det_x, det_y, det_z);

	if(det_max == det_x){
		x = det_x;
		y = xz*yz -xy*zz;
		z = xy*yz - xz*yy;
	}else if(det_max == det_y){
		x = xz*yz - xy*zz;
		y = det_y;
		z = xy*xz - yz*xx;
	}else{
		x = xy*yz - xz*yy;
		y = xy*xz - yz*xx;
		z = det_z;
	}
	cout<<"x: "<<x<<" y: "<<y<<" z: "<<z<<endl;
	return centroide;
}



point3 axisExtraction::getCentroide(model* g_model){
	double x = 0, y = 0, z=0;
	point3 sum = 0;
	for (int i = 0; i < g_model->nv; ++i)
	{
		sum = sum + g_model->positions[i];
	}
	point3 centroide = sum*(1.0/g_model->nv);

	return centroide;
}

	
double areaMesh(model* g_model){
	double A = 0;
	for (int i = 0; i < g_model->nf; ++i)
	{
		point3 e1, e2, e3;
		e1.x = g_model->positions[g_model->facades[i].v2].x - g_model->positions[g_model->facades[i].v1].x;
		e1.y = g_model->positions[g_model->facades[i].v2].y - g_model->positions[g_model->facades[i].v1].y;
		e1.z = g_model->positions[g_model->facades[i].v2].z - g_model->positions[g_model->facades[i].v1].z;
		
		e2.x = g_model->positions[g_model->facades[i].v3].x - g_model->positions[g_model->facades[i].v1].x;
		e2.y = g_model->positions[g_model->facades[i].v3].y - g_model->positions[g_model->facades[i].v1].y;
		e2.z = g_model->positions[g_model->facades[i].v3].z - g_model->positions[g_model->facades[i].v1].z;

		e3.x = e1.y*e2.z - e1.z*e2.y;
		e3.y = e1.z*e2.x - e1.x*e2.z;
		e3.z = e1.x*e2.y - e1.y*e2.x;

		A += 0.5*sqrt(e3.x*e3.x + e3.y*e3.y + e3.z*e3.z);

	}

	return A;
}


double axisExtraction::averageDistance(model* g_model){
	double area = areaMesh(g_model);
	cout<<"AREA TOTAL:  "<<areaMesh(g_model)<<endl;
	return 0.03*sqrt(area);
}

double axisExtraction::epsilon(model* g_model){
	double area = areaMesh(g_model);
	cout<<"AREA TOTAL:  "<<areaMesh(g_model)<<endl;
	return 0.03*sqrt(area);
}


void axisExtraction::cero_level(model* g_model){
	int i = 0;
	for(int n=0; n< g_model->nv; n++)
    {
        if( abs(g_model->positions[n].z) < 0.01 ){
        	//cout<<"0- level"<<endl;
        	i++;
        }
    }
    
} 
