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
#include "axisAlignment.h"
#define N 200
using namespace std;




double min3(double a , double b, double c){
	double min = a;
	if(min>b) min = b;
	if(min>c) min = c;
	return min;
}

double min2(double a , double b){
	double min = a;
	if(min>b) min = b;
	return min;
}


void axisAlignment::printData(model* g_model)
{


}

double distance(point3 a, point3 b){
    return sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y)  + (a.z-b.z)*(a.z-b.z)  );
}


double landa(point3 m, point3 n){
	if(m == 0 || n == 0) return 0.5;
	else return distance(m,n);
}


double editDistance(vector<point3> a, vector<point3> b)
{
	int d[N][N];
	int e = 0;

	d[0][0] = landa(a[0],b[0]);

	for (int i = 1; i < N; ++i)
	{
		for (int j = 1; j < N; ++j)
		{
			d[i][j] = min3 (   
							( d[i-1][j-1]+ landa(a[i],b[j]) ) ,       
							( d[i-1][j]+ min2(landa(a[i],b[j]),landa(a[i],e))    ) ,
							( d[i][j-1] + min2(landa(a[i],b[j] ),(landa(e,b[j]) )) )
							);
		}	
	}
	return d[N-1][N-1];
}


vector<curveStruct> axisAlignment::find_curve(vector<curveStruct> curves1, vector<curveStruct> curves2 )
{

	double maxQuality = -1;
	double aux;
	vector<curveStruct> c;
	curveStruct c1;
	curveStruct c2;
	for (int a = 0; a < curves1.size(); ++a)
	{
		for (int b = 0; b < curves2.size(); ++b)
		{
			aux = curves1[a].q*curves2[b].q;//*(1.0/(editDistance(curves1[a].curve, curves2[b].curve) ));
			if(maxQuality<aux){
				c1 = curves1[a];
				c2 = curves2[b];
			}

		}
	}
	c.push_back(c1);
	c.push_back(c2);
	return c;

}
