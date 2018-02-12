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

#define CLEAR_COLOR 0.529f, 0.807f, 0.92f, 1.0f

typedef Angel::vec3  color3;
typedef Angel::vec3  point3;
typedef Angel::vec2  point2;

int windowWidth = 600;
int windowHeight = 600;
GLuint program;
GLuint floor_buffer;
GLuint floor_buffer2;
GLuint EBO;
GLuint EBO2;


GLfloat angle = 0.0;

GLuint  u_tView;
GLfloat  fovy = 45.0;

vec3 init_eye(7.0, 3.0, -10.0);
vec3 eye = init_eye;
vec3 _cameraDir( -7.0f, -3.0f, 10.0f );

GLfloat  aspect; 

const int floor_NumVertices = 6; //(1 face)*(2 triangles/face)*(3 vertices/triangle)
point3 floor_points[floor_NumVertices]; // positions for all vertices
color3 floor_colors[floor_NumVertices]; // colors for all vertices
point3 floor_normal[floor_NumVertices];

float pasoX = 0;
float pasoZ = 0;
float pasoY = 0;



/*
struct  vectex
{
	float x;
	float y;
	float z;
};


struct facade
{
	GLuint v1;
	GLuint v2;
	GLuint v3;
};


struct model{
	point3 * positions;
	facade *facades;
	color3 * colors;

	int nv, nf;

	model(string fileName, int numImg){
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
		cout<<"nv:   "<<nv<<endl;

		delPos2 = readLine.find(" ",delPos1);
		nf = atoi(readLine.substr(delPos1,delPos2+1).c_str());
		cout<<"nf:   "<<nf<<endl;

		positions = new point3[nv];
		
		for(int n=0; n<nv; n++){
			getline(in,readLine);
			delPos1 = readLine.find(" ",0);
			positions[n].x = atof(readLine.substr(0,delPos1).c_str()) +3*numImg;
			//cout<<"positions x:   "<<positions[n].x<<endl;
			delPos2 = readLine.find(" ", delPos1+1);
			positions[n].y = atof(readLine.substr(delPos1,delPos2 ).c_str());
			//cout<<"positions y:   "<<positions[n].y<<endl;
			delPos3 = readLine.find(" ", delPos2+1);
			positions[n].z = atof(readLine.substr(delPos2,delPos3 ).c_str());
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
	}

};

*/


model* g_model;

model* g_model2;

void init(){
    //model g_model("386.off");
    g_model  = new model( "001.off" , 0);
    g_model2 = new model( "000.off" , 1); 

    // Create and initialize a vertex buffer object for floor, to be used in display()
    glGenBuffers(1, &floor_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, floor_buffer);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(floor_points),NULL, GL_STATIC_DRAW);

    glBufferData(GL_ARRAY_BUFFER, sizeof(point3) * g_model->nv,g_model->positions, GL_STATIC_DRAW);

    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(facade) * g_model->nf, g_model->facades, GL_STATIC_DRAW); 
    //glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeof(g_model->facades), g_model->facades);



    glGenBuffers(1, &floor_buffer2);
    glBindBuffer(GL_ARRAY_BUFFER, floor_buffer2);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(floor_points),NULL, GL_STATIC_DRAW);

    glBufferData(GL_ARRAY_BUFFER, sizeof(point3) * g_model2->nv,g_model2->positions, GL_STATIC_DRAW);

    glGenBuffers(1, &EBO2);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO2);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(facade) * g_model2->nf, g_model2->facades, GL_STATIC_DRAW); 





    program = InitShader("../Shaders/vshader42.glsl", "../Shaders/fshader42.glsl");    
    glEnable( GL_DEPTH_TEST );
    glClearColor( CLEAR_COLOR); 
    glLineWidth(2.0);
}


void drawObj(GLuint vertexBuffer, GLuint indexBuffer, int numFaces){
	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
	GLuint vPosition = glGetAttribLocation(program, "vPosition");
    glEnableVertexAttribArray(vPosition);
    glVertexAttribPointer(vPosition, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(0) );
    //glDrawArrays(GL_TRIANGLES, 0, num_vertices);
    glDrawElements(GL_TRIANGLES, numFaces*3, GL_UNSIGNED_INT, 0);

    glDisableVertexAttribArray(vPosition);
}


void display(){
    GLuint  u_tModel;  
    GLuint  u_tProj;  
    u_tView;
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glUseProgram(program);
    u_tModel = glGetUniformLocation( program, "u_tModel" );
    u_tProj =  glGetUniformLocation( program, "u_tProj" );
    u_tView =  glGetUniformLocation( program, "u_tView" );
	mat4  matrixProj = Perspective(fovy, 1, 0.1f, 100.0f );
    glUniformMatrix4fv(u_tProj, 1, GL_TRUE, matrixProj); // GL_TRUE: matrix is row-major

    glUniform4fv(glGetUniformLocation(program,"CameraEye" ),1 , eye);

    vec3 _cameraTarget = eye + _cameraDir;
    vec3 _worldUp( 0.0f, 1.0f, 0.0f );

    mat4  matrixView= LookAt( eye,_cameraTarget,_worldUp );
    glUniformMatrix4fv( u_tView, 1, GL_TRUE, matrixView );
    mat4 matrixModel= mat4(1.0f);
    matrixModel = Translate(0.0, 0.0, 0.0) * Scale (1.0, 1.0, 1.0) * Rotate(angle, 0.0, 2.0, 0.0);


    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModel); // GL_TRUE: matrix is row-major
    
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawObj(floor_buffer, EBO, g_model->nf);


    drawObj(floor_buffer2, EBO2, g_model2->nf);

    glutSwapBuffers();
}


void reshape(int width, int height){
    glViewport(0, 0, width, height);
    aspect = (GLfloat) width  / (GLfloat) height;
    glutPostRedisplay();
}


void idle(){

}


void keyboard(unsigned char key, int x, int y){
    switch(key) {
	    case 'X': eye[0] += 1.0; _cameraDir[0] -= 1.0; break;
        case 'x': eye[0] -= 1.0; _cameraDir[0] += 1.0; break;
        case 'Y': eye[1] += 1.0; _cameraDir[1] -= 1.0; break;
        case 'y': eye[1] -= 1.0; _cameraDir[1] += 1.0; break;
        case 'Z': eye[2] += 1.0; _cameraDir[2] -= 1.0; break;
        case 'z': eye[2] -= 1.0; _cameraDir[2] += 1.0; break;
    }
    glutPostRedisplay();
}


int main(int argc, char **argv){

	glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("Main Window");
    int err = glewInit();
    if (GLEW_OK != err){ 
        printf("Error: glewInit failed: %s\n", (char*) glewGetErrorString(err)); 
        exit(1);
    }

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);
    //addControl();
    init();
    glutMainLoop();
    return 0;

}