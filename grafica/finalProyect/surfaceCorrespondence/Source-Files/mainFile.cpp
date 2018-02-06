//g++ rotate-cube-new.cpp InitShader.cpp -o rotate-cube-new -lglut -lGLEW -lGL -lGLU && ./rotate-cube-new
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

#define CLEAR_COLOR 0.529f, 0.807f, 0.92f, 1.0f

typedef Angel::vec3  color3;
typedef Angel::vec3  point3;
typedef Angel::vec2  point2;

int windowWidth = 600;
int windowHeight = 600;
GLuint program;
GLuint floor_buffer;
GLuint EBO;


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


void floor()
{
    floor_colors[0] = color3(0,1,0); floor_points[0] = point3(5,0,8);   floor_normal[0] = point3(0,1,0);
    floor_colors[1] = color3(0,1,0); floor_points[1] = point3(5,0,-4);  floor_normal[1] = point3(0,1,0);
    floor_colors[2] = color3(0,1,0); floor_points[2] = point3(-5,0,-4); floor_normal[2] = point3(0,1,0);

    floor_colors[3] = color3(0,1,0); floor_points[3] = point3(-5,0,8);  floor_normal[3] = point3(0,1,0);
    floor_colors[4] = color3(0,1,0); floor_points[4] = point3(-5,0,-4); floor_normal[4] = point3(0,1,0);
    floor_colors[5] = color3(0,1,0); floor_points[5] = point3(5,0,8);   floor_normal[5] = point3(0,1,0);
}


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

	model(string fileName){
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
			positions[n].x = atof(readLine.substr(0,delPos1).c_str());
			//cout<<"positions x:   "<<positions[n].x<<endl;
			delPos2 = readLine.find(" ", delPos1+1);
			positions[n].y = atof(readLine.substr(delPos1,delPos2 ).c_str());
			//cout<<"positions y:   "<<positions[n].y<<endl;
			delPos3 = readLine.find(" ", delPos2+1);
			positions[n].z = atof(readLine.substr(delPos2,delPos3 ).c_str());
			//cout<<"positions z:   "<<positions[n].z<<endl;
		}		 

		facades = new facade[nf];
		cout<<"nf:  "<<nf<<endl;
		for (int n = 0; n < 20; ++n)
		{
			getline(in,readLine);
			delPos1 = readLine.find(" ",0);
			delPos2 = readLine.find(" ",delPos1+1);
			facades[n].v1 = atoi(readLine.substr(delPos1,delPos2 ).c_str());
			delPos3 = readLine.find(" ",delPos2+1);
			facades[n].v2 = atoi(readLine.substr(delPos2,delPos3 ).c_str());
			delPos4 = readLine.find(" ",delPos3+1);
			facades[n].v3 = atoi(readLine.substr(delPos3,delPos4 ).c_str());
			cout<<"facades v1 : "<<facades[n].v1<<endl;
			cout<<"facades v2 : "<<facades[n].v2<<endl;
		}




	}



};

model* g_model;

void init(){
    floor();     

    //model g_model("386.off");

    g_model = new model( "386.off" );

    cout << "m1.nv: " << g_model->nv << endl;
    cout << "m1.nf: " << g_model->nf << endl;

    // Create and initialize a vertex buffer object for floor, to be used in display()
    glGenBuffers(1, &floor_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, floor_buffer);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(floor_points),NULL, GL_STATIC_DRAW);
    //glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(floor_points), floor_points);
    glBufferData(GL_ARRAY_BUFFER, sizeof(point3) * g_model->nv,g_model->positions, GL_STATIC_DRAW);
    //glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3) * g_model->nv, g_model->positions);



    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(facade) * g_model->nf, g_model->facades, GL_STATIC_DRAW); 
    //glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeof(g_model->facades), g_model->facades);


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

    cout<<"draw.."<<endl;

    glDisableVertexAttribArray(vPosition);

    //glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
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
    //drawObj(floor_buffer, floor_NumVertices);  // draw the floor
    //drawObj(floor_buffer, 2087);
    drawObj(floor_buffer, EBO, g_model->nf);
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
    /*switch(key) {
	    case 033: // Escape Key
	    case 'X': eye[0] += 1.0;  glUniform1f( glGetUniformLocation(program, "_pasoX"), pasoX ); pasoX+=1.0;  break;// _cameraDir[0] -= 1.0; break;
	    case 'x': eye[0] -= 1.0;  glUniform1f( glGetUniformLocation(program, "_pasoX"), pasoX ); pasoX-=1.0;  break;// _cameraDir[0] += 1.0; break;
	    //case 'Y': eye[1] += 1.0;  glUniform1f( glGetUniformLocation(program2, "_pasoY"), pasoZ ); pasoY+=1.0;  break;// _cameraDir[1] -= 1.0; break;
	    //case 'y': eye[1] -= 1.0;  if(eye[1]<1) {eye[1] += 1.0;pasoY+=1.0;} glUniform1f( glGetUniformLocation(program2, "_pasoY"), pasoZ ); pasoY-=1.0;  break;// _cameraDir[1] += 1.0; break;
	    case 'Z': eye[2] += 1.0;  glUniform1f( glGetUniformLocation(program, "_pasoZ"), pasoZ ); pasoZ+=1.0;  break;// _cameraDir[2] -= 1.0; break;
	    case 'z': eye[2] -= 1.0;  glUniform1f( glGetUniformLocation(program, "_pasoZ"), pasoZ ); pasoZ-=1.0;  break;// _cameraDir[2] += 1.0; break;
    }*/
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