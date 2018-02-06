
#include <GL/glew.h>

#include <GL/glut.h>

//#include "LoadShaders.h"

#include <iostream>

#include "Angel-yjc.h"


using namespace std;

#include <fstream>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#include <string>
#include <vector>
#include <fstream>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>

#define CLEAR_COLOR 0.529f, 0.807f, 0.92f, 1.0f
#define BUFFER_OFFSET( offset )   ((GLvoid*) (offset))


GLuint VertexArrayID;
GLuint programID;
GLuint MatrixID ;


glm::mat4 Projection;
glm::mat4 View;
glm::mat4 Model;
glm::mat4 mv;

static int winwidth = 680,winheight = 560;
static int mx,my;
static int flag=0;
static float rotx=0.0f, roty=0.0f;



glm::vec3 init_eye(7.0, 3.0, -10.0);
glm::vec3 eye = init_eye;               // current viewer position
glm::vec3 _cameraDir( -7.0f, -3.0f, 10.0f );


GLfloat  aspect;


GLuint floor_buffer;
glm::vec3 floor_points[6]; // positions for all vertices
glm::vec3 floor_colors[6]; // colors for all vertices


GLuint program; 

void floor()
{
    floor_colors[0] = glm::vec3(0,1,0); floor_points[0] = glm::vec3(5,0,8);   //floor_normal[0] = glm::vec3(0,1,0);
    floor_colors[1] = glm::vec3(0,1,0); floor_points[1] = glm::vec3(5,0,-4);  //floor_normal[1] = glm::vec3(0,1,0);
    floor_colors[2] = glm::vec3(0,1,0); floor_points[2] = glm::vec3(-5,0,-4); //floor_normal[2] = glm::vec3(0,1,0);

    floor_colors[3] = glm::vec3(0,1,0); floor_points[3] = glm::vec3(-5,0,8);  //floor_normal[3] = glm::vec3(0,1,0);
    floor_colors[4] = glm::vec3(0,1,0); floor_points[4] = glm::vec3(-5,0,-4); //floor_normal[4] = glm::vec3(0,1,0);
    floor_colors[5] = glm::vec3(0,1,0); floor_points[5] = glm::vec3(5,0,8);   //floor_normal[5] = glm::vec3(0,1,0);
}


void idle(){


	glutPostRedisplay();
}

void draw(){
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );





	glutSwapBuffers();

}


void mousefunc(int button,int state,int x,int y)
{
	if (button==GLUT_LEFT_BUTTON) {
		if (state==GLUT_DOWN) {
			flag=1;
		} else {
			flag=0;
		}
	}
}

void motionfunc(int x,int y)
{
	if (flag>0) {
		cout << "??? 1" << endl;
		if (flag>1) {
			rotx+=360.0f*(x-mx)/winwidth;
			roty+=360.0f*(y-my)/winheight;

		}

		mx=x;
		my=y;

		draw();

		flag=2;
	}
}

void keyboardfunc(unsigned char key,int x,int y)
{
	if (key=='q' || key==27) exit(0);
}




void init(){
	floor();

	glGenBuffers(1, &floor_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, floor_buffer);
    cout<<"sizeof(floor_points):  "<<sizeof(floor_points)<<endl;
    glBufferData(GL_ARRAY_BUFFER, sizeof(floor_points) + sizeof(floor_colors) ,NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(floor_points), floor_points);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(floor_points) , sizeof(floor_colors),floor_colors);

    program = InitShader( "vertexShader.vsh", "fragmentShader.fsh" );
    glEnable( GL_DEPTH_TEST );
    glClearColor( CLEAR_COLOR); 
    glLineWidth(2.0);

}

void drawObj(GLuint buffer, int num_vertices){
	glBindBuffer(GL_ARRAY_BUFFER, buffer);

	GLuint vPosition = glGetAttribLocation(program, "vPosition");
    glEnableVertexAttribArray(vPosition);
    glVertexAttribPointer(vPosition, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(0) );

    //cout<<"sizeof(glm::vec3)*6 :  "<<sizeof(glm::vec3)*6<<endl;

    glm::vec3 lightColor(1.0f, 1.0f, 1.0f);
	glUniform3fv( glGetUniformLocation(program, "lightColor"),1,  glm::value_ptr(lightColor));


    GLuint vColor = glGetAttribLocation(program, "vColor"); 
    glEnableVertexAttribArray(vColor);
    glVertexAttribPointer(vColor, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(  sizeof(glm::vec3)*6 )); 

    glDrawArrays(GL_TRIANGLES, 0, num_vertices);

    glDisableVertexAttribArray(vPosition);
    glDisableVertexAttribArray(vColor);

    cout<<"drawObj:  "<<endl;

}


void display( void ){
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glUseProgram(program); // Use the shader program




	Projection = glm::perspective(glm::radians(45.0f), 1.0f, 0.1f, 100.0f);

	//glUniformMatrix4fv(glGetUniformLocation(programID, "projection"), 1, GL_FALSE, glm::value_ptr( Projection ));


	/*View = glm::lookAt(
    glm::vec3(4,3,10), // Camera is at (4,3,3), in World Space
    glm::vec3(0,0,0), // and looks at the origin
    glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
    );*/

	//glUniformMatrix4fv(glGetUniformLocation(programID, "projection"), 1, GL_FALSE, glm::value_ptr( View ));


	glm::vec3 _cameraTarget = eye + _cameraDir;
    glm::vec3 _worldUp( 0.0f, 1.0f, 0.0f );


	Model = glm::mat4(1.0f);
	//Model = Translate(0.0, 0.0, 0.0) * Scale (1.0, 1.0, 1.0) * Rotate(0.2, 0.0, 2.0, 0.0);
	View= lookAt( eye,_cameraTarget,_worldUp );

	//mv = Projection * View;


	


	glUniformMatrix4fv(glGetUniformLocation(program, "model"), 1, GL_FALSE, glm::value_ptr( Model ));
	glUniformMatrix4fv(glGetUniformLocation(program, "view"), 1, GL_FALSE, glm::value_ptr( View ));
	glUniformMatrix4fv(glGetUniformLocation(program, "projection"), 1, GL_FALSE, glm::value_ptr( Projection ));



    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    drawObj(floor_buffer, 6);

    glutSwapBuffers();

}

void reshape(int width, int height){
    glViewport(0, 0, width, height);
    aspect = (GLfloat) width  / (GLfloat) height;
    glutPostRedisplay();
}



int main(int argc, char** argv){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE);
	glutInitWindowSize(600,600);
	glutInitWindowPosition(50,50);
	glutCreateWindow("Ejercicio-2");

	
	cout<<"run..."<<endl;

	//glEnableVertexAttribArray(1);

	int err;
	err = glewInit();
    if (GLEW_OK != err){ 
        printf("Error: glewInit failed: %s\n", (char*) glewGetErrorString(err)); 
        exit(1);
    }

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
    init();
	
	//glutDisplayFunc(draw);

	//glutMouseFunc(mousefunc);
   	//glutMotionFunc(motionfunc);
   	//glutKeyboardFunc(keyboardfunc);


	//glutIdleFunc(idle);


	//iniciar();

	glutMainLoop();

}
