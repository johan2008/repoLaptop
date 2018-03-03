//g++ -std=c++11 mainFile.cpp InitShader.cpp model.cpp axisExtraction.cpp axisAlignment.cpp -o mainFile -lglut -lGLEW -lGL -lGLU && ./mainFile
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
#include "axisAlignment.h"
#include <map>
using namespace std;

#define CLEAR_COLOR 0.529f, 0.807f, 0.92f, 1.0f
#define PI 3.141592654

double longitud_curve(vector< point3 > positions);
double avg_curve(vector<point3> positions);

typedef Angel::vec3  color3;
typedef Angel::vec3  point3;
typedef Angel::vec2  point2;

point3 centroide_model1;
point3 centroide_model2;

vector<vector<point3>> curves_0_levels1;
vector<vector<point3>> curves_0_levels2;

int windowWidth = 600;
int windowHeight = 600;
GLuint program;
GLuint model_buffer;
GLuint model_buffer2;
GLuint line_buffer1;
GLuint line_buffer2;
GLuint line_buffer3;
GLuint cero_level_buffer1[10];
GLuint cero_level_buffer2[10];
GLuint EBO;
GLuint EBO2;

/*
struct curveStruct{
    vector<point3> curve;
    double q;
};
*/

vector<curveStruct> model_curve1;
vector<curveStruct> model_curve2;


GLfloat angle = 0.0;
GLuint  u_tView;
GLfloat  fovy = 45.0;

vec3 init_eye(5.0, 3.0, 0.0);
vec3 eye = init_eye;
vec3 _cameraDir( -5.0f, -3.0f, 0.0f );

GLfloat  aspect; 

const int floor_NumVertices = 6; //(1 face)*(2 triangles/face)*(3 vertices/triangle)
point3 floor_points[floor_NumVertices]; // positions for all vertices
color3 floor_colors[floor_NumVertices]; // colors for all vertices
point3 floor_normal[floor_NumVertices];

float pasoX = 0;
float pasoZ = 0;
float pasoY = 0;

vector<color3> _colors;
vector<color3> _colors2;
vector<color3> _colors_cero_level;
    
model* g_model;
model* g_model2;


//BLUE AXI-Z
point3 Line1[2] = {
    point3( 0.0, 0.0,  0.0),
    point3( 00.0, 0.0, 10.0)
};
color3 vertex_colors_Line1[2] = {
    color3( 0.0, 0.0, 0.0),  // 
    color3( 0.0, 0.0, 0.0)  // 
};

//RED AXI-X
point3 Line2[2] = {
    point3( 0.0, 0.0,  0.0),
    point3( 2.0, 0.0, 00.0)
};
color3 vertex_colors_Line2[2] = {
    color3( 0.0, 0.0, 0.0),  // 
    color3( 0.0, 0.0, 0.0)  // 
};

//GREEN AXI-Y
point3 Line3[2] = {
    point3( 0.0 , 0.0,  0.0),
    point3( 00.0,20.0, 0.0)
};
color3 vertex_colors_Line3[2] = {
    color3( 0.0, 0.0, 0.0),  // black
    color3( 0.0, 0.0, 0.0)  // 
};

double distancePlano(){

}

void colorSurface(model* g_model, vector<color3> &_colors){
    for(int n=0; n< g_model->nv; n++)
    {
        if( abs(g_model->positions[n].z) <0.0085 ){
            //cout<<"entro"<<endl;
            //<<g_model->positions[n]<<endl;
            point3 color_element;
            color_element.x = 0.0;
            color_element.y = 0.0;
            color_element.z = 1.0;
            _colors.push_back(color_element);
        }else if(abs(g_model->positions[n].z) <0.082){
            point3 color_element;
            color_element.x = 0.0;
            color_element.y = 0.0;
            color_element.z = 0.0;
            _colors.push_back(color_element);
        }else if(abs(g_model->positions[n].z) <0.182){
            point3 color_element;
            color_element.x = 0.0;
            color_element.y = 0.0;
            color_element.z = 0.0;
            _colors.push_back(color_element);
        }
        else{
            point3 color_element;
            color_element.x = 0.0;
            color_element.y = 0.0;
            color_element.z = 0.0 ;
            _colors.push_back(color_element);   
        }


            
    }
}


void colorAxi(vector<color3> positions, vector<color3> &_colors){
    for(int n=0; n< positions.size(); n++)
    {
        //if( abs(g_model->positions[n].z) <0.0085 ){
            //cout<<"entro"<<endl;
            //<<g_model->positions[n]<<endl;
            point3 color_element;
            color_element.x = 255.0;
            color_element.y = 255.0;
            color_element.z = 0.0;
            _colors.push_back(color_element);
            cout<<"color"<<endl;
        //}
        /*else if(abs(g_model->positions[n].z) <0.082){
            point3 color_element;
            color_element.x = 1.0;
            color_element.y = 0.0;
            color_element.z = 0.0;
            _colors.push_back(color_element);
        }else if(abs(g_model->positions[n].z) <0.182){
            point3 color_element;
            color_element.x = 0.0;
            color_element.y = 1.0;
            color_element.z = 0.0;
            _colors.push_back(color_element);
        }
        else{
            point3 color_element;
            color_element.x = 0.0;
            color_element.y = 0.0;
            color_element.z = 1.0 ;
            _colors.push_back(color_element);   
        }
        */

            
    }
}


vector<point3> positions;
color3 *colors;


double qCurve(vector<point3> positions){
    return longitud_curve(positions)*(avg_curve(positions) -0.5)/0.5;
}


void position_cero_level(model* g_model, vector<curveStruct> &model_curve ){
    int i=0;
    int count = 0;
    int angle=90;        //45° for example 
    
    for (int i = 0; i < 10; )
    {

        vector<point3> positions_curve;
        vector<color3> colors_curve;

        angle -=i*5;
        double landa = 0;
        float radians=(2*3.1416*angle)/360; 
        float cosine=(float)cos(radians); 
        float sine=(float)sin(radians);

        for(int n=0; n< g_model->nv; n++)
        {
            if( abs( g_model->positions[n].z*cosine + g_model->positions[n].z*sine ) <0.0085 ){
                count++;
                landa +=abs( g_model->positions[n].z*cosine + g_model->positions[n].z*sine ); 
                point3 pos;
                pos.x = g_model->positions[n].x;
                pos.y = g_model->positions[n].y;
                pos.z = g_model->positions[n].z;
                positions.push_back(pos);
                positions_curve.push_back(pos);
                //i++;
            }  
        }
        for(int n=0; n< positions.size(); n++)
        {
            point3 color_element;
            color_element.x = 255;
            color_element.y = 255;
            color_element.z = 0;
            colors_curve.push_back(color_element);
        }
        //cout<<"count:  "<<count<<endl;
        //cout<<" landa final: " <<landa/count<<endl;
        landa = 0;
        count = 0;
        
        //curves_0_levels.push_back(positions);
        //cout<<"position size:  "<<positions_curve.size()<<endl;
        if(positions_curve.size()<200) continue;
        else i++;
        double q = qCurve(positions_curve);
        curveStruct c = {positions_curve,colors_curve,q};
        model_curve.push_back(c);
    }

}


double distance_function(point3 a, point3 b){
    return sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y)  + (a.z-b.z)*(a.z-b.z)  );
}


double longitud_curve(vector< point3 > positions)
{
    double distance_curve = 0;
    for (int i = 0; i < positions.size()-1; ++i)
    {
        distance_curve += distance_function(positions[i], positions[i+1] );
    }
    return distance_curve;
}


double avg_curve(vector<point3> positions)
{
    double distance = 0;
    for (int i = 0; i < positions.size(); ++i)
    {
        distance += abs(positions[i].z);
    }
    return distance;
}


double qAxis(){
    return longitud_curve(positions)*(avg_curve(positions) -0.5)/0.5;
}


void init(){
    g_model  = new model( "../Data/SCAPE/Meshes/mesh000.off" , 0, 1);
    g_model2 = new model( "../Data/SCAPE/Meshes/mesh001.off" , 1, 1); 

    //g_model  = new model( "cat0.off" , 0);
    //g_model  = new model( "../Data/TOSCA/Meshes/wolf0.off" , 0,150);
    //g_model2 = new model( "../Data/TOSCA/Meshes/wolf0.off" , 1,150); 



    axisExtraction axisE1;// = new axisExtraction();
    axisAlignment axisAl;
    axisE1.printData(g_model);  
    axisE1.getCentroide(g_model);   
    axisE1.averageDistance(g_model);
    axisE1.cero_level(g_model);
    axisE1.create_plane(g_model);

    position_cero_level(g_model,model_curve1);
    position_cero_level(g_model2,model_curve2);

    //cout<<"ccc "<<model_curve2[0].curve[0]<<endl;
    //longitud_curve(positions);
    //cout<<"size1:  "<<model_curve1.size()<<endl;
    //cout<<"size1:  "<<model_curve2.size()<<endl;
    vector<curveStruct> c = axisAl.find_curve(model_curve1, model_curve2);
    //cout<<"axisAlignment:   :: "<<axisAl.find_curve(model_curve1, model_curve2)<<endl;

    //cout<< "axisAlignment::::  " <<c.size()<<endl;
    //cout<< "axisAlignment::::  " <<c[0].curve[0]<<endl;
    //_colors = c[0].


    //cout<<"curve and quality: " << model_curve1[0].q<<endl;
    //cout<<"curve and quality: " << model_curve1[1].q<<endl;
    //cout<<"epsilon:  "<<axisE1.epsilon(g_model)<<endl;
    //cout<<"epsilon:  "<<axisE1.epsilon(g_model2)<<endl;


    centroide_model1 = axisE1.getCentroide(g_model);
    centroide_model2 = axisE1.getCentroide(g_model2);

    //colorAxi(c[0].curve, _colors);
    //colorSurface(g_model,_colors);
    //colorSurface(g_model2, _colors2);

    // Create and initialize a vertex buffer object for floor, to be used in display()
    glGenBuffers(1, &model_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, model_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(point3)*g_model->nv  +sizeof(color3)*_colors.size(),NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3) * g_model->nv, g_model->positions);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(point3) * g_model->nv,g_model->positions, GL_STATIC_DRAW);

    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3) *g_model->nv, sizeof(color3) * _colors.size(),_colors.data());

    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(facade) * g_model->nf, g_model->facades, GL_STATIC_DRAW); 
    //glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeof(g_model->facades), g_model->facades);


    glGenBuffers(1, &model_buffer2);
    glBindBuffer(GL_ARRAY_BUFFER, model_buffer2);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(floor_points),NULL, GL_STATIC_DRAW);

    glBufferData(GL_ARRAY_BUFFER, sizeof(point3) * g_model2->nv+sizeof(color3)*_colors2.size(),NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3) * g_model2->nv, g_model2->positions);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3) *g_model2->nv, sizeof(color3) * _colors2.size(),_colors2.data());

    glGenBuffers(1, &EBO2);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO2);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(facade) * g_model2->nf, g_model2->facades, GL_STATIC_DRAW); 


    //line
    glGenBuffers(1, &line_buffer1);
    glBindBuffer(GL_ARRAY_BUFFER, line_buffer1);
    glBufferData(GL_ARRAY_BUFFER, sizeof(point3)*2 + sizeof(color3)*2,NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3)*2, Line1);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3)*2, sizeof(color3)*2,vertex_colors_Line1);

    glGenBuffers(1, &line_buffer2);
    glBindBuffer(GL_ARRAY_BUFFER, line_buffer2);
    glBufferData(GL_ARRAY_BUFFER, sizeof(point3)*2 + sizeof(color3)*2,NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3)*2, Line2);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3)*2, sizeof(color3)*2,vertex_colors_Line2);

    glGenBuffers(1, &line_buffer3);
    glBindBuffer(GL_ARRAY_BUFFER, line_buffer3);
    glBufferData(GL_ARRAY_BUFFER, sizeof(point3)*2 + sizeof(color3)*2,NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3)*2, Line3);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3)*2, sizeof(color3)*2,vertex_colors_Line3);


/*    glGenBuffers(1, &cero_level_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, cero_level_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(point3)*positions.size() + sizeof(point3)*positions.size() ,NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3)*positions.size(), positions.data());
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3)*positions.size(), sizeof(color3)*positions.size(),_colors_cero_level.data());
*/


    for (int i = 0; i < model_curve1.size(); ++i)
    {
        glGenBuffers(1, &cero_level_buffer1[i]);
        glBindBuffer(GL_ARRAY_BUFFER, cero_level_buffer1[i]);
        cout<<"MODEL CURVE: "<<model_curve1[i].curve.size()<<endl;
        cout<<"MODEL CURVE: "<<model_curve1[i].color.size()<<endl;
        glBufferData(GL_ARRAY_BUFFER, sizeof(point3)*model_curve1[i].curve.size() + sizeof(point3)*model_curve1[i].curve.size() ,NULL, GL_STATIC_DRAW);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3)*model_curve1[i].curve.size(), model_curve1[i].curve.data());
        glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3)*model_curve1[i].curve.size(), sizeof(point3)*model_curve1[i].curve.size(),model_curve1[i].color.data());

    }


    for (int i = 0; i < model_curve2.size(); ++i)
    {
        glGenBuffers(1, &cero_level_buffer2[i]);
        glBindBuffer(GL_ARRAY_BUFFER, cero_level_buffer2[i]);
        //cout<<"MODEL CURVE: "<<model_curve1[i].curve.size()<<endl;
        //cout<<"MODEL CURVE: "<<model_curve1[i].color.size()<<endl;
        glBufferData(GL_ARRAY_BUFFER, sizeof(point3)*model_curve2[i].curve.size() + sizeof(point3)*model_curve2[i].curve.size() ,NULL, GL_STATIC_DRAW);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3)*model_curve2[i].curve.size(), model_curve2[i].curve.data());
        glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3)*model_curve2[i].curve.size(), sizeof(point3)*model_curve2[i].curve.size(),model_curve2[i].color.data());

    }
    


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

    GLuint vColor = glGetAttribLocation(program, "vColor"); 
    glEnableVertexAttribArray(vColor);
    glVertexAttribPointer(vColor, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(sizeof(point3)*g_model->nv  )); 

    glDrawElements(GL_TRIANGLES, numFaces*3, GL_UNSIGNED_INT, 0);

    glDisableVertexAttribArray(vPosition);
    glDisableVertexAttribArray(vColor);
}


void drawObj2(GLuint buffer, int num_vertices){
    //--- Activate the vertex buffer object to be drawn ---//
    glBindBuffer(GL_ARRAY_BUFFER, buffer);

    /*----- Set up vertex attribute arrays for each vertex attribute -----*/
    GLuint vPosition = glGetAttribLocation(program, "vPosition");
    glEnableVertexAttribArray(vPosition);
    glVertexAttribPointer(vPosition, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(0) );

    GLuint vColor = glGetAttribLocation(program, "vColor"); 
    glEnableVertexAttribArray(vColor);
    glVertexAttribPointer(vColor, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(sizeof(point3) * num_vertices) ); 
      // the offset is the (total) size of the previous vertex attribute array(s)

    /* Draw a sequence of geometric objs (triangles) from the vertex buffer
       (using the attributes specified in each enabled vertex attribute array) */
    glDrawArrays(GL_LINE_STRIP, 0, num_vertices);
    //glDrawArrays(GL_LINE_STRIP,0,vector1.size());
    /*--- Disable each vertex attribute array being enabled ---*/
    glDisableVertexAttribArray(vPosition);
    glDisableVertexAttribArray(vColor);
}


void drawCeroLevel(GLuint buffer, int num_vertices){
    //--- Activate the vertex buffer object to be drawn ---//
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    GLuint vPosition = glGetAttribLocation(program, "vPosition");
    glEnableVertexAttribArray(vPosition);
    glVertexAttribPointer(vPosition, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(0) );

    GLuint vColor = glGetAttribLocation(program, "vColor"); 
    glEnableVertexAttribArray(vColor);
    glVertexAttribPointer(vColor, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(sizeof(point3) * num_vertices) ); 
      // the offset is the (total) size of the previous vertex attribute array(s)

    /* Draw a sequence of geometric objs (triangles) from the vertex buffer
       (using the attributes specified in each enabled vertex attribute array) */
    glDrawArrays(GL_LINE_STRIP, 0, num_vertices);
    //glDrawArrays(GL_LINE_STRIP,0,vector1.size());
    /*--- Disable each vertex attribute array being enabled ---*/
    cout<<"draw cero level"<<endl;
    glDisableVertexAttribArray(vPosition);
    glDisableVertexAttribArray(vColor);
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
    cout<<"CENTROIDE: "<<centroide_model1<<endl;
    mat4 matrixModelModelo= mat4(1.0f);
    matrixModelModelo = Translate(-centroide_model1.x, -centroide_model1.y, -centroide_model1.z) * Scale (1.0, 1.0, 1.0) * Rotate(angle, 0.0, 2.0, 0.0);
    //matrixModelModelo = Translate(-centroide_model1.x, -centroide_model1.y, -centroide_model1.z) * Scale (1.0, 1.0, 1.0) * Rotate(angle, 0.0, 2.0, 0.0);
    //matrixModelModelo = Translate(0, 0, 0) * Scale (1.0, 1.0, 1.0) * Rotate(angle, 0.0, 2.0, 0.0);

    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModelModelo); // GL_TRUE: matrix is row-major
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawObj(model_buffer, EBO, g_model->nf);

    drawObj(model_buffer2, EBO2, g_model2->nf);


    //draw lines
    glUniform1f( glGetUniformLocation(program, "flag"), 0.0 );
    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModel); // GL_TRUE: matrix is row-major
    glUniform1f( glGetUniformLocation(program, "texture_mapped_ground"), 0.0 );
           // Wireframe floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    //drawObj2(line_buffer1, 2);  // draw the floor


    glUniform1f( glGetUniformLocation(program, "flag"), 0.0 );
    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModel); // GL_TRUE: matrix is row-major
    glUniform1f( glGetUniformLocation(program, "texture_mapped_ground"), 0.0 );
            // Wireframe floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    //drawObj2(line_buffer2, 2);  // draw the floor

    glUniform1f( glGetUniformLocation(program, "flag"), 0.0 );
    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModel); // GL_TRUE: matrix is row-major
    glUniform1f( glGetUniformLocation(program, "texture_mapped_ground"), 0.0 );
             // Wireframe floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    //drawObj2(line_buffer3, 2);  // draw the floor


    for (int i = 0; i < model_curve1.size(); ++i)
    {
        glUniform1f( glGetUniformLocation(program, "flag"), 0.0 );
        glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModelModelo); // GL_TRUE: matrix is row-major
        glUniform1f( glGetUniformLocation(program, "texture_mapped_ground"), 0.0 );
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        //cout<<"size:  "<<positions.size()<<endl;
        drawCeroLevel(cero_level_buffer1[i], model_curve1[i].curve.size());
    }

    for (int i = 0; i < model_curve2.size(); ++i)
    {
        glUniform1f( glGetUniformLocation(program, "flag"), 0.0 );
        glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModelModelo); // GL_TRUE: matrix is row-major
        glUniform1f( glGetUniformLocation(program, "texture_mapped_ground"), 0.0 );
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        //cout<<"size:  "<<positions.size()<<endl;
        drawCeroLevel(cero_level_buffer2[i], model_curve2[i].curve.size());
    }

      // 



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
    glutCreateWindow("Final Proyect");
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