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
//#include "PerlinNoise.h"
#include "texture.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
using namespace std;

#define CLEAR_COLOR 0.529f, 0.807f, 0.92f, 1.0f

typedef Angel::vec3  color3;
typedef Angel::vec3  point3;
typedef Angel::vec2  point2;

int windowWidth = 600;
int windowHeight = 600;
int num_vertices_terrain;
GLuint program;
GLuint program2;
GLuint floor_buffer;
GLuint GLterrain;


unsigned int skyboxVAO, skyboxVBO;


GLuint GLCube;

static GLuint texName;


GLfloat angle = 0.0;

GLuint  u_tView;
GLfloat  fovy = 45.0;

//vec3 init_eye(7.0, 3.0, -10.0);
vec3 init_eye(0.0, 5.0, 0.0);

vec3 eye = init_eye;
//vec3 _cameraDir( -7.0f, -3.0f, 10.0f );
vec3 _cameraDir( 7.0, 1.0, 0.0 );

GLfloat  aspect; 

const int floor_NumVertices = 6; //(1 face)*(2 triangles/face)*(3 vertices/triangle)
point3 floor_points[floor_NumVertices]; // positions for all vertices
color3 floor_colors[floor_NumVertices]; // colors for all vertices
point3 floor_normal[floor_NumVertices];


const int floor_NumVertices2 = 100; //(1 face)*(2 triangles/face)*(3 vertices/triangle)
point3 floor_points2[floor_NumVertices]; // positions for all vertices
color3 floor_colors2[floor_NumVertices]; // colors for all vertices
point3 floor_normal2[floor_NumVertices];


float pasoX = 0;
float pasoZ = 0;
float pasoY = 0;

vec2 quad_texCoord[6] = {
    vec2(0.0, 0.0),  // for d
    vec2(0.0, 1.0),  // for c
    vec2(1.0, 1.0),  // for a 
    vec2(1.0, 1.0),  // for b
    vec2(1.0, 0.0),  // for a
    vec2(0.0, 0.0),  // for c
};



float skyboxVertices[] = {
        // positions          
        -15.0f,  15.0f, -15.0f,
        -15.0f, -15.0f, -15.0f,
         15.0f, -15.0f, -15.0f,
         15.0f, -15.0f, -15.0f,
         15.0f,  15.0f, -15.0f,
        -15.0f,  15.0f, -15.0f,

        -15.0f, -15.0f,  15.0f,
        -15.0f, -15.0f, -15.0f,
        -15.0f,  15.0f, -15.0f,
        -15.0f,  15.0f, -15.0f,
        -15.0f,  15.0f,  15.0f,
        -15.0f, -15.0f,  15.0f,

         15.0f, -15.0f, -15.0f,
         15.0f, -15.0f,  15.0f,
         15.0f,  15.0f,  15.0f,
         15.0f,  15.0f,  15.0f,
         15.0f,  15.0f, -15.0f,
         15.0f, -15.0f, -15.0f,

        -15.0f, -15.0f,  15.0f,
        -15.0f,  15.0f,  15.0f,
         15.0f,  15.0f,  15.0f,
         15.0f,  15.0f,  15.0f,
         15.0f, -15.0f,  15.0f,
        -15.0f, -15.0f,  15.0f,

        -15.0f,  15.0f, -15.0f,
         15.0f,  15.0f, -15.0f,
         15.0f,  15.0f,  15.0f,
         15.0f,  15.0f,  15.0f,
        -15.0f,  15.0f,  15.0f,
        -15.0f,  15.0f, -15.0f,

        -15.0f, -15.0f, -15.0f,
        -15.0f, -15.0f,  15.0f,
         15.0f, -15.0f, -15.0f,
         15.0f, -15.0f, -15.0f,
        -15.0f, -15.0f,  15.0f,
         15.0f, -15.0f,  15.0f
    };


//floor_points[0] = point(x,y,noise(x,y))


int rows = 100;
int cols = 100;

double terrain[100][100];
double terrain2[100][100];
double flying = 0.1;


struct JTerrain{
	double width;
	double depth;
	int rows;
	int cols;
	point2 pos;

	int i;

	point3 *positions;
	point3 *normals;
	vec2 *texture_terrain;

	JTerrain(double width, double depth, int rows, int cols,point2 pos){
		//cout<<" - "<<pos.x<<endl;
		this->width = width;
		this->depth = depth;
		this->rows  = rows;
		this->cols  = cols;
		this->pos   = pos;
		//this->pos = new point3[6*rows*cols];
		//cambios criticos......................................
		this->positions = new point3[6*rows*cols*16];
		this->texture_terrain = new vec2[6*rows*cols*16];
		unsigned int seed = 237;
		//PerlinNoise pn(seed);


		//cout<<"POS X: "<<rows<<endl;
		//cout<<"POS Y: "<<cols<<endl;

		this->i = 0;
		double auxX , auxY, auxZ;
		for (int z = 0; z < rows; ++z)
		{
			for (int x = 0; x < cols; ++x)
			{
				unsigned int seed = 237;
				//PerlinNoise pn(seed);
				//cout<<" -> " << pos.x+(width/cols)*x <<endl;
				
				auxX = pos.x+(width/cols)*x;
				auxZ = pos.y+(depth/rows)*z;
				positions[ i+0 ] = point3( auxX , 0, auxZ);
				texture_terrain[i+0] = vec2(0.0,1.0); 
    			positions[ i+1 ] = point3( auxX + double(width/cols), 0, auxZ);  
    			texture_terrain[i+1] = vec2(1.0,1.0); 
    			positions[ i+2 ] = point3( auxX , 0, auxZ + (double)(depth/rows) );
    			texture_terrain[i+2] = vec2(0.0,0.0); 
    			positions[ i+3 ] = point3( auxX , 0, auxZ + (double)(depth/rows) );
    			texture_terrain[i+3] = vec2(0.0,0.0); 
    			positions[ i+4 ] = point3( auxX +(double)(width/cols), 0, auxZ + (double)(depth/rows));
    			texture_terrain[i+4] = vec2(1.0,0.0); 
    			positions[ i+5 ] = point3( auxX + (double)(width/cols), 0, auxZ);
				texture_terrain[i+5] = vec2(1.0,1.0); 
				i += 6;
			}
		}

	}



	void addLayer(){
		int i = this->i;
		double auxX , auxY, auxZ;
		double width = this->width*2;
		double depth = this->depth*2;
		int rows = this->rows;
		int cols = this->cols;
		point2 pos(this->pos.x- width/4 , this->pos.y- depth/4 );
		this->pos.x = this->pos.x- width/4;
		this->pos.y = this->pos.y- depth/4;
		//this->positions = new point3[6*rows*cols*4*4];
		unsigned int seed = 237;
		//PerlinNoise pn(seed);


		int u=0;
		for (int z = 0; z < rows; ++z)
		{
			for (int x = 0; x < cols; ++x)
			{
				if( (x >=cols && z >= rows)  && (x < cols*3 && z >= rows)   ) {
					if(z>=rows*3) ;
					else{
						//cout<<"continue......................................."<<endl;
						//cout<<"u "<<++u<< "   z:  "<<z<<" x:  " <<x<<endl;
						continue;
					}
				}
				//cout<<"afuera..............."<<endl;
				//cout<<"u "<<++u<<endl;
				unsigned int seed = 237;
				//PerlinNoise pn(seed);
				//cout<<" -> " << pos.x+(width/cols)*x <<endl;
				//cout<<"break1 . "<<endl;
				auxX = pos.x+(double)(width/cols)*x;
				auxZ = pos.y+(double)(depth/rows)*z;
				//cout<<"break2 . "<< i<<endl;
				positions[ i+0 ] = point3( auxX , 0, auxZ);
				texture_terrain[i+0] = vec2(0.0,1.0);
    			positions[ i+1 ] = point3( auxX + (width/cols), 0, auxZ);  
    			texture_terrain[i+1] = vec2(1.0,1.0);
    			positions[ i+2 ] = point3( auxX , 0, auxZ + (depth/rows) );
    			texture_terrain[i+2] = vec2(0.0,0.0);
    			positions[ i+3 ] = point3( auxX , 0, auxZ + (depth/rows) );
    			texture_terrain[i+3] = vec2(0.0,0.0);
    			positions[ i+4 ] = point3( auxX +(width/cols), 0, auxZ + (depth/rows));
    			texture_terrain[i+4] = vec2(1.0,0.0);
    			positions[ i+5 ] = point3( auxX + (width/cols), 0, auxZ);
				texture_terrain[i+5] = vec2(1.0,1.0); 
				//cout<<"break3 . "<<endl;
				i += 6;
			}
		}
		//this->rows = this->rows*2;
		//this->cols = this->cols*2;
		this->width = this->width*2;
		this->depth = this->depth*2;
		this->i    = i;
	}



};




unsigned int loadCubemap(vector<std::string> faces)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);

    int width, height, nrChannels;
    for (unsigned int i = 0; i < faces.size(); i++)
    {
        unsigned char *data = stbi_load(faces[i].c_str(), &width, &height, &nrChannels, 0);
        if (data)
        {
            glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 
                         0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data
            );
            stbi_image_free(data);
        }
        else
        {
            stbi_image_free(data);
        }
    }
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    return textureID;
}  


void init2(){

	program = InitShader("../Shaders/vshader.glsl", "../Shaders/fshader.glsl");
    program2 = InitShader("cube.vs", "cube.fs");


    glEnable(GL_DEPTH_TEST);

    glGenVertexArrays(1, &skyboxVAO);
    glGenBuffers(1, &skyboxVBO);
    glBindVertexArray(skyboxVAO);
    glBindBuffer(GL_ARRAY_BUFFER, skyboxVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(skyboxVertices), &skyboxVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindVertexArray(0);


	//point2 posInitial(-5,-4);
	point2 posInitial(-16,-16);
	//JTerrain jTerrain(10, 12, 3, 5, posInitial);
	JTerrain jTerrain(32, 32, 64, 64, posInitial);
	num_vertices_terrain = 6*jTerrain.rows*jTerrain.cols*13*13;
	jTerrain.addLayer();
	jTerrain.addLayer();
	jTerrain.addLayer();
	num_vertices_terrain = jTerrain.i; //6*jTerrain.rows*jTerrain.cols;//*4*4;

	

	GLuint floorNormalTexture[6];
	glGenTextures(2, floorNormalTexture);
	glActiveTexture( GL_TEXTURE0 );
	glBindTexture(GL_TEXTURE_2D, floorNormalTexture[0]);


	


	//int width, height, components;
	//unsigned *image = read_texture("image1.jpg", &width, &height, &components);

	int width, height, components;
	unsigned char *image = stbi_load("closeupmossnature.jpg", &width, &height, &components, 0);


	glUniform1i( glGetUniformLocation(program, "texture_2D"), 0 );
    

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
	//free(image);











	//GLuint floorNormalTexture2;
	//glGenTextures(1, &floorNormalTexture2);
	glActiveTexture( GL_TEXTURE1 );
	glBindTexture(GL_TEXTURE_2D, floorNormalTexture[1]);




	int width2, height2, components2;
	unsigned char *image2 = stbi_load("mossnatureground.jpg", &width2, &height2, &components2, 0);
	
	//free(image2);

	glUniform1i( glGetUniformLocation(program, "texture2"), 1 );

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width2, height2, 0, GL_RGB, GL_UNSIGNED_BYTE, image2);



	glActiveTexture( GL_TEXTURE2 );
	glBindTexture(GL_TEXTURE_2D, floorNormalTexture[2]);

	int width3, height3, components3;
	unsigned char *image3 = stbi_load("grassweednature.jpg", &width3, &height3, &components3, 0);
	
	//free(image2);

	glUniform1i( glGetUniformLocation(program, "texture3"), 1 );

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width3, height3, 0, GL_RGB, GL_UNSIGNED_BYTE, image3);





	glActiveTexture( GL_TEXTURE3 );
	glBindTexture(GL_TEXTURE_2D, floorNormalTexture[3]);

	int width4, height4, components4;
	unsigned char *image4 = stbi_load("cobblestonedirtcobble.jpg", &width4, &height4, &components4, 0);
	
	//free(image2);

	glUniform1i( glGetUniformLocation(program, "texture4"), 1 );

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width4, height4, 0, GL_RGB, GL_UNSIGNED_BYTE, image4);









	glActiveTexture( GL_TEXTURE4 );
	glBindTexture(GL_TEXTURE_2D, floorNormalTexture[4]);

	int width5, height5, components5;
	unsigned char *image5 = stbi_load("groundsandfrozen.jpg", &width5, &height5, &components5, 0);
	
	//free(image2);

	glUniform1i( glGetUniformLocation(program, "texture5"), 1 );

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width5, height5, 0, GL_RGB, GL_UNSIGNED_BYTE, image5);




	glActiveTexture( GL_TEXTURE5 );
	glBindTexture(GL_TEXTURE_2D, floorNormalTexture[3]);

	int width6, height6, components6;
	unsigned char *image6 = stbi_load("footprintsnowground.jpg", &width6, &height6, &components6, 0);
	
	//free(image2);

	glUniform1i( glGetUniformLocation(program, "texture6"), 1 );

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width6, height6, 0, GL_RGB, GL_UNSIGNED_BYTE, image6);





	unsigned int programVAO;
	glGenVertexArrays(1, &programVAO);

	glGenBuffers(1, &GLterrain);
	glBindVertexArray(programVAO);
    glBindBuffer(GL_ARRAY_BUFFER, GLterrain);
    cout<<"VALOR I::::::::::::::::::::::::  "<<jTerrain.i<<endl;
    //cout<<"VALOR SIZE:  "<<6*jTerrain.rows*jTerrain.cols*2<<endl;

    
    glBufferData(GL_ARRAY_BUFFER, sizeof(point3)*jTerrain.i + sizeof(vec2)*jTerrain.i ,NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3)*jTerrain.i, jTerrain.positions);

    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3)*jTerrain.i , sizeof(vec2)*jTerrain.i, jTerrain.texture_terrain);


    glEnableVertexAttribArray(1);

    //glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3)*jTerrain.i , sizeof(vec2)*jTerrain.i, jTerrain.texture_terrain);
    //glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3)*jTerrain.i + sizeof(vec2)*jTerrain.i , sizeof(skyboxVertices), &skyboxVertices);




    //glUseProgram(program2);

	

	vector<std::string> faces
	{
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga"
	};
	unsigned int cubemapTexture = loadCubemap(faces);


	/*glUseProgram(program2);
	glUniform1i( glGetUniformLocation(program2, "skybox"), 0 );
	glBindVertexArray(skyboxVAO);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glBindVertexArray(0);
    glDepthFunc(GL_LESS); // set depth function back to default
	*/





    glEnable( GL_DEPTH_TEST );
    glClearColor( CLEAR_COLOR); 
    glLineWidth(2.0);

}




void drawPlane2(GLuint buffer, int num_vertices){

    GLuint vPosition = glGetAttribLocation(program, "vPosition");
    glEnableVertexAttribArray(vPosition);
    glVertexAttribPointer(vPosition, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(0) );

    //GLuint vColor = glGetAttribLocation(program, "vColor"); 
    //glEnableVertexAttribArray(vColor);
    //glVertexAttribPointer(vColor, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(sizeof(point3)*num_vertices  ));



    GLuint vTexCoord = glGetAttribLocation( program, "vTexCoord" ); 
    glEnableVertexAttribArray( vTexCoord );
    glVertexAttribPointer( vTexCoord, 2, GL_FLOAT, GL_TRUE, 0,BUFFER_OFFSET( sizeof(point3)*num_vertices) ); 






    glDrawArrays(GL_TRIANGLES, 0, num_vertices);

    glDisableVertexAttribArray(vPosition);
    //glDisableVertexAttribArray(vColor);
    glDisableVertexAttribArray(vTexCoord);




    	vector<std::string> faces
	{
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga"
	};
	unsigned int cubemapTexture = loadCubemap(faces);


    glUseProgram(program2);
	glUniform1i( glGetUniformLocation(program2, "skybox"), 0 );
	glBindVertexArray(skyboxVAO);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glBindVertexArray(0);
    glDepthFunc(GL_LESS); // set depth function back to default




}




void drawCube(GLuint buffer, int num_vertices){

    GLuint vPosition = glGetAttribLocation(program2, "vPosition");
    glEnableVertexAttribArray(vPosition);
    glVertexAttribPointer(vPosition, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(0) );




    //GLuint vTexCoord = glGetAttribLocation( program, "vTexCoord" ); 
    //glEnableVertexAttribArray( vTexCoord );
    //cout<<"NUMERO DE VERTICES::: "<<num_vertices<<endl;
    //glVertexAttribPointer( vTexCoord, 2, GL_FLOAT, GL_TRUE, 0,BUFFER_OFFSET( sizeof(point3)*num_vertices) ); 

    	vector<std::string> faces
	{
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga",
	    "ame_desert/red_up.tga"
	};
	unsigned int cubemapTexture = loadCubemap(faces);


	glUniform1i( glGetUniformLocation(program2, "skybox"), 0 );
	glBindVertexArray(skyboxVAO);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glBindVertexArray(0);
    glDepthFunc(GL_LESS); // set depth function back to default





    //glDrawArrays(GL_TRIANGLES, 0, num_vertices);

    glDisableVertexAttribArray(vPosition);
    //glDisableVertexAttribArray(vColor);
    //glDisableVertexAttribArray(vTexCoord);








}


void display2(){

	GLuint  u_tModel;  
    GLuint  u_tProj;



    GLuint  u_tModel2;  
    GLuint  u_tProj2;
    GLuint  u_view2;



	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glUseProgram(program);
    

    u_tModel =  glGetUniformLocation( program, "u_tModel" );
    u_tProj  =  glGetUniformLocation( program, "u_tProj" );
    u_tView  =  glGetUniformLocation( program, "u_tView" );




    vec3 _cameraTarget = eye + _cameraDir;
    vec3 _worldUp( 0.0f, 1.0f, 0.0f );

    mat4  matrixView= LookAt( eye,_cameraTarget,_worldUp );
    glUniformMatrix4fv( u_tView, 1, GL_TRUE, matrixView );

    mat4  matrixProj = Perspective(fovy, 1, 0.1f, 100.0f );
    glUniformMatrix4fv(u_tProj, 1, GL_TRUE, matrixProj); 


    mat4 matrixModel= mat4(1.0f);
    matrixModel = Translate(0.0, 0.0, 0.0) * Scale (1.0, 1.0, 1.0) * Rotate(angle, 0.0, 2.0, 0.0);
    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModel);




    //ourShader.use();
    glUniform1i( glGetUniformLocation(program, "texture_2D"), 0 );
    glUniform1i( glGetUniformLocation(program, "texture2"), 1 );
    glUniform1i( glGetUniformLocation(program, "texture3"), 2 );
    glUniform1i( glGetUniformLocation(program, "texture4"), 3 );
    glUniform1i( glGetUniformLocation(program, "texture5"), 4 );
    glUniform1i( glGetUniformLocation(program, "texture6"), 5 );

    glUniform1i( glGetUniformLocation(program, "skybox"), 0 );

    //ourShader.setInt("texture2",1);

    glUniform1f( glGetUniformLocation(program, "_pasoX"), pasoX );
    glUniform1f( glGetUniformLocation(program, "_pasoZ"), pasoZ );

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawPlane2(GLterrain, num_vertices_terrain);





    //cube
    glDepthFunc(GL_LEQUAL);
    glUseProgram(program2);
    u_tModel2 =  glGetUniformLocation( program2, "u_tModel" );
    u_tProj2  =  glGetUniformLocation( program2, "projection" );
    u_view2  =  glGetUniformLocation( program2, "view" );


    mat4  matrixView2= LookAt( eye,_cameraTarget,_worldUp );
    glUniformMatrix4fv( u_tView, 1, GL_TRUE, matrixView2 );

    mat4  matrixProj2 = Perspective(fovy, 1, 0.1f, 100.0f );
    glUniformMatrix4fv(u_tProj, 1, GL_TRUE, matrixProj2); 


    mat4 matrixModel2= mat4(1.0f);
    matrixModel2 = Translate(0.0, 0.0, 0.0) * Scale (1.0, 1.0, 1.0) * Rotate(angle, 0.0, 2.0, 0.0);
    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModel2);
    //drawCube(program2,36);


	glUniform1f( glGetUniformLocation(program2, "_pasoX"), pasoX );
	glUniform1f( glGetUniformLocation(program2, "_pasoZ"), pasoZ );
	glUniform1f( glGetUniformLocation(program2, "_pasoY"), pasoY );


    	vector<std::string> faces
	{
	    "ame_desert/desertsky_rt.tga",
	    "ame_desert/desertsky_lf.tga",
	    "ame_desert/desertsky_dn.tga",
	    "ame_desert/desertsky_bk.tga",
	    "ame_desert/desertsky_bk.tga",
	    "ame_desert/desertsky_ft.tga"
	};
	unsigned int cubemapTexture = loadCubemap(faces);


    //glUseProgram(program2);
	glUniform1i( glGetUniformLocation(program2, "skybox"), 0 );
	glBindVertexArray(skyboxVAO);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glBindVertexArray(0);
    glDepthFunc(GL_LESS); // set depth function back to default



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
	    case 033: // Escape Key
	    case 'X': eye[0] += 1.0;  glUniform1f( glGetUniformLocation(program, "_pasoX"), pasoX ); pasoX+=1.0;  break;// _cameraDir[0] -= 1.0; break;
	    case 'x': eye[0] -= 1.0;  glUniform1f( glGetUniformLocation(program, "_pasoX"), pasoX ); pasoX-=1.0;  break;// _cameraDir[0] += 1.0; break;
	    case 'Y': eye[1] += 1.0;  glUniform1f( glGetUniformLocation(program2, "_pasoY"), pasoZ ); pasoY+=1.0;  break;// _cameraDir[1] -= 1.0; break;
	    case 'y': eye[1] -= 1.0;  if(eye[1]<1) {eye[1] += 1.0;pasoY+=1.0;} glUniform1f( glGetUniformLocation(program2, "_pasoY"), pasoZ ); pasoY-=1.0;  break;// _cameraDir[1] += 1.0; break;
	    case 'Z': eye[2] += 1.0;  glUniform1f( glGetUniformLocation(program, "_pasoZ"), pasoZ ); pasoZ+=1.0;  break;// _cameraDir[2] -= 1.0; break;
	    case 'z': eye[2] -= 1.0;  glUniform1f( glGetUniformLocation(program, "_pasoZ"), pasoZ ); pasoZ-=1.0;  break;// _cameraDir[2] += 1.0; break;
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


    glutDisplayFunc(display2);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);
    //addControl();
    init2();
    glutMainLoop();
    return 0;


}