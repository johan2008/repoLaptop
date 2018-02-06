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
using namespace std;

#define CLEAR_COLOR 0.529f, 0.807f, 0.92f, 1.0f
#define anglevelocity 2 //angle increment for every frame


static GLuint texName;


#define ImageWidth  64
#define ImageHeight 64
GLubyte Image[ImageHeight][ImageWidth][4];

#define stripeImageWidth 32
GLubyte stripeImage[4*stripeImageWidth];


vec2 quad_texCoord[6] = {
    vec2(1.25, 0.0),  // for d
    vec2(1.25, 1.5),  // for c
    vec2(0.0, 1.5),  // for a 
    vec2(0.0, 0.0),  // for b
    vec2(0.0, 1.5),  // for a
    vec2(1.25, 0.0),  // for c
};


void image_set_up(void){
 int i, j, c; 
 
 /* --- Generate checkerboard image to the image array ---*/
  for (i = 0; i < ImageHeight; i++)
    for (j = 0; j < ImageWidth; j++)
      {
       c = (((i & 0x8) == 0) ^ ((j & 0x8) == 0));

       if (c == 1) /* white */
    {
         c = 255;  
     Image[i][j][0] = (GLubyte) c;
         Image[i][j][1] = (GLubyte) c;
         Image[i][j][2] = (GLubyte) c;
    }
       else  /* green */
    {
         Image[i][j][0] = (GLubyte) 0;
         Image[i][j][1] = (GLubyte) 150;
         Image[i][j][2] = (GLubyte) 0;
    }

       Image[i][j][3] = (GLubyte) 255;
      }

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

/*--- Generate 1D stripe image to array stripeImage[] ---*/
  for (j = 0; j < stripeImageWidth; j++) {
     /* When j <= 4, the color is (255, 0, 0),   i.e., red stripe/line.
        When j > 4,  the color is (255, 255, 0), i.e., yellow remaining texture
      */
    stripeImage[4*j] = (GLubyte)    255;
    stripeImage[4*j+1] = (GLubyte) ((j>4) ? 255 : 0);
    stripeImage[4*j+2] = (GLubyte) 0; 
    stripeImage[4*j+3] = (GLubyte) 255;
  }

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
/*----------- End 1D stripe image ----------------*/

/*--- texture mapping set-up is to be done in 
      init() (set up texture objects),
      display() (activate the texture object to be used, etc.)
      and in shaders.
 ---*/

} /* end function */


typedef Angel::vec3  color3;
typedef Angel::vec3  point3;

typedef Angel::vec4 point4;
typedef Angel::vec4 color4;


GLuint Angel::InitShader(const char* vShaderFile, const char* fShaderFile);

GLuint program;       /* shader program object id */
GLuint cube_buffer;   /* vertex buffer object id for cube */
GLuint sphere_buffer;
GLuint floor_buffer;  /* vertex buffer object id for floor */
GLuint line_buffer1;
GLuint line_buffer2;
GLuint line_buffer3;

GLuint fog_mode = 0;

GLuint  u_tView;

// Projection transformation parameters
GLfloat  fovy = 45.0;  // Field-of-view in Y direction angle (in degrees)
GLfloat  aspect;       // Viewport aspect ratio
GLfloat  zNear = 0.5, zFar = 100.0;

GLfloat angle = 0.0; // rotation angle in degrees
GLfloat angleRotation = 0.0;
//vec4 init_eye(3.0, 2.0, 0.0, 1.0); // initial viewer position
vec3 init_eye(7.0, 3.0, -10.0);
vec3 eye = init_eye;               // current viewer position
vec3 _cameraDir( -7.0f, -3.0f, 10.0f );

int flagShading =0;
int enableLighting = 0;
int windowWidth = 600;
int windowHeight = 600;
int animationFlag = 1; // 1: animation; 0: non-animation. Toggled by key 'a' or 'A'
int texture_mapped_ground = 0;

int moveFlag = 0;

int cubeFlag = 1;   // 1: solid cube; 0: wireframe cube. Toggled by key 'c' or 'C'
int floorFlag = 1;  // 1: solid floor; 0: wireframe floor. Toggled by key 'f' or 'F'
int sphereFlag = 1;
int endflag = 0;//flag for the ending of a mode
int lightSource = 0;

float newx=0; // translate x, after calculation d
float newz=0;//translate z, after calculating  d

float angle2 =0;

float pointbeginx ;//start point x
float pointbeginz;//start point z
float pointendx;//end point x
float pointendz;//end point z

float anglespeed=0;//angle increment for every frame. Starts off at 0, until it is set to the anglevelocity. They are the same thing.


float vec[3];//vector for line segment
float vecnorm[3];//vector for normal to line segment

float rotaxis[3];//rotation axis calculate from cross sproduct of the 2 above vectors


int mode = 0;//flag for mode 0
int angleflag = 0;//flag to initialize angle to 0
float rad = 57.2957795;



mat4 rot(1);

int smooth = 0;

int flagSmooth = 0;
const int cube_NumVertices = 36; //(6 faces)*(2 triangles/face)*(3 vertices/triangle)
point3 cube_points[100]; 
color3 cube_colors[100];

const int floor_NumVertices = 6; //(1 face)*(2 triangles/face)*(3 vertices/triangle)
point3 floor_points[floor_NumVertices]; // positions for all vertices
color3 floor_colors[floor_NumVertices]; // colors for all vertices

point3 floor_normal[floor_NumVertices]; 

point3 sphere_points[24]; // positions for all vertices
color3 sphere_colors[24]; // colors for all vertices


double scalar_product(vector<double> a, vector<double> b){
    double product = 0;
    for (int i = 0; i <= a.size()-1; i++)
        for (int i = 0; i <= b.size()-1; i++)
            product = product + (a[i])*(b[i]);
    return product;
}


vector<point3> _vertices;
vector<color3> _colors;
vector<color3> _normas;

point3 vertices[8] = {
    point3( -0.5, -0.5,  0.5),
    point3( -0.5,  0.5,  0.5),
    point3(  0.5,  0.5,  0.5),
    point3(  0.5, -0.5,  0.5),
    point3( -0.5, -0.5, -0.5),
    point3( -0.5,  0.5, -0.5),
    point3(  0.5,  0.5, -0.5),
    point3(  0.5, -0.5, -0.5)
};


color3 vertex_colors[8] = {
    color3( 0.0, 0.0, 0.0),  // black
    color3( 1.0, 0.0, 0.0),  // red
    color3( 1.0, 1.0, 0.0),  // yellow
    color3( 0.0, 1.0, 0.0),  // green
    color3( 0.0, 0.0, 1.0),  // blue
    color3( 1.0, 0.0, 1.0),  // magenta
    color3( 1.0, 1.0, 1.0),  // white
    color3( 0.0, 1.0, 1.0)   // cyan
};

//blue
point3 Line1[2] = {
    point3( 0.0, 0.0,  0.0),
    point3( 00.0, 0.0, 10.0)
};
color3 vertex_colors_Line1[2] = {
    color3( 0.0, 0.0, 1.0),  // black
    color3( 0.0, 0.0, 1.0)  // red
};


point3 Line2[2] = {
    point3( 0.0, 0.0,  0.0),
    point3( 20.0, 0.0, 00.0)
};
color3 vertex_colors_Line2[2] = {
    color3( 1.0, 0.0, 0.0),  // black
    color3( 1.0, 0.0, 0.0)  // red
};


point3 Line3[2] = {
    point3( 0.0 , 0.0,  0.0),
    point3( 00.0,20.0, 0.0)
};
color3 vertex_colors_Line3[2] = {
    color3( 1.0, 0.0, 1.0),  // black
    color3( 1.0, 0.0, 1.0)  // red
};


int Index = 0; // YJC: This must be a global variable since quad() is called
               //      multiple times and Index should then go up to 36 for
               //      the 36 vertices and colors

// quad(): generate two triangles for each face and assign colors to the vertices
void quad( int a, int b, int c, int d )
{
    cube_colors[Index] = vertex_colors[a]; cube_points[Index] = vertices[a]; Index++;
    cube_colors[Index] = vertex_colors[b]; cube_points[Index] = vertices[b]; Index++;
    cube_colors[Index] = vertex_colors[c]; cube_points[Index] = vertices[c]; Index++;

    cube_colors[Index] = vertex_colors[c]; cube_points[Index] = vertices[c]; Index++;
    cube_colors[Index] = vertex_colors[d]; cube_points[Index] = vertices[d]; Index++;
    cube_colors[Index] = vertex_colors[a]; cube_points[Index] = vertices[a]; Index++;
}


void colorcube(){
    quad( 1, 0, 3, 2 );
    quad( 2, 3, 7, 6 );
    quad( 3, 0, 4, 7 );
    quad( 6, 5, 1, 2 );
    quad( 4, 5, 6, 7 );
    quad( 5, 4, 0, 1 );
}


//-------------------------------
// generate 2 triangles: 6 vertices and 6 colors
void floor()
{
    floor_colors[0] = color3(0,1,0); floor_points[0] = point3(5,0,8);   floor_normal[0] = point3(0,1,0);
    floor_colors[1] = color3(0,1,0); floor_points[1] = point3(5,0,-4);  floor_normal[1] = point3(0,1,0);
    floor_colors[2] = color3(0,1,0); floor_points[2] = point3(-5,0,-4); floor_normal[2] = point3(0,1,0);

    floor_colors[3] = color3(0,1,0); floor_points[3] = point3(-5,0,8);  floor_normal[3] = point3(0,1,0);
    floor_colors[4] = color3(0,1,0); floor_points[4] = point3(-5,0,-4); floor_normal[4] = point3(0,1,0);
    floor_colors[5] = color3(0,1,0); floor_points[5] = point3(5,0,8);   floor_normal[5] = point3(0,1,0);
}


void setStartPoint(float x, float y, float z){
    pointbeginx=x;
    pointbeginz=z;
}

void setEndPoint(float x, float y, float z){
    pointendx=x;
    pointendz=z;
}


float getAngle(float beginx,float beginy, float endx, float endy){
    float inter = (endy-beginy)/(endx-beginx);
    return atan(inter)* rad ;
}  


void createFromFile( string filename ){
    ifstream _fileHandle(filename);
    string _line;
    int _numFaces ,_numVertices;
    _fileHandle>>_numFaces;
    int _indxCounter = 0;
    for ( int q = 0; q < _numFaces; q++ )
    {
        _fileHandle >>_numVertices;
        for ( int p = 0; p < _numVertices; p++ )
        {
            getline( _fileHandle, _line );
            point3 _elemns;// = _split( _line );
            _fileHandle >>_elemns.x>>_elemns.y>>_elemns.z;
            _vertices.push_back(_elemns);
        }
        int size = _vertices.size();
        point3 normVertex = normalize( cross(_vertices[size-2] - _vertices[size-3], _vertices[size-1] - _vertices[size-2]) );
        _normas.push_back( normVertex );
        _normas.push_back( normVertex );
        _normas.push_back( normVertex );


    }
    for(int n=0; n<_vertices.size(); n++)
    {
        point3 color_element;
        color_element.x = 1.0;
        color_element.y = 0.84;
        color_element.z = 0.0;
        _colors.push_back(color_element);
    }
}


// OpenGL initialization
void init(){
    colorcube();

    image_set_up();
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glGenTextures(1, &texName);      // Generate texture obj name(s)

    glActiveTexture( GL_TEXTURE0 );  // Set the active texture unit to be 0 
    glBindTexture(GL_TEXTURE_2D, texName); // Bind the texture to this texture unit

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, ImageWidth, ImageHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, Image);


    //sphere
    glGenBuffers(1,&sphere_buffer);
    glBindBuffer(GL_ARRAY_BUFFER,sphere_buffer);

    glBufferData(GL_ARRAY_BUFFER, sizeof(point3)*_vertices.size() + sizeof(color3)*_vertices.size() +sizeof(color3)*_vertices.size(),NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3) * _vertices.size(), _vertices.data());
    //glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3) * _vertices.size(), sizeof(color3) * _vertices.size(),vertex_colors_buffer);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3) * _vertices.size(), sizeof(color3) * _vertices.size(),_colors.data());

    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3) * _vertices.size()+ sizeof(point3) * _vertices.size(), sizeof(point3) * _vertices.size(),_normas.data());


    



    floor();     
    // Create and initialize a vertex buffer object for floor, to be used in display()
    glGenBuffers(1, &floor_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, floor_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(floor_points) + sizeof(floor_colors) +  sizeof(floor_normal)+ sizeof(quad_texCoord),NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(floor_points), floor_points);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(floor_points) , sizeof(floor_colors),floor_colors);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(floor_colors) + sizeof(floor_points), sizeof(floor_normal),floor_normal);

    glBufferSubData(GL_ARRAY_BUFFER, sizeof(floor_colors) + sizeof(floor_points) + sizeof(floor_normal), sizeof(quad_texCoord), quad_texCoord);


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

    // Load shaders and create a shader program (to be used in display())
    program = InitShader("../Shaders/vshader42.glsl", "../Shaders/fshader42.glsl");
    
    glEnable( GL_DEPTH_TEST );
    glClearColor( CLEAR_COLOR); 
    glLineWidth(2.0);
}
//----------------------------------------------------------------------------
// drawObj(buffer, num_vertices):
//   draw the object that is associated with the vertex buffer object "buffer"
//   and has "num_vertices" vertices.
//
void drawObj(GLuint buffer, int num_vertices){
    color4 light_global( 1.0, 1.0, 1.0, 1.0 );

    point4 light_position( -0.638813, -0.238667, 0.738211 , 0.0);

    //mat3 viewRot = mat3(u_tView);
    //mat3 invTransp = transpose(viewRot);
    //point4 light_direction = invTransp*light_position;


    color4 light_ambient( 0.0, 0.0, 0.0, 1.0 );
    color4 light_diffuse( 0.8, 0.8, 0.8, 1.0 );
    color4 light_specular( 0.2, 0.2, 0.2, 1.0 );
    //plane material
    color4 ambientColorPlane( 0.2, 0.2, 0.2, 1.0 );
    color4 difusaColorPlane( 0.0, 1.0, 0.0, 1.0 );
    color4 specularColorPlane( 0.0, 0.0, 0.0, 1.0 );
    float  material_shininess = 125.0;

    color4 Ambient_product = light_ambient * ambientColorPlane + light_global * ambientColorPlane;
    color4 Diffuse_product = light_diffuse * difusaColorPlane;
    color4 Specular_product = light_specular * specularColorPlane;


    glUniform4fv(glGetUniformLocation(program,"AmbientProduct" ),1 , Ambient_product);
    glUniform4fv(glGetUniformLocation(program,"DiffuseProduct" ),1 , Diffuse_product);
    glUniform4fv(glGetUniformLocation(program,"SpecularProduct" ),1 , Specular_product);
    glUniform1f( glGetUniformLocation(program, "Shininess"), material_shininess );


    glUniform3fv( glGetUniformLocation(program, "fL"),1,  light_position);

    //--- Activate the vertex buffer object to be drawn ---//
    glBindBuffer(GL_ARRAY_BUFFER, buffer);

    /*----- Set up vertex attribute arrays for each vertex attribute -----*/
    GLuint vPosition = glGetAttribLocation(program, "vPosition");
    glEnableVertexAttribArray(vPosition);
    glVertexAttribPointer(vPosition, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(0) );

    GLuint vColor = glGetAttribLocation(program, "vColor"); 
    glEnableVertexAttribArray(vColor);
    glVertexAttribPointer(vColor, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(sizeof(point3)*num_vertices  )); 


    
    GLuint vNormal = glGetAttribLocation(program, "aNormal"); 
    glEnableVertexAttribArray(vNormal);
    glVertexAttribPointer(vNormal, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(sizeof(point3) * num_vertices + sizeof(color3)*num_vertices )  ); 
    

    GLuint vTexCoord = glGetAttribLocation( program, "vTexCoord" ); 
    glEnableVertexAttribArray( vTexCoord );
    glVertexAttribPointer( vTexCoord, 2, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET( sizeof(point3)*num_vertices + sizeof(color3)*num_vertices +sizeof(color3)*num_vertices) ); 


      // the offset is the (total) size of the previous vertex attribute array(s)
    glDrawArrays(GL_TRIANGLES, 0, num_vertices);

    glDisableVertexAttribArray(vPosition);
    glDisableVertexAttribArray(vColor);
    glDisableVertexAttribArray(vNormal);
    glDisableVertexAttribArray(vTexCoord);

}



void drawObjSphera(GLuint buffer, int num_vertices){
    //--- Activate the vertex buffer object to be drawn ---//

    point4 light_position( 0.1, 0.0, -1.0, 0.0 );
    color4 light_ambient( 0.0, 0.0, 0.0, 1.0 );
    color4 light_diffuse( 0.8, 0.8, 0.8, 1.0 );
    color4 light_specular( 0.2, 0.2, 0.2, 1.0 );

    //Sphere material
    color4 ambientColorSphere( 0.2, 0.2, 0.2, 1.0 );
    color4 difusaColorSphere( 1.0, 0.84, 0.0, 1.0 );
    color4 specularColorSphere( 1.0, 0.84, 0.0, 1.0 );
    float  material_shininess = 125.0;

    color4 Ambient_product = light_ambient * ambientColorSphere;
    color4 Diffuse_product = light_diffuse * difusaColorSphere;
    color4 Specular_product = light_specular * specularColorSphere;

    glUniform4fv(glGetUniformLocation(program,"AmbientProduct" ),1 , Ambient_product);
    glUniform4fv(glGetUniformLocation(program,"DiffuseProduct" ),1 , Diffuse_product);
    glUniform4fv(glGetUniformLocation(program,"SpecularProduct" ),1 , Specular_product);
    glUniform4fv( glGetUniformLocation(program, "LightPosition"),1, light_position );
    glUniform1f( glGetUniformLocation(program, "Shininess"), material_shininess );


    glBindBuffer(GL_ARRAY_BUFFER, buffer);

    /*----- Set up vertex attribute arrays for each vertex attribute -----*/
    GLuint vPosition = glGetAttribLocation(program, "vPosition");
    glEnableVertexAttribArray(vPosition);
    glVertexAttribPointer(vPosition, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(0) );

    GLuint vColor = glGetAttribLocation(program, "vColor"); 
    glEnableVertexAttribArray(vColor);
    glVertexAttribPointer(vColor, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(sizeof(point3) * num_vertices) ); 


    GLuint vNormal = glGetAttribLocation(program, "aNormal"); 
    glEnableVertexAttribArray(vNormal);
    glVertexAttribPointer(vNormal, 3, GL_FLOAT, GL_FALSE, 0,BUFFER_OFFSET(sizeof(point3) * num_vertices + sizeof(color3)*num_vertices )  ); 
    

    glDrawArrays(GL_TRIANGLES, 0, num_vertices);

    /*--- Disable each vertex attribute array being enabled ---*/
    glDisableVertexAttribArray(vPosition);
    glDisableVertexAttribArray(vColor);
    glDisableVertexAttribArray(vNormal);
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



void display( void ){

    GLuint  u_tModel;  
    GLuint  u_tProj;  
    u_tView;

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glUseProgram(program); // Use the shader program


    if(flagShading == 0){
        glUniform1f( glGetUniformLocation(program, "flag"), 0.0 );
    }else{
        glUniform1f( glGetUniformLocation(program, "flag"), 3.0 );
    }

    if(texture_mapped_ground == 0){
        glUniform1f( glGetUniformLocation(program, "texture_mapped_ground"), 0.0 );
    }else{
        glUniform1f( glGetUniformLocation(program, "texture_mapped_ground"), 3.0 );
    }



    if(flagSmooth == 0){
        glUniform1f( glGetUniformLocation(program, "flagS"), 0.0 );
    }else{
        glUniform1f( glGetUniformLocation(program, "flagS"), 3.0 );
    }



    glUniform1i( glGetUniformLocation(program, "texture_2D"), 0 );
    

    u_tModel = glGetUniformLocation( program, "u_tModel" );
    u_tProj =  glGetUniformLocation( program, "u_tProj" );
    u_tView =  glGetUniformLocation( program, "u_tView" );

    mat4  matrixProj = Perspective(fovy, 1, 0.1f, 100.0f );

    glUniformMatrix4fv(u_tProj, 1, GL_TRUE, matrixProj); // GL_TRUE: matrix is row-major



    color4 fogColor = (0.7, 0.7, 0.7, 0.5);


    glUniform1i( glGetUniformLocation(program, "fog_mode"), fog_mode);

    glUniform4fv(glGetUniformLocation(program,"FogColor" ),1 , fogColor);
    glUniform4fv(glGetUniformLocation(program,"CameraEye" ),1 , eye);

    
    vec3 _cameraTarget = eye + _cameraDir;
    vec3 _worldUp( 0.0f, 1.0f, 0.0f );

    mat4  matrixView= LookAt( eye,_cameraTarget,_worldUp );
    glUniformMatrix4fv( u_tView, 1, GL_TRUE, matrixView );


    mat4 matrixModel= mat4(1.0f);
    matrixModel = Translate(0.0, 0.0, 0.0) * Scale (1.0, 1.0, 1.0) * Rotate(angle, 0.0, 2.0, 0.0);




    mat4 matrixModelSphere= mat4(1.0f);
    if(moveFlag == 0){
        matrixModelSphere = Translate(3.0, 1.0, 5.0) * Scale (1.0, 1.0, 1.0) * Rotate(angle, 0.0, 2.0, 0.0);
    }
    else{
        matrixModelSphere = Translate(newx, 1.0, newz) * Scale (1.0, 1.0, 1.0) ;

        
        rot = Rotate(1,rotaxis[0], rotaxis[1], rotaxis[2]) * rot;

        matrixModelSphere = matrixModelSphere*rot;//;*Rotate(angleRotation,rotaxis[0], rotaxis[1], rotaxis[2]);



    }

      //mat3 auxModel = transpose(inverse( mat3(matrixModelSphere)));
      //glUniformMatrix4fv(auxModel, 1, GL_TRUE, auxModel);
    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModelSphere); // GL_TRUE: matrix is row-major
    glUniform1i( glGetUniformLocation(program, "Texture_app_flag"), 0);
    if (sphereFlag == 0) // Filled floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    else              // Wireframe floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawObjSphera(sphere_buffer, _vertices.size());  // draw the sphere



    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModel); // GL_TRUE: matrix is row-major
    //glUniform1i( glGetUniformLocation(program, "texture_2D"), 0 );
    glUniform1i( glGetUniformLocation(program, "Texture_app_flag"), 1);
    if (floorFlag == 1) // Filled floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    else              // Wireframe floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawObj(floor_buffer, floor_NumVertices);  // draw the floor



    //draw lines
    glUniform1f( glGetUniformLocation(program, "flag"), 0.0 );
    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModel); // GL_TRUE: matrix is row-major
    glUniform1f( glGetUniformLocation(program, "texture_mapped_ground"), 0.0 );
    if (floorFlag == 1) // Filled floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    else              // Wireframe floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawObj2(line_buffer1, 2);  // draw the floor


    glUniform1f( glGetUniformLocation(program, "flag"), 0.0 );
    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModel); // GL_TRUE: matrix is row-major
    glUniform1f( glGetUniformLocation(program, "texture_mapped_ground"), 0.0 );
    if (floorFlag == 1) // Filled floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    else              // Wireframe floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawObj2(line_buffer2, 2);  // draw the floor

    glUniform1f( glGetUniformLocation(program, "flag"), 0.0 );
    glUniformMatrix4fv(u_tModel, 1, GL_TRUE, matrixModel); // GL_TRUE: matrix is row-major
    glUniform1f( glGetUniformLocation(program, "texture_mapped_ground"), 0.0 );
    if (floorFlag == 1) // Filled floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    else              // Wireframe floor
       glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawObj2(line_buffer3, 2);  // draw the floor


    glutSwapBuffers();
}


void idle (void){
    //angle += 0.02;//POCO
    if(moveFlag == 0){
    }
    else{    
        angle += 0.0;    //YJC: change this value to adjust the cube rotation speed.
        angleRotation += 2;

        float aux;

        float theta;//angle that is arctan of x and y components
        float d;//distance
        float pi = 3.14159265;



        angle2+=anglespeed =  1 ;//angle increment for every frame
            //reset angle if a new mode is beginning
        if (angleflag==1){
            angle2=0;
        }

        //from a to b
        if (mode ==0){
            setStartPoint(3,1,5);
            setEndPoint(-1,1,-4);
            //setting vector values used to calculate cross product and axis of rotation



            vec[0]=(3- (-1));
            vec[1]=(1-1);
            vec[2]=(5- (-4));

            aux = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

            vec[0]=vec[0]/aux;
            vec[1]=vec[1]/aux;
            vec[2]=vec[2]/aux;



        }
        //from b to c
        else if (mode == 1){
            setStartPoint(-1,1,-4);
            setEndPoint(3.5,1,-2.5);
            //setting vector values used to calculate cross product and axis of rotation
            vec[0]=(-1- 3.5);
            vec[1]=(1-1);
            vec[2]=(-4-(-2.5));
        }
        //from c to a
        else if (mode == 2){
            cout<<"set mode2:  "<<endl;
            setStartPoint(3.5,1,-2.5);
            setEndPoint(3,1,5);
            //setting vector values used to calculate cross product and axis of rotation
            vec[0]=(3.5-3);
            vec[1]=(1-1);
            vec[2]=(-2.5-5);
        }


        vecnorm[0]=0;
        vecnorm[1]=1;
        vecnorm[2]=0;
        //Cross product to find rotation axes
        rotaxis[0]=vecnorm[2]*vec[1]-vecnorm[1]*vec[2];
        rotaxis[1]=vecnorm[2]*vec[0]-vecnorm[0]*vec[2];
        rotaxis[2]=vecnorm[1]*vec[0]-vecnorm[0]*vec[1];

        d = angle2*(2*pi*1)/360;


        if (endflag==1){
            mode+=1;

            cout<<"MODE:  "<< mode <<endl;
            endflag=0;
            if(mode==2) setStartPoint(3.5,1,-2.5);
            if(mode==1) setStartPoint(-1,1,-4);
            //if(mode==0) setStartPoint(3.5,1,-2.5);
            if (mode>2){
            mode=0;
            }
        }


        if (mode==0){
            angleflag=0;
            theta = getAngle(pointbeginx,pointbeginz,pointendx,pointendz);

            cout<<"pointbeginx:  "<<pointbeginx<<endl;
            cout<<"pointbeginz:  "<<pointbeginz<<endl;
            cout<<"angle2:   "<<angle2<<endl;
            cout<<"d*sin(theta*pi/180):  "<<d*sin(theta*pi/180)<<endl;
            newz= pointbeginz-d*sin(theta*pi/180);
            newx=pointbeginx -d*cos(theta*pi/180);

            //end case
            
            cout<<"newz:   "<< newz <<endl;
            cout<<"newx:   "<< newx <<endl;

            cout<<"mode0::::::::::::::::::::::::  "<<mode<<endl;
            if (newx<=-1){
                newx = -1;
                newz = -4;
                endflag =1;
                angleflag =1;
            }
        }

        else if (mode == 1){
            angleflag=0;
            theta = getAngle(pointbeginx,pointbeginz,pointendx,pointendz);
            newz=pointbeginz + d*sin(theta*pi/180);
            newx=pointbeginx + d*cos(theta*pi/180);
                //end case
            cout<<"newz:   "<< newz <<endl;
            cout<<"newx:   "<< newx <<endl;

            cout<<"mode1::::::::::::::::::::::::::::::::::  "<<mode<<endl;
            if (newx>=3.5){
                newx = 3.5;
                newz = -2.5;
                endflag =1;
                angleflag =1;
            }
        }
        //if we are going from c to a
        else if (mode == 2){
            angleflag=0;
            theta = getAngle(pointbeginx,pointbeginz,pointendx,pointendz);
            newz=pointbeginz - d*sin(theta*pi/180);
            newx=pointbeginx - d*cos(theta*pi/180);

            cout<<"newz:   "<< newz <<endl;
            cout<<"newx:   "<< newx <<endl;

            cout<<"mode2:::::::::::::::::::::::::::  "<<mode<<endl;
            if (newx<=3)
            {
                newx=3;
                newz=5;
                endflag=1;
                angleflag=1;
            }
        }
    }

    glutPostRedisplay();
}


void keyboard(unsigned char key, int x, int y){
    switch(key) {
    case 033: // Escape Key

    case 'b': case 'B':
        moveFlag = 1;
        break;
    case 'q': case 'Q':
        exit( EXIT_SUCCESS );
        break;

    case 'X': eye[0] += 1.0; _cameraDir[0] -= 1.0; break;
    case 'x': eye[0] -= 1.0; _cameraDir[0] += 1.0; break;
    case 'Y': eye[1] += 1.0; _cameraDir[1] -= 1.0; break;
    case 'y': eye[1] -= 1.0; _cameraDir[1] += 1.0; break;
    case 'Z': eye[2] += 1.0; _cameraDir[2] -= 1.0; break;
    case 'z': eye[2] -= 1.0; _cameraDir[2] += 1.0; break;

        case 'a': case 'A': // Toggle between animation and non-animation
        animationFlag = 1 -  animationFlag;
            if (animationFlag == 1) glutIdleFunc(idle);
            else                    glutIdleFunc(NULL);
            break;
       
        case 'c': case 'C': // Toggle between filled and wireframe cube
        cubeFlag = 1 -  cubeFlag;   
            break;

        case 'f': case 'F': // Toggle between filled and wireframe floor
        floorFlag = 1 -  floorFlag; 
            break;

    case ' ':  // reset to initial viewer/eye position
        eye = init_eye;
        break;
    }
    glutPostRedisplay();
}


void reshape(int width, int height){
    glViewport(0, 0, width, height);
    aspect = (GLfloat) width  / (GLfloat) height;
    glutPostRedisplay();
}


void myMouseFunc(int button, int state, int x, int y){
    if(button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN){
        cout<<"GLUT_RIGHT_BUTTON"<<endl;
        if(moveFlag){
            animationFlag = 1 -  animationFlag;
            if (animationFlag == 1) glutIdleFunc(idle);
            else                    glutIdleFunc(NULL);
        }  
    }
}


void menu(int num){
  if(num == 0){
    glutDestroyWindow(1);
    exit(0);
  }else{
  }
  glutPostRedisplay();
} 


void myMenu(int id){
    switch(id){
    case 1:
        moveFlag = 0;
        angle2=1  ;
        cout<<"case1:   "<<endl;
        glutIdleFunc(idle);
        break;
    case 2:
        exit( EXIT_SUCCESS );
        break;
    case 3:
        break;
    }    
    glutPostRedisplay();
}


void shadowMenu(int id){
    int num =0;
    if(num == 0){
        glutDestroyWindow(1);
        exit(0);
    }else{
    }
    glutPostRedisplay();
}


void lightMenu(int id){
    //sphere.fill_flag = true;

    cout<<"id::::::::    "<<id<<endl;
    switch(id){
    case 1:
        enableLighting = 0;
        flagShading =0;
        break;
    case 2:
        enableLighting = 1;
        flagShading =1;        
        break;
    
    }
}


void wireMenu(int id){
    switch(id){
    case 1:
        sphereFlag = 0;
        break;
    case 2:
        sphereFlag = 1;
        break;    
    }
}


void lightSourceMenu(int id){
    switch(id){
    case 1:
        lightSource = 0;
        break;
    case 2:
        lightSource = 1;
        break;
    
    }
}


void shadingGenerate(){
    glBindBuffer(GL_ARRAY_BUFFER,sphere_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(point3)*_vertices.size() + sizeof(color3)*_vertices.size() +sizeof(color3)*_vertices.size(),NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(point3) * _vertices.size(), _vertices.data());
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3) * _vertices.size(), sizeof(color3) * _vertices.size(),_colors.data());
    if(flagSmooth == 0) glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3) * _vertices.size()+ sizeof(point3) * _vertices.size(), sizeof(point3) * _vertices.size(),_normas.data());
    else  glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3) * _vertices.size()+ sizeof(point3) * _vertices.size(), sizeof(point3) * _vertices.size(),_vertices.data());
}


void shadingMenu(int id){
    switch(id){
    case 1:
        sphereFlag = 0;
        flagSmooth = 0;
        shadingGenerate();
        break;
    case 2:
        sphereFlag = 0;
        flagSmooth = 1;
        shadingGenerate();
        break;    
    }
}


void fogMenu(int id){
    fog_mode = id-1;
    cout<<"fog_mode:  "<<fog_mode <<endl;
}


void textureMappedGroundMenu(int id){
    switch(id){
    case 1:
        texture_mapped_ground = 0;
        break;
    case 2:
        texture_mapped_ground = 1;
        break;
    }
}


void addControl(){
    glutKeyboardFunc(keyboard);
    glutMouseFunc(myMouseFunc);
    GLuint subShadowMenu = glutCreateMenu(shadowMenu);
    glutAddMenuEntry("No",1);
    glutAddMenuEntry("Yes",2);
    GLuint subLightMenu = glutCreateMenu(lightMenu);
    glutAddMenuEntry("No",1);
    glutAddMenuEntry("Yes",2);
    GLuint subWireMenu = glutCreateMenu(wireMenu);
    glutAddMenuEntry("No",1);
    glutAddMenuEntry("Yes",2);
    GLuint subShadingMenu = glutCreateMenu(shadingMenu);
    glutAddMenuEntry("Flat shading",1);
    glutAddMenuEntry("Smooth shading",2);
    GLuint subLSMenu = glutCreateMenu(lightSourceMenu);
    glutAddMenuEntry("Spot Light",1);
    glutAddMenuEntry("Point Source",2);
    GLuint subFogMenu = glutCreateMenu(fogMenu);
    glutAddMenuEntry("No Fog",1);
    glutAddMenuEntry("Linear Fog",2);
    glutAddMenuEntry("Exponential Fog",3);
    glutAddMenuEntry("Exponential Square Fog",4);
    GLuint subGroundTexMenu = glutCreateMenu(textureMappedGroundMenu);
    glutAddMenuEntry("No", 1);
    glutAddMenuEntry("Yes", 2);
    glutCreateMenu(myMenu);
    glutAddMenuEntry("Default View Port",1);
    glutAddSubMenu("Enable Lighting", subLightMenu);
    glutAddSubMenu("Shading", subShadingMenu);
    glutAddSubMenu("Lighting", subLSMenu);
    glutAddSubMenu("Wire Frame", subWireMenu);
    glutAddSubMenu("Fog", subFogMenu);
    glutAddSubMenu("Texture Mapped Ground", subGroundTexMenu);
    glutAddMenuEntry("Quit",2);
    glutAttachMenu(GLUT_LEFT_BUTTON);
} 

//------------------------------------------------------------------------------
int main(int argc, char **argv){ 

    if(argc>1){
        string name = "../image/"+  std::string(argv[1]);
        createFromFile(name);
    }

    int err;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("Color Cube");
    /* Call glewInit() and error checking */
    err = glewInit();
    if (GLEW_OK != err){ 
        printf("Error: glewInit failed: %s\n", (char*) glewGetErrorString(err)); 
        exit(1);
    }

    rotaxis[0] = 1;
    rotaxis[1] = 1;
    rotaxis[2] = 1;
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
    //glutKeyboardFunc(keyboard);
    //glutMouseFunc(myMouseFunc);
    addControl();
    init();
    glutMainLoop();
    return 0;
}
