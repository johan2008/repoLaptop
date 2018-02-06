/***************************
 * File: vshader42.glsl:
 *   A simple vertex shader.
 *
 * - Vertex attributes (positions & colors) for all vertices are sent
 *   to the GPU via a vertex buffer object created in the OpenGL program.
 *
 * - This vertex shader uses the Model-View and Projection matrices passed
 *   on from the OpenGL program as uniform variables of type mat4.
 ***************************/

 #version 120  // YJC: Comment/un-comment this line to resolve compilation errors
                 //      due to different settings of the default GLSL version





attribute  vec3 vPosition;
attribute  vec3 vColor;
attribute  vec3 aNormal;
attribute  vec2 vTexCoord;



varying vec4 color;
varying vec3 Normal;
varying vec3 FragPos; 

uniform mat4 model_view;
uniform mat4 projection;


uniform mat4 u_tModel;
uniform mat4 u_tView;
uniform mat4 u_tProj;
uniform mat4 auxModel;



varying mat4 u_tModelF;
varying mat4 u_tViewF;


uniform float flagS;

uniform vec4 LightPosition;

varying  vec3 fN;
varying  vec3 fE;
//varying  vec3 fL;

varying vec4 vertex ;

varying vec4 vPosition4;

varying vec2 texCoord;



void main() 
{


    u_tModelF = u_tModel;
    u_tViewF = u_tView;

    vPosition4 = vec4(vPosition.x, vPosition.y, vPosition.z, 1.0);


	//if(flagS>1){
	//	fN =  -vPosition4.xyz;
	//}else{
		fN =  -aNormal;	
	//}


	fN = vec3( u_tModel * vec4( fN, 0.0 ) );
    
    //fE = (u_tModel*vPosition4).xyz;

    fE = (u_tView*u_tModel*vPosition4).xyz;

    //fL = LightPosition.xyz;
    
	//fL = LightPosition.xyz - vPosition4.xyz;
    
    gl_Position = u_tProj * u_tView*u_tModel  *vPosition4;

    color = vec4(vColor.r, vColor.g, vColor.b, 1.0);


    vertex = vPosition4;

    texCoord = vTexCoord;


} 
