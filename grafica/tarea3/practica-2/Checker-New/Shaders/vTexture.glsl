/*--------------
Vertex Shader: filename "vTexture.glsl"
---------------*/
 #version 120    // YJC: Comment/un-comment this line to resolve compilation errors
                   //      due to different settings of the default GLSL version

varying  vec3 vPosition;
varying  vec2 vTexCoord;

uniform mat4 ModelView;
uniform mat4 Projection;

uniform vec4 uColor; // obj color (as a uniform variable)

varying vec4 color;
varying vec2 texCoord;

void main() 
{     
 vec4 vPosition4 = vec4(vPosition.x, vPosition.y, vPosition.z, 1.0);
 gl_Position = Projection * ModelView * vPosition4;

 color = uColor;
 texCoord = vTexCoord;

} 
