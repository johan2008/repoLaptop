 #version 120  // YJC: Comment/un-comment this line to resolve compilation errors
                 //      due to different settings of the default GLSL version

uniform mat4 u_tModel;
uniform mat4 u_tView;
uniform mat4 u_tProj;


varying vec4 vPosition4;

attribute  vec3 vPosition;

void main() 
{


    vPosition4 = vec4(vPosition.x, vPosition.y, vPosition.z, 1.0);
    gl_Position = u_tProj * u_tView*u_tModel  *vPosition4;


    //gl_Position = vec4(vPosition.x, vPosition.y, vPosition.z, 1.0);


} 
