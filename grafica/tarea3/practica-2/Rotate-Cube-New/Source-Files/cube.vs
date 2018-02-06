#version 120 
attribute  vec3 vPosition;
attribute  vec3 aNormal;

varying vec3 TexCoords;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;


uniform float _pasoX;
uniform float _pasoZ;
uniform float _pasoY;


varying vec3 Position;
varying vec3 Normal;


void main()
{
    TexCoords = vPosition;
    vec4 pos = projection * view * vec4(vPosition.x + _pasoX, vPosition.y +_pasoY, vPosition.z + _pasoZ,  1.0);

    //Normal = mat3(transpose(inverse(model)))*aNormal;
    //Position = vec3(model*vec4(vPosition,1.0));
    gl_Position = pos.xyww;
}  