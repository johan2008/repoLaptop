#version 120
varying vec4 FragColor;

varying vec3 TexCoords;

uniform samplerCube skybox;

varying vec3 Normal;
varying vec3 Position;

uniform vec3 cameraPos;


void main()
{    
	vec3 I = normalize(Position - cameraPos);
	vec3 R = reflect(I, normalize(Normal));
	//gl_FragColor = vec4(textureCube(skybox, R).rgb, 1.0);
    gl_FragColor = textureCube(skybox, TexCoords);

}