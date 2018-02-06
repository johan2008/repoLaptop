#version 120 


attribute vec3 vPosition;
attribute  vec3 vColor;

varying vec4 color;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

uniform vec3 lightColor;



void main(){
  // Output position of the vertex, in clip space : MVP * position
  //gl_Position =  MVP * vec4(vPosition,1);
  // = vec4(vPosition, 1);

  gl_Position = projection * view * model * vec4(vPosition, 1);

  //color = vec4(vColor.r, vColor.g, vColor.b, 1.0);
  color = vec4(lightColor*vec3(vColor.r, vColor.g, vColor.b), 1.0 );


  	//vec3 norm = normalize(Normal);
	//vec3 lightDir = normalize(lightPos - FragPos);  


}
