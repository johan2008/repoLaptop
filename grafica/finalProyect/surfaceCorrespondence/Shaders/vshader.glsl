 #version 120

attribute  vec3 vPosition;
attribute  vec3 vColor;
attribute  vec2 vTexCoord;


uniform mat4 u_tModel;
uniform mat4 u_tView;
uniform mat4 u_tProj;

varying vec4 vPosition4;
varying vec4 color;
varying vec2 texCoord;

uniform float _pasoX;
uniform float _pasox;
uniform float _pasoY;
uniform float _pasoy;
uniform float _pasoZ;
uniform float _pasoz;



varying vec3 textCoordsCubo;



vec3 permutation(vec3 x) { return mod(((x*34.0)+1.0)*x, 289.0); }

float pnoise(vec2 v){
  const vec4 C = vec4(0.21132, 0.36602,-0.57735, 0.02439);
  vec2 i  = floor(v + dot(v, C.yy) );
  vec2 x0 = v -   i + dot(i, C.xx);
  vec2 i1;
  i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
  vec4 x12 = x0.xyxy + C.xxzz;
  x12.xy -= i1;
  i = mod(i, 289.0);
  vec3 p = permutation( permutation( i.y + vec3(0.0, i1.y, 1.0 ))+ i.x + vec3(0.0, i1.x, 1.0 ));
  vec3 m = max(0.5 - vec3(dot(x0,x0), dot(x12.xy,x12.xy),dot(x12.zw,x12.zw)), 0.0);
  m = m*m ;
  m = m*m ;
  vec3 x = 2.0 * fract(p * C.www) - 1.0;
  vec3 h = abs(x) - 0.5;
  vec3 ox = floor(x + 0.5);
  vec3 a0 = x - ox;
  m *= 1.7928 - 0.853734 * ( a0*a0 + h*h );
  vec3 g;
  g.x  = a0.x  * x0.x  + h.x  * x0.y;
  g.yz = a0.yz * x12.xz + h.yz * x12.yw;
  return 130.0 * dot(m, g);
}





void main(){


	//vPosition4 = vec4(vPosition.x, vPosition.y, vPosition.z, 1.0);
	float pasoX = _pasoX;
	float pasox = _pasox;
	float pasoY = _pasoY;
	float pasoy = _pasoy;
	float pasoZ = _pasoZ;
	float pasoz = _pasoz;
	vec2 y = vec2( (vPosition.x + pasoX)/8.0,(vPosition.z+pasoZ)/8.0);
	
	//pasoX ++;
	vPosition4 = vec4(vPosition.x + _pasoX, pnoise(y)*2 , vPosition.z + _pasoZ, 1.0);

	gl_Position = u_tProj * u_tView*u_tModel  *vPosition4;

	color = vec4(vColor.r, vColor.g, vColor.b, 1.0);

	texCoord = vTexCoord;
	texCoord = vec2(vTexCoord.x, vTexCoord.y);


	//textCoordsCubo = vPosition;


}