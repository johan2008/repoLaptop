
#version 120


varying vec4 color;
varying  vec2 texCoord;
varying vec4 vPosition4;


uniform sampler2D texture_2D;
uniform sampler2D texture2;
uniform sampler2D texture3;
uniform sampler2D texture4;
uniform sampler2D texture5;
uniform sampler2D texture6;

//varying vec3 textureDir;
//uniform samplerCube cubemap;


varying vec3 textCoordsCubo;
uniform samplerCube cubemap;




void main(){



	//gl_FragColor = color;

	//gl_FragColor = texture2D( texture_2D, texCoord );

	//gl_FragColor = mix(texture2D(texture_2D,texCoord), texture2D(texture2,texCoord),0.2);

	//vec4 tex0, tex1, tex2;
	//tex0 = texture2D(texture_2D,texCoord);
	//tex1 = texture2D(texture1,texCoord);



	//gl_FragColor = mix( texture2D(texture2,texCoord) ,texture2D(texture_2D,texCoord), 0.1);




	if(vPosition4.y<0.05){
		gl_FragColor = texture2D( texture_2D, texCoord );
	}
	else if( 0.05<=vPosition4.y && vPosition4.y<0.3){
		gl_FragColor = texture2D( texture2, texCoord );
	}else if( 0.3<=vPosition4.y && vPosition4.y<0.35){
		gl_FragColor = texture2D( texture3, texCoord );
	}else if( 0.35<=vPosition4.y && vPosition4.y<0.4){
		gl_FragColor = texture2D( texture4, texCoord );
	}else if( 0.4 <=vPosition4.y && vPosition4.y<1.4){
		gl_FragColor = texture2D( texture5, texCoord );
	}
	else{
		gl_FragColor = texture2D( texture6, texCoord );
	}


	//gl_FragColor = texture2D( texture3, texCoord );
	//gl_FragColor = textureCube(cubemap, textCoordsCubo);
}
		