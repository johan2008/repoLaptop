

#version 120  // YJC: Comment/un-comment this line to resolve compilation errors
                 //      due to different settings of the default GLSL version

//in  vec4 color;
//out vec4 fColor;


varying vec4 color;
//varying vec4 fColor;


varying  vec3 fN;
//varying  vec3 fL;
varying  vec3 fE;


varying vec4 vPosition4;


uniform vec3 objectColor;
uniform vec3 lightColor;
uniform vec3 lightPos;
uniform vec3 viewPos;


uniform vec3 fL;//LIGHT DIRECTIONAL


uniform vec4 AmbientProduct, DiffuseProduct, SpecularProduct;
uniform float Shininess;
uniform vec4 LightPosition;
uniform int fog_mode;

uniform float flag;

uniform float texture_mapped_ground;



uniform float flagS;


varying vec3 Normal; 

varying vec3 FragPos; 


uniform vec4 CameraEye;
uniform vec4 FogColor;

varying vec4 vertex;

/*    
float getFogFactor(float d)
{
    const float FogMax = 20.0;
    const float FogMin = 10.0;

    if (d>=FogMax) return 1;
    if (d<=FogMin) return 0;

    return 1 - (FogMax - d) / (FogMax - FogMin);
}
*/

varying mat4 u_tModelF;
varying mat4 u_tViewF;


varying  vec2 texCoord;
uniform sampler2D texture_2D;
uniform int Texture_app_flag;

void main() 
{ 


	if(flag>1){

        float density = 0.09, z = length(u_tViewF*u_tModelF*vPosition4), f;


        vec3 Nor = normalize(fN);   
        vec3 E = normalize(fE);  
        vec3 L = normalize(fL); //luz directional

        vec3 H = normalize( L + E );
        
        vec4 ambient = AmbientProduct;

        float Kd = max(dot(L, Nor), 0.0);
        vec4 diffuse = Kd * DiffuseProduct;
        
        float Ks = pow(max(dot(Nor, H), 0.0), Shininess);
        vec4 specular = Ks*SpecularProduct;

        // discard the specular highlight if the light's behind the vertex
        if( dot(L, Nor) < 0.0 ) {
    		specular = vec4(0.0, 0.0, 0.0, 1.0);
        }

        gl_FragColor = (ambient + diffuse + specular);


        if(fog_mode == 0){
            f = 1.0;
        }
        else if(fog_mode == 1){
            f =  (18.0-z)/(18.0-0.0);
        }
        else if(fog_mode == 2){
            f = exp(-density*z);
        }
        else{
            f = exp(-pow(density*z,2));
        }
        f = clamp(f, 0.0, 1.0);
        
        //gl_FragColor = mix(gl_FragColor, FogColor, 1-f);
        

        //gl_FragColor = gl_FragColor* texture2D( texture_2D, texCoord );  



        if (Texture_app_flag == 0)
             gl_FragColor = mix(gl_FragColor, FogColor, 1-f);//color;
        //else if (Texture_app_flag == 1)
        //    gl_FragColor = texture2D( texture_2D, texCoord );
        else if(texture_mapped_ground>0)// Texture_app_flag == 2
            gl_FragColor = mix(gl_FragColor* texture2D( texture_2D, texCoord ), FogColor, 1-f);//     gl_FragColor * texture2D( texture_2D, texCoord ); 
        else
            gl_FragColor = mix(gl_FragColor, FogColor, 1-f);

        /* for light source 2
        L = normalize( pointLightPosition.xyz - pos );
        if(point_light == 0){
            Lf = normalize(pointLightDir.xyz);
        }
        H = normalize( L +  E );
        // get distance attenuation
        dis = length(pos-pointLightPosition.xyz);
        attenuation = 1/(ConstAtt + LinearAtt*dis + QuadAtt*dis*dis);
        if(point_light == 0){   
        //attenuation *= pow(1,pointLightExp);
            if(dot(-L,Lf)>cos(pointLightAng)){
                attenuation *= pow(dot(-L,Lf),pointLightExp);
            }
            else{
                attenuation = 0;
            }
        }

        ambient = attenuation * pointAmbientProduct;
        d = max( dot(L, N), 0.0 );
        diffuse = attenuation * d * pointDiffuseProduct;
        // get specular light
        // pow(x, y) result is undefined if x<0 or if x=0 and y<=0
        s = pow( max(dot(N, H), 0.0), Shininess );
    
        specular = attenuation*s*pointSpecularProduct;
        if( dot(L, N) < 0.0 ) {
            specular = vec4(0.0, 0.0, 0.0, 1.0);
        } 
    */





	}
	else{
	//gl_FragColor = vec4(result,1.0);

        if(texture_mapped_ground>0 ){
            if (Texture_app_flag == 0) gl_FragColor = color;
            else gl_FragColor = color*texture2D( texture_2D, texCoord );
        }else{
            gl_FragColor = color;
        }

	}	
} 

