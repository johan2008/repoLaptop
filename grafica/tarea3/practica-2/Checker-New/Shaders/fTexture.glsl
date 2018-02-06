/*--------------
Fragment Shader: filename "fTexture.glsl"
---------------*/
 #version 120   // YJC: Comment/un-comment this line to resolve compilation errors
                  //      due to different settings of the default GLSL version

varying  vec4 color;
varying  vec2 texCoord;

uniform sampler2D texture_2D; /* Note: If using multiple textures,
                                       each texture must be bound to a
                                       *different texture unit*, with the
                                       sampler uniform var set accordingly.
                                 The (fragment) shader can access *all texture units*
                                 simultaneously.
                              */
uniform int Texture_app_flag; // 0: no texture application: obj color
                              // 1: texutre color
                              // 2: (obj color) * (texture color)
//varying vec4 fColor;

void main() 
{ 
  if (Texture_app_flag == 0)
     gl_FragColor = color;
  else if (Texture_app_flag == 1)
     gl_FragColor = texture2D( texture_2D, texCoord );
  else // Texture_app_flag == 2
     gl_FragColor = color * texture2D( texture_2D, texCoord );  
} 

