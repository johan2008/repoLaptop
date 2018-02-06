/*--------------
Fragment Shader: filename "fTexture.glsl"
---------------*/
// #version 150   // YJC: Comment/un-comment this line to resolve compilation errors
                  //      due to different settings of the default GLSL version

/**
 * For Mac OS X, GLSL version is by default 120 (OpenGL 2.1 with compatibility profile).
 * We need to use legacy synatx in GLSL 120:
 * 1. Replace 'out' with 'varying' for variables sent from vertex shader.
 * 2. Use built-in variable 'gl_FragColor' for the final output color ('out vec4 fColor' is not needed).
 * 3. Replace 'texture' with 'texture2D' for the texture lookup function.
 */

/* GLSL 120 */
varying vec4 color;
varying vec2 texCoord;

/* GLSL 130+ */
// in vec4 color;
// in vec2 texCoord;
// out vec4 fColor;

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

void main() 
{ 
  /* GLSL 120 */
  if (Texture_app_flag == 0)
     gl_FragColor = color;
  else if (Texture_app_flag == 1)
     gl_FragColor = texture2D( texture_2D, texCoord );
  else // Texture_app_flag == 2
     gl_FragColor = color * texture2D( texture_2D, texCoord );

  /* GLSL 130+ */
//   if (Texture_app_flag == 0)
//      fColor = color;
//   else if (Texture_app_flag == 1)
//      fColor = texture( texture_2D, texCoord );
//   else // Texture_app_flag == 2
//      fColor = color * texture( texture_2D, texCoord );
} 

