

#version 120  // YJC: Comment/un-comment this line to resolve compilation errors
                 //      due to different settings of the default GLSL version

varying vec4 color;

void main() 
{ 

	gl_FragColor = color;
    //gl_FragColor = vec4(1.0f, 0.5f, 0.2f, 1.0f);
} 

