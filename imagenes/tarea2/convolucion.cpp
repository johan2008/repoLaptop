#include <fstream>
#include <iostream>
#include <ctime>
#include <string>
#include <array>
#include <vector>
#include <iterator>
#include <cmath>
#include <algorithm>    // std::sort


#include <omp.h>

using namespace std;


#define length(x) (sizeof(x)/sizeof(x[0]))



void showData(const char *fileName){
    unsigned char info[54];
    FILE *nameFile =  fopen(fileName, "rb");
    fread(info, sizeof(unsigned char), 54 , nameFile); //read 54 bytes header
    // extract image height and wight from header
    int width = *(int*)&info[18];
    int height = *(int*)&info[22];


    cout<<"show data"<<endl;
    cout<<"width:  "<< width<<endl;
    cout<<"height:  "<<height<<endl;

    int sizeImage = 3*width*height;
    cout<<"size of image: "<<sizeImage<<endl;

    int row_padded = (width*3 +3 ) & (~3);
    cout<<"row_padded:  "<< row_padded<<endl;



}




unsigned char* readBMP(char* filename)
{
    int i;
    FILE* f = fopen(filename, "rb");
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int width = *(int*)&info[18];
    int height = *(int*)&info[22];

    cout<<"width:  "<<width<<endl;
    cout<<"height:   "<<height<<endl;  


    int size = 3 * width * height;
    unsigned char* data = new unsigned char[size]; // allocate 3 bytes per pixel
    fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
    fclose(f);

    for(i = 0; i < size; i += 3)
    {
            unsigned char tmp = data[i];
            data[i] = data[i+2];
            data[i+2] = tmp;
    }


    int filesize2 = width*height*3;

    unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
    unsigned char bmppad[3] = {0,0,0};

    bmpfileheader[ 2] = (unsigned char)(filesize2    );
    bmpfileheader[ 3] = (unsigned char)(filesize2>> 8);
    bmpfileheader[ 4] = (unsigned char)(filesize2>>16);
    bmpfileheader[ 5] = (unsigned char)(filesize2>>24);

    bmpinfoheader[ 4] = (unsigned char)(       width    );
    bmpinfoheader[ 5] = (unsigned char)(       width>> 8);
    bmpinfoheader[ 6] = (unsigned char)(       width>>16);
    bmpinfoheader[ 7] = (unsigned char)(       width>>24);
    bmpinfoheader[ 8] = (unsigned char)(       height    );
    bmpinfoheader[ 9] = (unsigned char)(       height>> 8);
    bmpinfoheader[10] = (unsigned char)(       height>>16);
    bmpinfoheader[11] = (unsigned char)(       height>>24);


    //rotation
    int angle=15;        //45° for example 
    //Convert degrees to radians 
    float radians=(2*3.1416*angle)/360; 

    float cosine=(float)cos(radians); 
    float sine=(float)sin(radians); 
    unsigned char* dataRotate = new unsigned char [width * height * 3];

    int x, y;

    int valor;

    cout<<"w: "<<width<<endl;
    cout<<"h: "<<height<<endl;
    for(int cy = 0; cy < height; cy++)
    {
        for(int cx = 0; cx < width; cx++)
        {
            for (int c = 0; c < 3; c++)
            {
                x = round (cos(radians)*cx +sin(radians)*cy );               
                y = round (-sin(radians)*cx + cos(radians)*cy );

                if(x<0) x = 0;
                if(y<0) y = 0;

                cout<<"x:  " << x <<endl;
                cout<<"y:  " << y <<endl;

                valor = data[3*(x+y*width) + c];

                dataRotate[3*(cx+cy*width) + c] = valor;

            }
            
        }
    }



    FILE* f2; 
    f2 = fopen("imgRotate.bmp","wb");
    fwrite(bmpfileheader,1,14,f2);
    fwrite(bmpinfoheader,1,40,f2);
    /*for(int i=0; i<newHeight; i++)
    {
        fwrite(newData+(newWidth*(newHeight-i-1)*3),3,newWidth,f);
        fwrite(bmppad,1,(4-(newWidth*3)%4)%4,f);
    }*/

    for(int i=height -1; i>=0; i--)
    {
        fwrite(dataRotate+(width*(height-i-1)*3),3,width,f2);
        fwrite(bmppad,1,(4-(width*3)%4)%4,f2);
    }






    return data;
}




void rotateImage(char* filename){
    FILE* f = fopen(filename, "rb");
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int width = *(int*)&info[18];
    int height = *(int*)&info[22];

    int size = 3 * width * height;
    unsigned char* data = new unsigned char[size]; // allocate 3 bytes per pixel
    fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
    fclose(f);

    /*for(int i = 0; i < size; i += 3)
    {
            unsigned char tmp = data[i];
            data[i] = data[i+2];
            data[i+2] = tmp;
    }*/


    int filesize2 = width*height*3;

    unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
    unsigned char bmppad[3] = {0,0,0};

    bmpfileheader[ 2] = (unsigned char)(filesize2    );
    bmpfileheader[ 3] = (unsigned char)(filesize2>> 8);
    bmpfileheader[ 4] = (unsigned char)(filesize2>>16);
    bmpfileheader[ 5] = (unsigned char)(filesize2>>24);

    bmpinfoheader[ 4] = (unsigned char)(       width    );
    bmpinfoheader[ 5] = (unsigned char)(       width>> 8);
    bmpinfoheader[ 6] = (unsigned char)(       width>>16);
    bmpinfoheader[ 7] = (unsigned char)(       width>>24);
    bmpinfoheader[ 8] = (unsigned char)(       height    );
    bmpinfoheader[ 9] = (unsigned char)(       height>> 8);
    bmpinfoheader[10] = (unsigned char)(       height>>16);
    bmpinfoheader[11] = (unsigned char)(       height>>24);


    int angle=15;        //45° for example 
    //Convert degrees to radians 
    float radians=(2*3.1416*angle)/360; 

    float cosine=(float)cos(radians); 
    float sine=(float)sin(radians); 
    unsigned char* dataRotate = new unsigned char [width * height * 3];
    int x, y;
    int valor;


    double t0,t1;
    //t0 = omp_get_wtime();

    int cy, cx,c;
    //#pragma omp parallel for schedule(static) private(cy,cx) //reduction(+:suma)
    for(cy = 0; cy < height; cy++)
    {
        for(int cx = 0; cx < width; cx++)
        {
            //#pragma omp critical
            for ( c = 0; c < 3; c++)
            {
                x = round (cos(radians)*cx +sin(radians)*cy );               
                y = round (-sin(radians)*cx + cos(radians)*cy );

                if(x<0) x = 0;
                if(y<0) y = 0;

                valor = data[3*(x+y*width) + c];
                dataRotate[3*(cx+cy*width) + c] = valor;

            }



            
        }
    }
    //t1 = omp_get_wtime();
    cout<<" time in rotate image:  "<<(t1-t0)<<endl;


    FILE* f2; 
    f2 = fopen("imgRotate.bmp","wb");
    fwrite(bmpfileheader,1,14,f2);
    fwrite(bmpinfoheader,1,40,f2);
    /*for(int i=0; i<newHeight; i++)
    {
        fwrite(newData+(newWidth*(newHeight-i-1)*3),3,newWidth,f);
        fwrite(bmppad,1,(4-(newWidth*3)%4)%4,f);
    }*/

    for(int i=height -1; i>=0; i--)
    {
        fwrite(dataRotate+(width*(height-i-1)*3),3,width,f2);
        fwrite(bmppad,1,(4-(width*3)%4)%4,f2);
    }    


}



void imageScale(int scaleX, int scaleY, const char *fileName){
    int i;
    FILE* f = fopen(fileName, "rb");
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int width = *(int*)&info[18];
    int height = *(int*)&info[22];

    int size = 3 * width * height;
    unsigned char* data = new unsigned char[size]; // allocate 3 bytes per pixel
    fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
    fclose(f);

    int newWidth = width*scaleX;
    int newHeight = height*scaleY;
    unsigned char* newData = new unsigned char [newWidth * newHeight * 3];
    double scaleWidth =  (double)newWidth / (double)width;
    double scaleHeight = (double)newHeight / (double)height;



    double t0,t1;
    //t0 = omp_get_wtime();

    int cy, cx;
    //#pragma omp parallel for schedule(static) private(cy,cx) //reduction(+:suma)
    for( cy = 0; cy < newHeight; cy++)
    {
        for( cx = 0; cx < newWidth; cx++)
        {
            int pixel = (cy * (newWidth *3)) + (cx*3);
            int nearestMatch =  (    ((int)(cy / scaleHeight) * (width *3)) +   ((int)(cx / scaleWidth) *3)    );
                
            newData[pixel    ] =  data[nearestMatch    ];
            newData[pixel + 1] =  data[nearestMatch + 1];
            newData[pixel + 2] =  data[nearestMatch + 2];
        }
    }


    //t1 = omp_get_wtime();
    cout<<" time in scale image:  "<<(t1-t0)<<endl;



    int filesize = newWidth*newHeight*3;

    unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
    unsigned char bmppad[3] = {0,0,0};

    bmpfileheader[ 2] = (unsigned char)(filesize    );
    bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
    bmpfileheader[ 4] = (unsigned char)(filesize>>16);
    bmpfileheader[ 5] = (unsigned char)(filesize>>24);

    bmpinfoheader[ 4] = (unsigned char)(       newWidth    );
    bmpinfoheader[ 5] = (unsigned char)(       newWidth>> 8);
    bmpinfoheader[ 6] = (unsigned char)(       newWidth>>16);
    bmpinfoheader[ 7] = (unsigned char)(       newWidth>>24);
    bmpinfoheader[ 8] = (unsigned char)(       newHeight    );
    bmpinfoheader[ 9] = (unsigned char)(       newHeight>> 8);
    bmpinfoheader[10] = (unsigned char)(       newHeight>>16);
    bmpinfoheader[11] = (unsigned char)(       newHeight>>24);
    f = fopen("img.bmp","wb");
    fwrite(bmpfileheader,1,14,f);
    fwrite(bmpinfoheader,1,40,f);

    for(int i=newHeight -1; i>=0; i--)
    {
    	//cout<<"scale:"<<endl;
        fwrite(newData+(newWidth*(newHeight-i-1)*3),3,newWidth,f);
        fwrite(bmppad,1,(4-(newWidth*3)%4)%4,f);
    }




}

int  getMin(int  array[], int  size){
	int min = 0;
	for(int i = 0; i < size; i++)
		if(array[i] < min)
			min = array[i];
	return  min;
}

int convolucion(const char *fileName){
	FILE* f = fopen(fileName, "rb");
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int width = *(int*)&info[18];
    int height = *(int*)&info[22];

    int size = 3 * width * height;
    unsigned char* data = new unsigned char[size]; // allocate 3 bytes per pixel
    fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
    fclose(f);





    int kernel[3][3] = {
                        {0, 0, 0},
                        {0, 1, 0},
                        {0, 0, 0}
                        };

    int filesize2 = width*height*3;

    unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
    unsigned char bmppad[3] = {0,0,0};

    bmpfileheader[ 2] = (unsigned char)(filesize2    );
    bmpfileheader[ 3] = (unsigned char)(filesize2>> 8);
    bmpfileheader[ 4] = (unsigned char)(filesize2>>16);
    bmpfileheader[ 5] = (unsigned char)(filesize2>>24);

    bmpinfoheader[ 4] = (unsigned char)(       width    );
    bmpinfoheader[ 5] = (unsigned char)(       width>> 8);
    bmpinfoheader[ 6] = (unsigned char)(       width>>16);
    bmpinfoheader[ 7] = (unsigned char)(       width>>24);
    bmpinfoheader[ 8] = (unsigned char)(       height    );
    bmpinfoheader[ 9] = (unsigned char)(       height>> 8);
    bmpinfoheader[10] = (unsigned char)(       height>>16);
    bmpinfoheader[11] = (unsigned char)(       height>>24);


    int mitad, i,j,m,n,mm,nn,ii,jj, acumulador;


 	mitad = length(kernel) / 2;
 	//int acumulador = 0;



 	unsigned char* dataRotate = new unsigned char [width * height * 3];


 	//for (i = 0; i < length(image); ++i) // 

 	//cout<<"height convolucion:   "<<height<<endl;
 	//cout<<"width convolucion:    "<<width<<endl;

 	double t0,t1;
 	t0 = omp_get_wtime();

 	#pragma omp parallel for schedule(static) private(i,j) //reduction(+:suma)

 	for ( i = 0; i < height; ++i) // Filas
    
    {
        //for (j = 0; j < length(image); ++j) // Columnas
        for ( j = 0; j < width; ++j) // Columnas
        {
            // Variable acumuladora
            acumulador = 0;
            

            #pragma omp critical
            for ( int c = 0; c < 3; c++){


            for (m = 0; m < length(kernel); ++m) // Filas del Kernel
            {
                mm = length(kernel) - 1 - m; // Indice de la fila del kernel alrevez

                for (n = 0; n < length(kernel); ++n) // Columnas del kernel
                {
                    nn = length(kernel) - 1 - n; // Indice de la columna del kernel alrevez

                    
                    ii = i + (m - mitad);
                    jj = j + (n - mitad);

                    // validar limites de la imagen 00000
                    if (ii >= 0 && ii < height && jj >= 0 && jj < width)
                    {
                        //acumulador += data[ii][jj] * kernel[mm][nn];
                    	acumulador += data[3*(jj+ii*width) + c] * kernel[mm][nn];

                    }                        
                }
            }
            //dataRotate[i][j] = acumulador;


            //valor = data[3*(x+y*width) + c];

            dataRotate[3*(j+i*width) + c] = acumulador;

            }




        }
	}

	t1 = omp_get_wtime();
    cout<<" time in convolucion image:  "<<(t1-t0)<<endl;

	FILE* f2; 
    f2 = fopen("convolucion.bmp","wb");
    fwrite(bmpfileheader,1,14,f2);
    fwrite(bmpinfoheader,1,40,f2);
    /*for(int i=0; i<newHeight; i++)
    {
        fwrite(newData+(newWidth*(newHeight-i-1)*3),3,newWidth,f);
        fwrite(bmppad,1,(4-(newWidth*3)%4)%4,f);
    }*/

    for(int i=height -1; i>=0; i--)
    {
        fwrite(dataRotate+(width*(height-i-1)*3),3,width,f2);
        fwrite(bmppad,1,(4-(width*3)%4)%4,f2);
    }  



}


void ecualizacion(const char *fileName){
	FILE* f = fopen(fileName, "rb");
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int width = *(int*)&info[18];
    int height = *(int*)&info[22];

    int size = 3 * width * height;
    unsigned char* data = new unsigned char[size]; // allocate 3 bytes per pixel
    fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
    fclose(f);





    int kernel[3][3] = {
                        {0, 0, 0},
                        {0, 1, 0},
                        {0, 0, 0}
                        };

    int filesize2 = width*height*3;

    unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
    unsigned char bmppad[3] = {0,0,0};

    bmpfileheader[ 2] = (unsigned char)(filesize2    );
    bmpfileheader[ 3] = (unsigned char)(filesize2>> 8);
    bmpfileheader[ 4] = (unsigned char)(filesize2>>16);
    bmpfileheader[ 5] = (unsigned char)(filesize2>>24);

    bmpinfoheader[ 4] = (unsigned char)(       width    );
    bmpinfoheader[ 5] = (unsigned char)(       width>> 8);
    bmpinfoheader[ 6] = (unsigned char)(       width>>16);
    bmpinfoheader[ 7] = (unsigned char)(       width>>24);
    bmpinfoheader[ 8] = (unsigned char)(       height    );
    bmpinfoheader[ 9] = (unsigned char)(       height>> 8);
    bmpinfoheader[10] = (unsigned char)(       height>>16);
    bmpinfoheader[11] = (unsigned char)(       height>>24);

 	unsigned char* dataRotate = new unsigned char [width * height * 3];



    int hist[256];
    for(int i=0; i<256; i++) hist[i] = 0;

    for(int y=0; y<height; y++){
    	for (int x = 0; x < width; ++x)
    	{
    		/*for(int c=0; c<3; c++){
    			int ind = data[3*(x+y*width) + c];	
    			hist[ind] = hist[ind] +1;
    		}*/
    		int indR = data[3*(x+y*width) + 0];	
    		hist[indR] = hist[indR] +1;
    		int indG = data[3*(x+y*width) + 1];	
    		hist[indG] = hist[indG] +1;
    		int indB = data[3*(x+y*width) + 2];	
    		hist[indB] = hist[indB] +1;
    	}
    }

    int accum_hist[256];
    accum_hist[0] = hist[0];
    for (int i = 1; i < 256; ++i) accum_hist[i] = hist[i] + accum_hist[i-1];

    /*FILE* pFile;
	pFile = fopen("hist.txt","w");
	fprintf(pFile ,"pix ,Hist ,acHist ,ImSize: %dx %d\n",height , width);
	for (int i = 0; i < 256; ++i)
	{
		fprintf(pFile ," %d, %d, %d\n",i,hist[i],accum_hist[i]);
	}
	fclose(pFile);
	*/

	int accumMinimum = getMin(accum_hist, 256);
	int lookUp[256];
	for (int i = 0; i < 256; ++i)
	{
		lookUp[i] = floor (255*( accum_hist[i]-accumMinimum)/(height*width -accumMinimum));
	}

	int y, x;
	double t0,t1;
 	t0 = omp_get_wtime();

 	#pragma omp parallel for schedule(static) private(y,x) 
	for(y=0; y<height; y++){
    	for (x = 0; x < width; ++x)
    	{
    		#pragma omp critical
    		for(int c=0; c<3; c++){
    			//int ind = data[3*(x+y*width) + c];	
    			//hist[ind] = hist[ind] +1;
    			dataRotate[3*(x+y*width) + c] = lookUp[  data[3*(x+y*width) +c]  ];
    		}
    	}
    }

    t1 = omp_get_wtime();
    cout<<" time in ecualizacion image:  "<<(t1-t0)<<endl;

    FILE* f2; 
    f2 = fopen("hist.bmp","wb");
    fwrite(bmpfileheader,1,14,f2);
    fwrite(bmpinfoheader,1,40,f2);
    /*for(int i=0; i<newHeight; i++)
    {
        fwrite(newData+(newWidth*(newHeight-i-1)*3),3,newWidth,f);
        fwrite(bmppad,1,(4-(newWidth*3)%4)%4,f);
    }*/

    for(int i=height -1; i>=0; i--)
    {
        fwrite(dataRotate+(width*(height-i-1)*3),3,width,f2);
        fwrite(bmppad,1,(4-(width*3)%4)%4,f2);
    }  

}


int medianImage(const char *fileName){
    FILE* f = fopen(fileName, "rb");
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int width = *(int*)&info[18];
    int height = *(int*)&info[22];

    int size = 3 * width * height;
    unsigned char* data = new unsigned char[size]; // allocate 3 bytes per pixel
    fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
    fclose(f);





    int kernel[3][3] = {
                        {0, 0, 0},
                        {0, 1, 0},
                        {0, 0, 0}
                        };

    int filesize2 = width*height*3;

    unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
    unsigned char bmppad[3] = {0,0,0};

    bmpfileheader[ 2] = (unsigned char)(filesize2    );
    bmpfileheader[ 3] = (unsigned char)(filesize2>> 8);
    bmpfileheader[ 4] = (unsigned char)(filesize2>>16);
    bmpfileheader[ 5] = (unsigned char)(filesize2>>24);

    bmpinfoheader[ 4] = (unsigned char)(       width    );
    bmpinfoheader[ 5] = (unsigned char)(       width>> 8);
    bmpinfoheader[ 6] = (unsigned char)(       width>>16);
    bmpinfoheader[ 7] = (unsigned char)(       width>>24);
    bmpinfoheader[ 8] = (unsigned char)(       height    );
    bmpinfoheader[ 9] = (unsigned char)(       height>> 8);
    bmpinfoheader[10] = (unsigned char)(       height>>16);
    bmpinfoheader[11] = (unsigned char)(       height>>24);


    int mitad, i,j,m,n,mm,nn,ii,jj, acumulador;


    mitad = length(kernel) / 2;
    //int acumulador = 0;



    unsigned char* dataRotate = new unsigned char [width * height * 3];


    //for (i = 0; i < length(image); ++i) // 

    //cout<<"height convolucion:   "<<height<<endl;
    //cout<<"width convolucion:    "<<width<<endl;

    double t0,t1;
    t0 = omp_get_wtime();

    std::vector<int> aux;
    int median;

    #pragma omp parallel for schedule(static) private(i,j) //reduction(+:suma)

    for ( i = 0; i < height; ++i) // Filas
    
    {
        //for (j = 0; j < length(image); ++j) // Columnas
        for ( j = 0; j < width; ++j) // Columnas
        {
            // Variable acumuladora
            acumulador = 0;
            

            #pragma omp critical
            for ( int c = 0; c < 3; c++){


            for (m = 0; m < length(kernel); ++m) // Filas del Kernel
            {
                mm = length(kernel) - 1 - m; // Indice de la fila del kernel alrevez

                for (n = 0; n < length(kernel); ++n) // Columnas del kernel
                {
                    nn = length(kernel) - 1 - n; // Indice de la columna del kernel alrevez

                    
                    ii = i + (m - mitad);
                    jj = j + (n - mitad);

                    // validar limites de la imagen 00000
                    if (ii >= 0 && ii < height && jj >= 0 && jj < width  )
                    {
                        //acumulador += data[3*(jj+ii*width) + c] * kernel[mm][nn];
                        aux.push_back(data[3*(jj+ii*width) + c]);    //+= data[3*(jj+ii*width) + c] * kernel[mm][nn];

                        //cout<<"a   ";
                    }                        
                }
            }
            //cout<<"---------------------------------------------------------"<<endl;
            std::sort (aux.begin(), aux.begin() +aux.size());


            if(aux.size()%2 == 0){
                median = (aux[aux.size()/2] +  aux[aux.size()/2 -1])/2;     
            }else{
                median = aux[aux.size()/2 ]; 
            }


            //dataRotate[i][j] = acumulador;


            //valor = data[3*(x+y*width) + c];

            dataRotate[3*(j+i*width) + c] = acumulador;

            }




        }
    }

    t1 = omp_get_wtime();
    cout<<" time in convolucion image:  "<<(t1-t0)<<endl;

    FILE* f2; 
    f2 = fopen("mediana.bmp","wb");
    fwrite(bmpfileheader,1,14,f2);
    fwrite(bmpinfoheader,1,40,f2);
    /*for(int i=0; i<newHeight; i++)
    {
        fwrite(newData+(newWidth*(newHeight-i-1)*3),3,newWidth,f);
        fwrite(bmppad,1,(4-(newWidth*3)%4)%4,f);
    }*/

    for(int i=height -1; i>=0; i--)
    {
        fwrite(dataRotate+(width*(height-i-1)*3),3,width,f2);
        fwrite(bmppad,1,(4-(width*3)%4)%4,f2);
    }  



}




int main(){
    showData("bee.bmp");
    //readBMP("lena.bmp");
    

    imageScale(3,3,"bee.bmp");


    rotateImage("bee.bmp");

    convolucion("bee.bmp");

    ecualizacion("bee.bmp");

    medianImage("bee.bmp");

    return 0;
}