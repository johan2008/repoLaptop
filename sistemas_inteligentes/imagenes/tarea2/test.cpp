
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <time.h>

#define length(x) (sizeof(x)/sizeof(x[0]))

using namespace std;

int convolucion ()
{
    int mitad, i,j,m,n,mm,nn,ii,jj, acumulador;
    
    // Defino Kernel
    int kernel[3][3] = {
                        {1, 1, 1},
                        {1, 1, 1},
                        {1, 1, 1}
                    };
    // Defino Imagen
    int image[7][7] = { 
                        {9, 2, 7, 8, 1, 6, 5},
                        {5, 5, 6, 7, 8, 2, 5},
                        {5, 5, 6, 7, 8, 2, 5},
                        {5, 5, 6, 7, 8, 2, 5},
                        {1, 4, 2, 4, 3, 8, 1},
                        {8, 4, 9, 1, 3, 7, 3},
                        {3, 6, 0, 4, 5, 8, 7}
                    };     
    // Defino Imagen resultante
    int result[7][7]; 
    // Hallar la mitad del kernel para posicionar la matriz desde ahi
    
    mitad = length(kernel) / 2;
    
    /* Proceso de convolucion
     * Recorro la imagen en los dos primeros for, al igual que el kernel
     * en la variable mm hallo el indice de la fila del kernel alrevez, al
     * igual que la variable nn, almacena la columna del kernel alrevez,
     * las variables ii,jj son para almacenar la posicion de las imagenes tomando
     * en cuenta su limite exterior es decir i-1, j-1,la variable acumulador almacena el resultado
     * que luego es asignado en la posicion de la imagen resultante
     */   

    for (i = 0; i < length(image); ++i) // Filas
    {
        for (j = 0; j < length(image); ++j) // Columnas
        {
            // Variable acumuladora
            acumulador = 0;
            
            for (m = 0; m < length(kernel); ++m) // Filas del Kernel
            {
                mm = length(kernel) - 1 - m; // Indice de la fila del kernel alrevez

                for (n = 0; n < length(kernel); ++n) // Columnas del kernel
                {
                    nn = length(kernel) - 1 - n; // Indice de la columna del kernel alrevez

                    
                    ii = i + (m - mitad);
                    jj = j + (n - mitad);

                    // validar limites de la imagen 00000
                    if (ii >= 0 && ii < length(image) && jj >= 0 && jj < length(image))
                    {
                        acumulador += image[ii][jj] * kernel[mm][nn];
                    }                        
                }
            }
            result[i][j] = acumulador;
        }
    }
    
    // Muestro imagen final aplicando la convolucion
    cout<<endl<<endl<<"Imagen Aplicando convolucion: "<<endl;
    for(int i = 0; i < length(result); i++)
    {
        for(int j = 0; j < length(result); j++)
        {
            cout<<" "<<result[i][j];            
        }
        cout<<endl;
    }
    
    return 0;

}

int main(int argc, char** argv) {

    convolucion();
    
    return 0;
}