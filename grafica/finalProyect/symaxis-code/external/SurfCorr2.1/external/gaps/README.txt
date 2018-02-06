GAPS Users -

This directory contains all code for the GAPS software library.
There are several subdirectories:

    pkgs - source and include files for all packages (software libraries).
    apps - source files for several application and example programs. 
    makefiles - unix-style make file definitions
    vc - visual studio solution files
    lib - archive library (.lib) files.
    bin - executable files.

If you are using linux or cygwin and have gcc, OpenGL, and glut installed, 
you should be able to compile the code by typing "make clean; make"
in this directory.  The makefiles directory has the shared makefile
settings which could be edited for other operating systems and compilers.
To write a linux program that simply uses GAPS, then you should include
"-I XXX/gaps/pkgs" in your compile flags (CFLAGS) and "-L XXX/gaps/lib" 
in your link flags (LDFLAGS), where XXX is the directory where you 
installed the gaps software.  

If you are using Windows Visual Studio, then you should be able to 
open the solution file vc.sln in the vc subdirectory and then 
"Rebuild Solution."

The software is distributed under the MIT license (see LICENSE.txt)
and thus can be used for any purpose without warranty of any kind.

- Tom Funkhouser



