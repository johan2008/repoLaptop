OLD README FILE - NOT RELEVANT


INSTALLAON:
   Step 1 GSL:
	cd external/gsl-1.14/
	./configure
	make

   Step 2 BasicNumerical:
	qmake
	make

	or 
	make

This is basic library for numerical operations. 

TODO: include external libraries in this distribution!!!

It basically serves as a wrapper for GSL for dense matrices.
Currently Implemented:
	Comlplex / Real Matrices + Vectors
	Matrix (real): Inverse, LU, Cholesky, SVD
	FFT


INSTALLATION:
	gsl - http://www.gnu.org/software/gsl/	

	NOT REQUIRED: qt - http://qt.nokia.com/downloads
	
	NOT IMPLEMENTED: glpk - http://www.gnu.org/software/glpk/


Linear Programming:
	http://www.mathtools.net/C_C__/Optimization/


Sparse matrices are to come.
	http://www.cise.ufl.edu/research/sparse/
	Cholmode:
		http://www.cise.ufl.edu/research/sparse/cholmod/
	Matlab uses:	http://www.mathworks.com/access/
			/helpdesk/help/techdoc/ref/mldivide.html

	Other packages:
		http://www.mathtools.net/C_C__/Sparse/index.html

	SparceLib + IML

