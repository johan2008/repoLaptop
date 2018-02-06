#
# Make targets
#
ARCH=$(shell uname -m)
# make source code and copy executables to bin/internal
bin/internal:
	cd code; make
	rsync -rup code/bin/$(ARCH)/ bin/internal
	cd external/SurfCorr2.1; make
	rsync -rup external/SurfCorr2.1/bin/$(ARCH)/ bin/external
	cd bin/external; chmod 755 ./*
	cd scripts; chmod 755 BlendedInterpolation
	cd scripts; chmod 755 TestErrorAllVertices
	cd scripts; chmod 755 IntrinsicSymmetryBlended
	cd scripts; chmod 755 Testme

.PHONY: bin/internal

clean:
	rm -f bin/internal/*
	rm -f bin/external/*
	cd code; make clean
	cd external/SurfCorr2.1; make clean
