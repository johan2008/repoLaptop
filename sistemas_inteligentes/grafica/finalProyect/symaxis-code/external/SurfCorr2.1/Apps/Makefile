#
# Makefile for GAPS
#

opt:
	$(MAKE) target "TARGET=$@"

debug:
	$(MAKE) target "TARGET=$@"

clean:
	$(MAKE) target "TARGET=$@"

release:
	mkdir -p release
	cp README.txt release
	cp LICENSE.txt release
	cp Makefile release
	cp -r scripts release
	$(MAKE) target "TARGET=$@"

target: 
	cd SurfaceVis; $(MAKE) $(TARGET)





