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
	cd BasicRoutines; $(MAKE) $(TARGET)
	cd BasicNumerical; $(MAKE) $(TARGET)
	cd SampledSurface; $(MAKE) $(TARGET)
	cd GLUTSurfaceVisLib; $(MAKE) $(TARGET)







