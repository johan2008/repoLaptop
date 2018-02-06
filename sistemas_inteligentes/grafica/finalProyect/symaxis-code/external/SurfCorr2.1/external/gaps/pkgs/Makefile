#
# Makefile for GAPS pkgs
#

opt:
	$(MAKE) target "TARGET=$@"

debug:
	$(MAKE) target "TARGET=$@"

clean:
	$(MAKE) target "TARGET=$@"

release:
	mkdir -p ../release/pkgs
	cp Makefile ../release/pkgs
	$(MAKE) target "TARGET=$@"

target: 
	cd RNBasics; $(MAKE) $(TARGET)
	cd R2Shapes; $(MAKE) $(TARGET)
	cd R3Shapes; $(MAKE) $(TARGET)
	cd R3Graphics; $(MAKE) $(TARGET)
	cd jpeg; $(MAKE) $(TARGET)





