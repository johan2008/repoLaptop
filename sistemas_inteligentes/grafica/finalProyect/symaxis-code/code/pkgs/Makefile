#
# Makefile for GAPS pkgs
#


# 
# Library file names
#

ARCH=$(shell uname -m)
LIB_DIR=../lib/$(ARCH)
LIB=$(LIB_DIR)/libgaps.a
PKG_LIBS= \
  $(LIB_DIR)/libR3Graphics.a \
  $(LIB_DIR)/libR3Shapes.a \
  $(LIB_DIR)/libR2Shapes.a \
  $(LIB_DIR)/libRNBasics.a \
  $(LIB_DIR)/libjpeg.a \
  $(LIB_DIR)/libpng.a \
  $(LIB_DIR)/libRNExt.a \
  $(LIB_DIR)/libBasicDataStructure.a \
  $(LIB_DIR)/libBasicNumerical.a \
  $(LIB_DIR)/libANN.a


# 
# Make targets
#

opt:
	$(MAKE) target "TARGET=$@"
	rm -f $(LIB)
	ar ur $(LIB) $(PKG_LIBS)

debug:
	$(MAKE) target "TARGET=$@"
	rm -f $(LIB)
	ar ur $(LIB) $(PKG_LIBS)

clean:
	$(MAKE) target "TARGET=$@"
	rm -f $(LIB)

release:
	mkdir -p ../release/pkgs
	cp Makefile ../release/pkgs
	$(MAKE) target "TARGET=$@"

target: 
	cd RNBasics; $(MAKE) $(TARGET)
	cd R2Shapes; $(MAKE) $(TARGET)
	cd R3Shapes; $(MAKE) $(TARGET)
	cd R3Graphics; $(MAKE) $(TARGET)
	cd fglut; $(MAKE) $(TARGET)
	cd jpeg; $(MAKE) $(TARGET)
	cd png; $(MAKE) $(TARGET)
	cd RNExt; $(MAKE) $(TARGET)
	cd BasicDataStructure; $(MAKE) $(TARGET)
	cd BasicNumerical; $(MAKE) $(TARGET)
	cd ANN; $(MAKE) $(TARGET)
