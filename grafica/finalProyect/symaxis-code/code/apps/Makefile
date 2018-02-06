#
# Makefile for GAPS apps
#

opt:
	$(MAKE) target "TARGET=$@"

debug:
	$(MAKE) target "TARGET=$@"

clean:
	$(MAKE) target "TARGET=$@"

release:
	mkdir -p ../release/apps
	cp Makefile ../release/apps
	$(MAKE) target "TARGET=$@"

target: 
	cd anchors_to_canonic_align; $(MAKE) $(TARGET)
	cd axis2align2; $(MAKE) $(TARGET)
	cd axis2dist; $(MAKE) $(TARGET)
	cd axis2prp; $(MAKE) $(TARGET)
	cd axisalign_symtips_to_tipscor; $(MAKE) $(TARGET)
	cd axislength; $(MAKE) $(TARGET)
	cd axisrefine; $(MAKE) $(TARGET)
	cd cor2triplets; $(MAKE) $(TARGET)
	cd map2error; $(MAKE) $(TARGET)
	cd msh2msh; $(MAKE) $(TARGET)
	cd myprp2prp; $(MAKE) $(TARGET)
	cd prp_norm; $(MAKE) $(TARGET)
	cd prp2persistence; $(MAKE) $(TARGET)
	cd refinecoarse; $(MAKE) $(TARGET)
	cd symaxis2tipscor; $(MAKE) $(TARGET)
	cd symmap2tipscor; $(MAKE) $(TARGET)
	cd msh2prp; $(MAKE) $(TARGET)
	cd check_axis_consistency; $(MAKE) $(TARGET)
	cd symdist2axis; $(MAKE) $(TARGET)
	cd symdist2binary; $(MAKE) $(TARGET)
	cd symmap2dist; $(MAKE) $(TARGET)
	cd precisemap2vertexmap; $(MAKE) $(TARGET)
	cd intrinsic_triplets_selection; $(MAKE) $(TARGET)
	cd mshalign; $(MAKE) $(TARGET)
	cd pca_extrinsic_sym_fast; $(MAKE) $(TARGET)
	cd vprp2prp; $(MAKE) $(TARGET)
	cd vpts2pts; $(MAKE) $(TARGET)
	cd symmaprefine; $(MAKE) $(TARGET)
