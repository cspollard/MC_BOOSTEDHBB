RivetMC_BOOSTEDHBB.so:
#	rivet-buildplugin `/afs/phas.gla.ac.uk/user/a/amorton/rivet/fastjet-3.0.6/fastjet-config --cxxflags --libs --plugins`  MC_BOOSTEDHBB.cc
	rivet-buildplugin RivetMC_BOOSTEDHBB.so MC_BOOSTEDHBB.cc `fastjet-config --prefix`/lib/libVariableR.a 
