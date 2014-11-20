RivetMC_BOOSTEDHBB.so: MC_BOOSTEDHBB.cc MC_BOOSTEDHBB.hh
	rivet-buildplugin RivetMC_BOOSTEDHBB.so MC_BOOSTEDHBB.cc `fastjet-config --prefix`/lib/libVariableR.a
