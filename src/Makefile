include $(PTOOLS)/config/Makefile.common

MODULES  = PUnits PConst PFunctions PData PDataHiP PPalette

LIB      = ptools

HDF5CFLAGS    = -I$(HDF5CPP)/include
HDF5LIBS      = -L$(HDF5CPP)/lib -L$(SZIP)/lib -lhdf5_cpp -lhdf5 -lz

# c++-compiler
CXXFLAGS      += $(HDF5CFLAGS) -I$(INCDIR)
LIBS          := $(ROOTLIBS) $(HDF5LIBS) 

SHLIBFILE    := lib$(LIB).$(DllSuf)

RLIBDEP	      = libCore.so libCint.so libRIO.so libHist.so libTree.so libPhysics.so \
		libMathCore.so libhdf5_cpp.so

ALLDEPEND     = $(addsuffix .d,$(MODULES))

all: depend lib libmap 

depend: $(ALLDEPEND)

modules: $(addsuffix .o,$(MODULES))

lib: $(LIBDIR)/$(SHLIBFILE)

info:	
	@echo "CURDIR   = $(CURDIR)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "Library  : $(LIBDIR)/$(SHLIBFILE)"
	@echo "SOFLAGS  = $(SOFLAGS)"
	@echo "LDFLAGS  = $(LDFLAGS)"
	@echo "ROOTLIBS = $(ROOTLIBS)"


$(LIBDIR)/$(SHLIBFILE): $(addsuffix .o,$(MODULES)) $(LIB)Dict.o 
	@echo "Building shared library $@ ... "
	@mkdir -p $(LIBDIR)
ifeq ($(PLATFORM),macosx)
	$(LD) $(SOFLAGS) $@ $^ -o $@ $(LIBS)
else
	$(LD) $(SOFLAGS) $^ -o $@ $(LIBS)
endif
#	ln -fs $(SHLIBFILE) lib/lib$(LIB).1.so
	@echo "Done"

libmap: $(LIBDIR)/lib$(LIB).rootmap 

$(LIBDIR)/lib$(LIB).rootmap: $(LIBDIR)/lib$(LIB).so 
	@echo "Creating library rootmap $@"
	rlibmap -f -o $@ -d $(SHLIBFILE) $(RLIBDEP) -c $(INCDIR)/LinkDef.hh

$(LIB)Dict.cc: 
	@echo "Generating dictionary $@"
	rootcint -f $(LIB)Dict.cc -c -m64 -p $(addprefix $(INCDIR)/,$(addsuffix .hh,$(MODULES))) $(INCDIR)/LinkDef.hh 


%.d: %.c
	@echo Making dependency for file $< ...
	@set -e;\
	$(CXX) -M -DMKDEPEND $(CXXFLAGS) $< |\
	sed 's!$*\.o!& $@!' >$@;\
	[ -s $@ ] || rm -f $@

%.d: %.cc
	@echo Making dependency for file $< ...
	@set -e;\
	$(CXX) -M -DMKDEPEND $(CXXFLAGS) $< |\
	sed 's!$*\.o!& $@!' >$@;\
	[ -s $@ ] || rm -f $@

%.o: %.cc
	@echo "Compiling $<"
	$(CXX) -c $(CXXFLAGS) $<

%.o: %.C
	@echo "Compiling $<"
	$(CXX) -c $(CXXFLAGS) $<

dictclean:
	@rm -f $(LIB)Dict.*

clean:
	@rm -f $(LIBDIR)/lib$(LIB).so	
	@rm -f $(LIBDIR)/lib$(LIB).rootmap
	@rm -f $(LIB)Dict.*
	@rm -f $(INCDIR)/*Dict.h
	@rm -f *.o
	@rm -f *~ *.d
#	@rm -f $(addprefix $(INCDIR)/,$(addsuffix .hh,$(MODULES))) $(INCDIR)/hdf2root.hh $(INCDIR)/LinkDef.hh $(INCDIR)/PlasmaGlob.hh

distclean: clean

