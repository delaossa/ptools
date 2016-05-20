include $(PTOOLS)/config/Makefile.common

EXES = pcheck pcopy pdelete hdf2root looper plotBunch doEvolutions doSnapshot plotBeamsEnergy rawConvert plotBunchL ele2osi

HDF5CFLAGS    = -I$(HDF5CPP)/include
HDF5LIBS      = -L$(HDF5CPP)/lib -L$(SZIP)/lib -lhdf5_cpp -lhdf5 -lz

CXXFLAGS      += $(HDF5CFLAGS) -I$(INCDIR) -I.
LIBS          := $(ROOTLIBS) -lGpad -lGraf $(HDF5LIBS) -L$(LIBDIR) -lptools

.PHONY: bin clean 

all: bin 

bin: $(addsuffix .e,$(EXES))

%.o: %.cc
	@echo "Compiling $@"
	$(CXX) -c $(CXXFLAGS) $< 

%.e: %.o
	@mkdir -p $(BINDIR)
	@echo "Linking $<"
	$(LD) $(LDFLAGS) -o $@ $< $(LIBS) $(H2RLIBS)
	ln -sf $(EXESDIR)/$@ $(BINDIR)/$(subst .e,,$@)

clean:
	@echo "Cleaning..."
	@rm -f *.o
	@rm -f *~ *.d
	@for file in $(EXES); do (rm -f $(BINDIR)/$$file); done
	@for file in $(EXES); do (rm -f $(EXESDIR)/$$file.e); done

