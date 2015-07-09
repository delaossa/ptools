include $(PTOOLS)/config/Makefile.common

MACROS = $(wildcard Plot*.C) $(wildcard Do*.C) $(wildcard Get*.C)

HDF5CFLAGS    = -I$(HDF5CPP)/include
HDF5LIBS      = -L$(HDF5CPP)/lib -L$(SZIP)/lib -lhdf5_cpp -lhdf5 -lz

CXXFLAGS      += $(HDF5CFLAGS) -I$(INCDIR) -I.
LIBS          := $(ROOTLIBS) -lGpad -lGraf $(HDF5LIBS) -L$(LIBDIR) -lptools

all: $(MACROS:.C=.o)

%.o: %.C
	@echo "Compiling .C $<"
	$(CXX) -c $(CXXFLAGS) $<

clean:
	@rm -f *.o

distclean: clean
