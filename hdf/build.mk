all:

include ../rules.mk
include ../defines.mk

CXXOPTS  += -std=c++11 -pedantic -MMD -MP
INCLUDES += ${LIBPBDATA_INCLUDE} ${HDF5_INCLUDE} ${PBBAM_INCLUDE} ${HTSLIB_INCLUDE} ${BOOST_INCLUDE}
DEP_LIBS += ${LIBPBDATA_LIB} ${HDF5_LIB} ${PBBAM_LIB} ${HTSLIB_LIB} ${ZLIB_LIB}
LDFLAGS  += $(patsubst %,$(dirname %),${DEP_LIBS})

all: libpbihdf.a libpbihdf.so

sources := $(wildcard *.cpp)
objects := $(sources:.cpp=.o)
shared_objects := $(sources:.cpp=.shared.o)
dependencies := $(objects:.o=.d) $(shared_objects:.o=.d)

libpbihdf.a: $(objects)
	$(AR) $(ARFLAGS) $@ $^

libpbihdf.so: $(shared_objects) $(DEP_LIBS)
	$(CXX) $(LD_SHAREDFLAGS) -o $@ $^

clean: 
	@rm -f libpbihdf.a libpbihdf.so
	@rm -f $(objects) $(shared_objects) $(dependencies)

-include $(dependencies)
