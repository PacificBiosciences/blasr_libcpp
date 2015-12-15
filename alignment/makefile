all:

THISDIR:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
-include ${CURDIR}/defines.mk
include ${THISDIR}/../rules.mk

CXXOPTS  := -std=c++11 -pedantic -Wno-long-long -Wall -Wextra
INCLUDES += ${THISDIR} ${LIBPBIHDF_INC} ${LIBPBDATA_INC}
SYSINCLUDES = ${HDF5_INC} ${PBBAM_INC} ${HTSLIB_INC} ${BOOST_INC}
LIBS     += ${LIBPBDATA_LIB} ${LIBPBIHDF_LIB} ${HDF5_LIB} ${PBBAM_LIB} ${HTSLIB_LIB} ${ZLIB_LIB}
LDFLAGS  += $(patsubst %,-L%,${LIBS})
LDLIBS += -lpbdata

ifeq (${nohdf},)
  LDLIBS+= -lpbihdf
  #LDFLAGS+= -flat_namespace # so we do not need LDLIBS+= -lhdf5 -lhdf5_cpp
endif
# We might also need some -l* for pbbam, etc.

all: libblasr.a libblasr${SH_LIB_EXT}

paths := . simulator format files utils tuples statistics qvs suffixarray \
	datastructures/alignment datastructures/alignmentset datastructures/anchoring datastructures/tuplelists \
	algorithms/alignment algorithms/alignment/sdp algorithms/anchoring algorithms/compare algorithms/sorting \
	query
paths := $(patsubst %,${THISDIR}%,${paths})
sources := $(shell find ${THISDIR} -name '*.cpp')

ifdef nohdf
sources := $(filter-out ${THISDIR}files/% ${THISDIR}utils/FileOfFileNames.cpp ${THISDIR}format/SAMHeaderPrinter.cpp, $(sources))
endif

sources := $(notdir ${sources})
objects := $(sources:.cpp=.o)
shared_objects := $(sources:.cpp=.shared.o)
dependencies := $(objects:.o=.d) $(shared_objects:.o=.d)

vpath %.cpp ${paths}

libblasr.a: $(objects)
	$(AR) $(ARFLAGS) $@ $^

libblasr${SH_LIB_EXT}: $(shared_objects)

clean: 
	rm -f libblasr.a libblasr.so *.o *.d

-include $(dependencies)
depend: $(dependencies:.d=.depend)
