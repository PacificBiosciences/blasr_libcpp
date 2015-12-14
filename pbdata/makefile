all:

THISDIR:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
-include ${CURDIR}/defines.mk
include ${THISDIR}/../rules.mk

CXXOPTS  += -std=c++11 -pedantic -Wall -Wextra
# CURDIR should have libconfig.h
INCLUDES += ${CURDIR}
INCLUDES += ${THISDIR} matrix reads qvs metagenome saf utils
INCLUDES += ${LIBBLASR_INC}
SYSINCLUDES = ${PBBAM_INC} ${HTSLIB_INC} ${BOOST_INC}
LIBS     += ${PBBAM_LIB} ${HTSLIB_LIB}
LDFLAGS  += $(patsubst %,-L%,${LIBS})

all: libpbdata.a libpbdata${SH_LIB_EXT}

paths := . matrix reads metagenome qvs saf utils loadpulses alignment amos sam
paths := $(patsubst %,${THISDIR}%,${paths})
sources := $(shell find ${THISDIR} -name '*.cpp')
sources := $(notdir ${sources})
objects := $(sources:.cpp=.o)
shared_objects := $(sources:.cpp=.shared.o)
dependencies := $(objects:.o=.d) $(shared_objects:.o=.d)

vpath %.cpp ${paths}

libpbdata.a: $(objects)
	$(AR) $(ARFLAGS) $@ $^

libpbdata${SH_LIB_EXT}: $(shared_objects)

libconfig.h:
	cp -af ${LIBCONFIG_H} $@

clean: 
	rm -f libpbdata.a  libpbdata.so *.o *.d

-include $(dependencies)
depend: $(dependencies:.d=.depend)
