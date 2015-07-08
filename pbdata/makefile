all:

include ../rules.mk
include ../defines.mk

CXXOPTS  += -std=c++11 -pedantic -MMD -MP
CPPFLAGS += -I. -Imatrix -Ireads -Iqvs -Imetagenome -Isaf -Iutils -Ialignment
INCLUDES += ${PBBAM_INCLUDE} ${HTSLIB_INCLUDE} ${BOOST_INCLUDE}
DEP_LIBS += ${PBBAM_LIB} ${HTSLIB_LIB}
LDFLAGS  += $(patsubst %,$(dirname %),${DEP_LIBS})

all: libpbdata.a libpbdata.so

sources := $(wildcard *.cpp) \
	       $(wildcard matrix/*.cpp) \
	       $(wildcard reads/*.cpp) \
	       $(wildcard metagenome/*.cpp) \
	       $(wildcard qvs/*.cpp) \
	       $(wildcard saf/*.cpp) \
	       $(wildcard utils/*.cpp) \
	       $(wildcard loadpulses/*.cpp) \
	       $(wildcard alignment/*.cpp) \
	       $(wildcard amos/*.cpp) \
	       $(wildcard sam/*.cpp) 

objects := $(sources:.cpp=.o)
shared_objects := $(sources:.cpp=.shared.o)
dependencies := $(objects:.o=.d) $(shared_objects:.o=.d)

libpbdata.a: $(objects)
	$(AR) $(ARFLAGS) $@ $^

libpbdata.so: $(shared_objects) $(DEP_LIBS)
	$(CXX) $(LD_SHAREDFLAGS) -o $@ $^

clean: 
	@rm -f libpbdata.a  libpbdata.so
	@rm -f $(objects) $(shared_objects) $(dependencies)

-include $(dependencies)
