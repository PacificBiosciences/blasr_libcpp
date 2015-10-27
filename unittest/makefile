SHELL=bash
THISDIR:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
SRCDIR:=${THISDIR}
-include ${CURDIR}/defines.mk
include ${THISDIR}/../rules.mk

MKDIR := mkdir

null :=
space := $(null) $(null)

broken_test_sources := \
	${SRCDIR}/alignment/format/SAMHeaderPrinter_gtest.cpp \
	$(null)


gtest_sources  := $(GTEST_SRCDIR)/gtest/gtest-all.cc \
                  $(GTEST_SRCDIR)/gtest/gtest_main.cc \
		  $(null)

test_sources   := $(wildcard ${SRCDIR}/pbdata/*.cpp) \
                  $(wildcard ${SRCDIR}/pbdata/utils/*.cpp) \
                  $(wildcard ${SRCDIR}/pbdata/metagenome/*.cpp) \
                  $(wildcard ${SRCDIR}/pbdata/saf/*.cpp) \
                  $(wildcard ${SRCDIR}/pbdata/reads/*.cpp) \
                  $(wildcard ${SRCDIR}/pbdata/qvs/*.cpp)  \
                  \
                  $(wildcard ${SRCDIR}/hdf/*.cpp) \
                  \
                  $(wildcard ${SRCDIR}/alignment/*.cpp) \
                  $(wildcard ${SRCDIR}/alignment/utils/*.cpp) \
                  $(wildcard ${SRCDIR}/alignment/query/*.cpp) \
                  $(wildcard ${SRCDIR}/alignment/datastructures/alignment/*.cpp) \
                  $(wildcard ${SRCDIR}/alignment/files/*.cpp) \
                  $(wildcard ${SRCDIR}/alignment/format/*.cpp) \
                  $(null)

# Remove broken tests from the test_sources list
test_sources   := $(filter-out $(broken_test_sources),$(test_sources))

paths := alignment alignment/files alignment/datastructures/alignment alignment/utils alignment/format \
	pbdata pbdata/utils pbdata/metagenome pbdata/saf pbdata/reads pbdata/qvs \
	hdf alignment/query
paths := $(patsubst %,${SRCDIR}%,${paths}) ${GTEST_SRCDIR}/gtest
sources   := $(gtest_sources) $(test_sources)
sources   := $(notdir ${sources})
objects   := $(patsubst %.cc,%.o,$(filter %.cc,$(sources))) \
             $(patsubst %.cpp,%.o,$(filter %.cpp,$(sources))) \
             $(null)
dependencies:=$(objects:%.o=%.d)


INCLUDES+= \
    ${SRCDIR} \
    $(GTEST_INC) \
    $(LIBBLASR_INC) \
    $(LIBPBIHDF_INC) \
    $(LIBPBDATA_INC) \
    $(PBBAM_INC) \
    $(HTSLIB_INC) \
    $(HDF5_INC) \
    $(BOOST_INC) \
    $(null)

LIBS+= \
    $(LIBBLASR_LIB) \
    $(LIBPBIHDF_LIB) \
    $(LIBPBDATA_LIB) \
    $(PBBAM_LIB) \
    $(HTSLIB_LIB) \
    $(HDF5_LIB) \
    $(HDF5_CPP_LIB) \
    $(ZLIB_LIB) \
    $(GCC_LIB) \
    $(null)

ldlibs     := -lblasr -lpbihdf -lpbdata -lpbbam -lhts -lhdf5_cpp -lhdf5 -lz
sys_ldlibs := -lpthread -ldl -lrt

cxxopts  := -std=c++11 -Wno-div-by-zero 
cxxflags := -O3
cppflags := $(patsubst %,-I%,${includes})
ldflags  := $(patsubst %,-L%,${LIBS}) $(sys_ldflags)



override CPPFLAGS := $(cppflags) $(CPPFLAGS)
override CXXFLAGS := $(cxxflags) $(cxxopts) $(CXXFLAGS)
override LDLIBS   := $(ldlibs) $(sys_ldlibs) $(LDLIBS)
override LDFLAGS  := $(ldflags) $(LDFLAGS)

COMPILE.cpp   = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
COMPILE.cc    = $(COMPILE.cpp)
LINK.o        = $(CXX) $(LDFLAGS) $(TARGET_ARCH)

vpath %.cpp ${paths}
vpath %.cc ${paths}

all: test-runner

test-runner: $(objects)
	$(LINK.o) $^ $(LDLIBS) -o $@

gtest: test-runner
	LD_LIBRARY_PATH=$(subst $(space),:,$(strip $(LIBS))) ./$< --gtest_output=xml:./xml/all.xml

# Build objects
%.o: %.cpp
	$(COMPILE.cpp) -o $@ $<

%.o: %.cc
	$(COMPILE.cc) -o $@ $<

clean:
	$(RM) -r $(OUTDIR) *.o test-runner

-include ${dependencies}
depend: $(dependencies:.d=.depend)
