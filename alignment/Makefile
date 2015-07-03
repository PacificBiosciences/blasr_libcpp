
include ../common.mk

# To enable building a shared library, invoke as "make SHARED_LIB=true ..."
ifneq ($(SHARED_LIB),)
    # Generating shared library
    CXX_SHAREDFLAGS := -fPIC
    LD_SHAREDFLAGS  := -shared -fPIC
    TARGET_LIB      := libblasr.so
    # Developers should set these to appropriate defaults (other systems 
    # will override these on the command line):
    HDF5_LIB        := ../../../../prebuilt.out/prebuilt.out/hdf5/hdf5-1.8.12/centos-5/lib/libhdf5.so
    ZLIB_LIB        := ../../../../prebuilt.tmpsrc/zlib/zlib_1.2.8/_output/install/lib/libz.so
    HTSLIB_LIB      := ../../../staging/PostPrimary/pbbam/_output/install-build/lib/libpbbam.so
    PBBAM_LIB       := ../../../staging/PostPrimary/pbbam/third-party/htslib/_output/install-build/lib/libhts.so
    LIBPBDATA_LIB   := ../../../staging/PostPrimary/pbbam/third-party/htslib/_output/install-build/lib/libhts.so
else
    # Generating shared library
    CXX_SHAREDFLAGS :=
    LD_SHAREDFLAGS  :=
    TARGET_LIB      := libblasr.a
    HDF5_LIB        :=
    ZLIB_LIB        :=
    HTSLIB_LIB      :=
    PBBAM_LIB       :=
    LIBPBDATA_LIB   :=
endif

DEP_LIBS := $(HDF5_LIB) $(ZLIB_LIB) $(HTSLIB_LIB) $(PBBAM_LIB) $(PBDATA_LIB)

LIBPBDATA_INCLUDE := ../pbdata
LIBPBIHDF_INCLUDE := ../hdf
PBBAM_INCLUDE := $(PBBAM)/include
HTSLIB_INCLUDE ?= $(PBBAM)/third-party/htslib

INCLUDES = -I$(LIBPBDATA_INCLUDE) \
           -I$(LIBPBIHDF_INCLUDE) \
	   -I.
ifneq ($(HDF5_INC),)
INCLUDES += -I$(HDF5_INC)
else
HDF_HEADERS := hdf5-1.8.12-headers
INCLUDES += -I../hdf/$(HDF_HEADERS)/src -I../hdf/$(HDF_HEADERS)/c++/src
endif

ifneq ($(ZLIB_ROOT), notfound)
	INCLUDES += -I$(ZLIB_ROOT)/include
endif

ifeq ($(origin nopbbam), undefined)
    INCLUDES += -I$(PBBAM_INCLUDE) -I$(HTSLIB_INCLUDE) -I$(BOOST_INCLUDE)
endif

CXXOPTS := -std=c++11 -pedantic -Wno-long-long -MMD -MP

sources := $(wildcard algorithms/alignment/*.cpp) \
		   $(wildcard algorithms/alignment/sdp/*.cpp) \
		   $(wildcard algorithms/anchoring/*.cpp) \
		   $(wildcard algorithms/compare/*.cpp) \
		   $(wildcard algorithms/sorting/*.cpp) \
		   $(wildcard datastructures/alignment/*.cpp) \
		   $(wildcard datastructures/alignmentset/*.cpp) \
		   $(wildcard datastructures/anchoring/*.cpp) \
		   $(wildcard datastructures/tuplelists/*.cpp) \
		   $(wildcard suffixarray/*.cpp) \
		   $(wildcard qvs/*.cpp) \
		   $(wildcard statistics/*.cpp) \
		   $(wildcard tuples/*.cpp) \
		   $(wildcard utils/*.cpp) \
		   $(wildcard files/*.cpp) \
		   $(wildcard format/*.cpp) \
		   $(wildcard simulator/*.cpp) \
		   $(wildcard *.cpp) 

ifdef nohdf
sources := $(filter-out files/% utils/FileOfFileNames.cpp, $(sources))
endif

objects := $(sources:.cpp=.o)
shared_objects := $(sources:.cpp=.shared.o)
dependencies := $(objects:.o=.d) $(shared_objects:.o=.d)

all : CXXFLAGS ?= -O3

debug : CXXFLAGS ?= -g -ggdb -fno-inline

profile : CXXFLAGS ?= -Os -pg

g: CXXFLAGS = -g -ggdb -fno-inline -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer 

all debug profile g: $(TARGET_LIB)

libblasr.a: $(objects)
	$(AR_pp) $(ARFLAGS) $@ $^

libblasr.so: $(shared_objects) $(DEP_LIBS)
	$(CXX_pp) $(LD_SHAREDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(LEGACY) $(INCLUDES) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.o) $(@:%.o=%.d)" -c $< -o $@

%.shared.o: %.cpp
	$(CXX_pp) $(CXX_SHAREDFLAGS) $(CXXOPTS) $(CXXFLAGS) $(LEGACY) $(INCLUDES) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.o) $(@:%.o=%.d)" -c $< -o $@

# .INTERMEDIATE: $(objects)

clean: 
	@rm -f libblasr.a libblasr.so
	@rm -f $(objects) $(shared_objects) $(dependencies)
	@find . -type f -name \*.o -delete
	@find . -type f -name \*.d -delete

-include $(dependencies)
