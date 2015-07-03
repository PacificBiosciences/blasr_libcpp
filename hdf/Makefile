
include ../common.mk

# To enable building a shared library, invoke as "make SHARED_LIB=true ..."
ifneq ($(SHARED_LIB),)
    # Generating shared library
    CXX_SHAREDFLAGS := -fPIC
    LD_SHAREDFLAGS  := -shared -fPIC
    TARGET_LIB      := libpbihdf.so
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
    TARGET_LIB      := libpbihdf.a
    HDF5_LIB        :=
    ZLIB_LIB        :=
    HTSLIB_LIB      :=
    PBBAM_LIB       :=
    LIBPBDATA_LIB   :=
endif

DEP_LIBS := $(HDF5_LIB) $(ZLIB_LIB) $(HTSLIB_LIB) $(PBBAM_LIB) $(PBDATA_LIB)

# FIXME: remove PBDATA_INCLUDE and assign directly to LIBPBDATA_INCLUDE.
#        Use only LIBPBDATA_INCLUDE to be consistent with the libblasr Makefile
PBDATA_INCLUDE := ../pbdata
LIBPBDATA_INCLUDE := $(PBDATA_INCLUDE)
PBBAM_INCLUDE := $(PBBAM)/include
HTSLIB_INCLUDE ?= $(PBBAM)/third-party/htslib

INCLUDES = -I$(LIBPBDATA_INCLUDE)


ifeq ($(origin nopbbam), undefined)
    INCLUDES += -I$(PBBAM_INCLUDE) -I$(HTSLIB_INCLUDE) -I$(BOOST_INCLUDE)
endif

ifneq ($(ZLIB_ROOT), notfound)
	INCLUDES += -I$(ZLIB_ROOT)/include
endif

CXXOPTS := -std=c++11 -pedantic -MMD -MP
sources := $(wildcard *.cpp)
objects := $(sources:.cpp=.o)
shared_objects := $(sources:.cpp=.shared.o)
dependencies := $(objects:.o=.d) $(shared_objects:.o=.d)

ifneq ($(HDF5_INC),)
HDF_HEADERS :=
INCLUDES += -I$(HDF5_INC)
else
HDF_HEADERS := hdf5-1.8.12-headers
INCLUDES += -I./$(HDF_HEADERS)/src -I./$(HDF_HEADERS)/c++/src
endif

all : CXXFLAGS ?= -O3

debug : CXXFLAGS ?= -g -ggdb -fno-inline

profile : CXXFLAGS ?= -Os -pg

g : CXXFLAGS ?= -g -ggdb -fno-inline -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer 

all debug profile g: $(HDF_HEADERS) $(TARGET_LIB)

libpbihdf.a: $(objects)
	$(AR_pp) $(ARFLAGS) $@ $^

libpbihdf.so: $(shared_objects) $(DEP_LIBS)
	$(CXX) $(LD_SHAREDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCLUDES) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.o) $(@:%.o=%.d)" -c $< -o $@

%.shared.o: %.cpp
	$(CXX) $(CXX_SHAREDFLAGS) $(CXXOPTS) $(CXXFLAGS) $(INCLUDES) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.o) $(@:%.o=%.d)" -c $< -o $@

$(HDF_HEADERS): 
	curl -k -L https://www.dropbox.com/s/8971bcyy5o42rxb/hdf5-1.8.12-headers.tar.bz2\?dl\=0 | tar xjf - 

# .INTERMEDIATE: $(objects)

clean: 
	@rm -f libpbihdf.a libpbihdf.so
	@rm -f $(objects) $(shared_objects) $(dependencies)

-include $(dependencies)
