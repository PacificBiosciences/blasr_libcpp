ARFLAGS         := rc
CXX_SHAREDFLAGS := -fPIC
#LD_SHAREDFLAGS  := -dynamiclib -fPIC
CPPFLAGS        += $(patsubst %,-I%,${INCLUDES})
CPPFLAGS        += $(patsubst %,-isystem%,${SYSINCLUDES})
CFLAGS          += -fno-common
LDFLAGS         += ${EXTRA_LDFLAGS}


%.a:
	${AR} ${ARFLAGS} $@ $^

%.so:
	${CXX} -shared ${LDFLAGS} -o $@ -Wl,-soname,$@ $^ ${LDLIBS}

%.dylib:
	${CXX} -dynamiclib ${LDFLAGS} -o $@ -Wl,-install_name,$@ $^ ${LDLIBS}

%.o: %.cpp
	${CXX} ${CXXOPTS} ${CXXFLAGS} ${CPPFLAGS} -c $< -o $@

%.shared.o: %.cpp
	${CXX} ${CXXOPTS} ${CXXFLAGS} ${CPPFLAGS} ${CXX_SHAREDFLAGS} -c $< -o $@

%.depend: %.cpp
	${CXX} ${CXXOPTS} ${CXXFLAGS} ${CPPFLAGS} -MM -MP -MG -MT $(@:.depend=.o) -MF $(@:.depend=.d) $<

%.shared.depend: %.cpp
	${CXX} ${CXXOPTS} ${CXXFLAGS} ${CPPFLAGS} -MM -MP -MG -MT $(@:.depend=.o) -MF $(@:.depend=.d) $<
