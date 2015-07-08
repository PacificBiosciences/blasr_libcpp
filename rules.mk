ARFLAGS         := rc
CXX_SHAREDFLAGS := -fPIC
LD_SHAREDFLAGS  := -shared -fPIC
CPPFLAGS        += $(patsubst %,-I%,${INCLUDES})


%.a:
	${AR} $(ARFLAGS) $@ $^

%.so:
	${CXX} ${LD_SHAREDFLAGS} ${LDFLAGS} -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(LEGACY) $(CPPFLAGS) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.o) $(@:%.o=%.d)" -c $< -o $@

%.shared.o: %.cpp
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(LEGACY) $(CPPFLAGS) $(CXX_SHAREDFLAGS) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.o) $(@:%.o=%.d)" -c $< -o $@
