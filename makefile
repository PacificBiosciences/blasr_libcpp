SHELL=/bin/bash

.PHONY: all debug profile gtest clean cleanall

# $Change: 140182 $ 

all:
	${MAKE} -C pbdata -f build.mk all
	${MAKE} -C hdf -f build.mk all
	${MAKE} -C alignment -f build.mk all

debug:
	${MAKE} subs

profile:
	${MAKE} subs

g:
	${MAKE} subs

gtest:
	${MAKE} -C unittest -f normal.Makefile gtest

clean:
	${MAKE} -C pbdata -f build.mk clean
	${MAKE} -C hdf -f build.mk clean
	${MAKE} -C alignment -f build.mk clean
	#@${MAKE} -C unittest clean

cleanall: clean

# Note: hdf/ will not build if nohdf=1 is configured.
