SHELL=/bin/bash

.PHONY: all debug profile gtest clean cleanall

# $Change: 140182 $ 

all:
	${MAKE} -C pbdata all
	${MAKE} -C hdf all
	${MAKE} -C alignment all

debug:
	${MAKE} subs

profile:
	${MAKE} subs

g:
	${MAKE} subs

gtest:
	${MAKE} -C unittest -f normal.Makefile gtest

clean:
	${MAKE} -C pbdata clean
	${MAKE} -C hdf clean
	${MAKE} -C alignment clean
	#@${MAKE} -C unittest clean

cleanall: clean

# Note: hdf/ will not build if nohdf=1 is configured.
