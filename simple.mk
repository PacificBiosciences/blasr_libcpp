SHELL          = bash
G_BUILDOS_CMD := bash -c 'set -e; set -o pipefail; id=$$(lsb_release -si | tr "[:upper:]" "[:lower:]"); rel=$$(lsb_release -sr); case $$id in ubuntu) printf "$$id-%04d\n" $${rel/./};; centos) echo "$$id-$${rel%%.*}";; *) echo "$$id-$$rel";; esac' 2>/dev/null
OS_STRING     ?= $(shell $(G_BUILDOS_CMD))

# magic for non-verbose builds
V ?= 0

CXX_0 = @echo "  CXX	$@"; $(CXX)
CXX_1 = $(CXX)
CXX_pp = $(CXX_$(V))

AR_0 = @echo "  AR	$@"; $(AR)
AR_1 = $(AR)
AR_pp = $(AR_$(V))

ARFLAGS := rc
