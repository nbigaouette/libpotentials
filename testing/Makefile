#################################################################
# Main makefile
# Please edit this makefile to adapt to your project.
# Type "make help" for usage
#################################################################

# To run with a maximum of 500 MiB:
# softlimit -a 500000000 ./mdgit
# softlimit is part of http://cr.yp.to/daemontools.html

# To check memory usage:
# while [ 1 ]; do pmap -d `pidof project_name` | grep writeable | sed "s|K||g" | awk '{print ""$4" KiB    "$4/1024" MiB    "$4/1024/1024" GiB"}' ; sleep 0.1 ;done 2> /dev/null

# Project options
# Intel ICC sets LIB. Reset it here to make sure it's empty. Else code won't compile.
LIB             :=
BIN              = libpotential_testing
SRCDIRS          = src
TESTDIRS         = unit_testing
SRCEXT           = cpp
HEADEXT          = hpp
# LANGUAGE         = C
LANGUAGE         = CPP

# Include the generic rules
include makefiles/Makefile.rules

.PHONY: s show
s: show
show:
	./plot.py

#################################################################
# Project specific options
CFLAGS          +=
LDFLAGS         +=

### Floats type: Use single precision or double precision?
### By default, it's double precision.
CFLAGS          += -DFLOATTYPE_SINGLE

LINK_PREFERED    = shared

# To include a library:
$(eval $(call Flags_template,stdcout,StdCout.hpp,ssh://optimusprime.selfip.net/git/nicolas/stdcout.git))
$(eval $(call Flags_template,potentials,LibPotentials.hpp,ssh://optimusprime.selfip.net/git/nicolas/libpotentials.git))
$(eval $(call Flags_template,memory,Memory.hpp,ssh://optimusprime.selfip.net/git/nicolas/memory.git))
# Or just the header(s):
# $(eval $(call CFLAGS_template,stdcout,StdCout.hpp,ssh://optimusprime.selfip.net/git/nicolas/stdcout.git))
# Or just the library:
# $(eval $(call LDFLAGS_template,stdcout,StdCout.hpp,ssh://optimusprime.selfip.net/git/nicolas/stdcout.git))

# Include the library stuff
# include makefiles/Makefile.library

############ End of file ########################################
