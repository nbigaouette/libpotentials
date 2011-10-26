#################################################################
# Main makefile
# Type "make help" for usage
#
# To compile (optimized) and install static and shared library
# with each avaible compilers, type:
# for c in pgi intel gcc 'sun studio12'; do make c $c optimized static shared install; done
#################################################################

# Project options
LIB              = potentials
SRCDIRS          = src
SRCEXT           = cpp
HEADEXT          = hpp
LANGUAGE         = CPP

# Include the generic rules
include makefiles/Makefile.rules

### Floats type: Use single precision or double precision?
### By default, it's double precision.
CFLAGS          += -DFLOATTYPE_SINGLE

# Don't install all headers

#################################################################
# Project specific options

$(eval $(call Flags_template,stdcout,StdCout.hpp,ssh://optimusprime.selfip.net/git/nicolas/stdcout.git))
$(eval $(call Flags_template,memory,Memory.hpp,ssh://optimusprime.selfip.net/git/nicolas/memory.git))

# Project is a library. Include the makefile for build and install.
include makefiles/Makefile.library


############ End of file ########################################
