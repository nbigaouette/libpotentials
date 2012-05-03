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
BIN              = $(LIB)
SRCDIRS          = src
SRCEXT           = cpp
HEADEXT          = hpp
TESTDIRS         = unit_testing
HEADERS          = $(wildcard $(addsuffix *.$(HEADEXT),$(addsuffix /, $(SRCDIRS)) ) )

LANGUAGE         = CPP
output_dir       = output

# Include the generic rules
include makefiles/Makefile.rules

#################################################################
# Project specific options

### Floats type: Use single precision or double precision?
### By default, it's double precision.
CFLAGS          += -DFLOATTYPE_SINGLE

### Use cubic splines to smooth the HermanSkillman field/potential?
CFLAGS          += -DHSSPLINE

$(eval $(call Flags_template,stdcout,StdCout.hpp,ssh://optimusprime.selfip.net/git/nicolas/stdcout.git))
$(eval $(call Flags_template,memory,Memory.hpp,ssh://optimusprime.selfip.net/git/nicolas/memory.git))
$(eval $(call FLAGS_template,assert,Assert.hpp,ssh://optimusprime.selfip.net/git/nicolas/assert.git))

# Project is a library. Include the makefile for build and install.
include makefiles/Makefile.library


############ End of file ########################################
