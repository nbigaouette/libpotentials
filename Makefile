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
HEADERS          = $(wildcard $(addsuffix *.$(HEADEXT),$(addsuffix /, $(SRCDIRS)) ) )
LANGUAGE         = CPP

# Include the generic rules
include makefiles/Makefile.rules

### Floats type: Use single precision or double precision?
### By default, it's double precision.
CFLAGS          += -DFLOATTYPE_SINGLE

#################################################################
# Project specific options

$(eval $(call Flags_template,stdcout,StdCout.hpp,ssh://optimusprime.selfip.net/git/nicolas/stdcout.git))
$(eval $(call Flags_template,memory,Memory.hpp,ssh://optimusprime.selfip.net/git/nicolas/memory.git))

# Project is a library. Include the makefile for build and install.
include makefiles/Makefile.library

.PHONY: version
version: src/Version.hpp
src/Version.hpp: force
	echo "#ifndef INC_LIBPOTENTIALS_VERSION_hpp" > src/Version.hpp
	echo "#define INC_LIBPOTENTIALS_VERSION_hpp" >> src/Version.hpp
	echo "namespace libpotentials {" >> src/Version.hpp
	echo "    const char *const build_time = \"`date`\";" >> src/Version.hpp
	echo "    const char *const build_sha = \"$(GIT_BRANCH)\";" >> src/Version.hpp
	echo "    const char *const build_branch = \"`$(GIT) rev-parse HEAD`\";" >> src/Version.hpp
	echo "}" >> src/Version.hpp
	echo "#endif // #ifndef INC_LIBPOTENTIALS_VERSION_hpp" >> src/Version.hpp

############ End of file ########################################
