# Include common Make variable definitions.
include Make.common

# source files
srcfiles 	:= $(wildcard src/*.C)

# header files
headerfiles 	:= $(wildcard include/*.h)

# object files
objects		:= $(patsubst %.C, %.o, $(srcfiles))

# libtool object (.lo) files
libtool_objects	:= $(patsubst %.C, %.lo, $(srcfiles))

# src dependency files
src_depend	:= $(patsubst %.C, %.d, $(srcfiles))

# driver source codes
driver_src	:= $(wildcard drivers/*.C)

# driver libtool object files
driver_libtool_objects := $(patsubst %.C, %.lo, $(driver_src))

# driver executable files
drivers		:= $(patsubst %.C, %, $(driver_src))

# driver dependency files
drivers_depend	:= $(patsubst %.C, %.d, $(driver_src))

# An absolute PATH to the current directory.  This is needed for the -rpath stuff that libtool does.
MPQ_DIR ?= $(shell pwd)

LIBNAME = libmpquad
MPQ_LIB = -L./lib -lmpquad

# FIXME: configure for nlopt.
NLOPT_DIR := /home/jwpeterson/software/libmesh_install/nlopt
NLOPT_INCLUDE := -I$(NLOPT_DIR)/include
NLOPT_LIBS := -L$(NLOPT_DIR)/lib -lnlopt

ALL_INCLUDES=$(GMP_INCLUDE) $(MPFR_INCLUDE) $(GMPFRXX_INCLUDE) $(MPQ_INCLUDE) $(NLOPT_INCLUDE)

# A note on static library linking (http://stackoverflow.com/questions/45135/linker-order-gcc)
# If any [static] library A depends on symbols defined in library B,
# then library A should appear first in the list supplied to the
# linker.
ALL_LIBS=$(MPQ_LIB) $(GMPFRXX_LIBS) $(MPFR_LIBS) $(GMP_LIBS) $(NLOPT_LIBS)

# Flags to turn on extra debugging and print routines.
#EXTRA_FLAGS=-g -DDEBUG
EXTRA_FLAGS=-Wall -O2 -std=c++11

# The plus signs on these recursive make calls were required for Linux, otherwise I was getting:
# warning: jobserver unavailable: using -j1.  Add `+' to parent make rule.
# Not sure how standard this is beyond GNU make...
all:
	+make -C gmpfrxx
	+make $(drivers)


# Use GNU libtool for linking.  Note: if the name of the output file
# is .la, it will try to build a libtool library archive, but you can
# also give it a .a and it won't try to build a shared library at all.
# It seems that you must use libtool objects (.lo files) when linking
# with libtool otherwise libtool generates an errant command which
# doesn't make any sense:
#
# ar cru /Users/petejw/projects/mp-quadrature/lib/.libs/libmpquad.a
#
# i.e. there have to be some object files following the name of the
# archive for this to work.
$(MPQ_DIR)/lib/$(LIBNAME).la: $(libtool_objects)
	@echo "Linking Library "$@"..."
	@mkdir -p lib
	@./libtool --quiet --tag=CXX --mode=link $(CXX) -o $(MPQ_DIR)/lib/$(LIBNAME).la $(libtool_objects) -rpath $(MPQ_DIR)/lib
	@./libtool --quiet --mode=install install -c ./lib/$(LIBNAME).la $(MPQ_DIR)/lib

# -MMD Like -MD except mention only user header files, not system header files.
# -MP  This option instructs CPP to add a phony target for each
#      dependency other than the main file, causing each to depend on
#      nothing.  These dummy rules work around errors make gives if you
#      remove header files without updating the Makefile to match.
# -MF  When used with -M or -MM, specifies a file to write the dependencies to.
# -MT  target = Change the target of the rule emitted by dependency
#      generation.  By default CPP takes the name of the main input file,
#      deletes any directory components and any file suffix such as .c,
#      and appends the platform's usual object suffix.
%.lo: %.C
	@echo "Compiling $<..."
	@./libtool --quiet --tag=CXX --mode=compile $(CXX) -MMD -MP -MT $@ $(EXTRA_FLAGS) $(ALL_INCLUDES) -c $< -o $@

# Link target for driver programs.
drivers/%: drivers/%.lo $(MPQ_DIR)/lib/$(LIBNAME).la
	@echo "Compiling driver program $@..."
	@./libtool --quiet --tag=CXX --mode=link $(CXX) -MMD -MP -MT $@ $(EXTRA_FLAGS) $(ALL_INCLUDES) -o $@ $< $(ALL_LIBS)

echo:
	@echo $(src_depend) $(drivers_depend)
#@echo $(objects)
#make -C gmpfrxx echo
#@echo $(GMP_LIBS)
#@echo $(GMP_INCLUDE)

clean:
	@make -C gmpfrxx clean
	@echo "Cleaning in mp-quadrature..."
	@rm -rf *~ src/*~ $(objects) ./lib/* ./lib/.libs src/.libs drivers/.libs $(drivers) drivers/*.o $(src_depend) $(drivers_depend) $(libtool_objects) $(driver_libtool_objects)

# Include dependency rules we generated for all the sources
-include $(src_depend) $(drivers_depend)
