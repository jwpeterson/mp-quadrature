# Include common Make variable definitions.
include Make.common

# source files
srcfiles 	:= $(wildcard src/*.C)

# object files
objects		:= $(patsubst %.C, %.o, $(srcfiles))

# driver source codes
driver_src	:= $(wildcard drivers/*.C)

# driver executable files
drivers		:= $(patsubst %.C, %, $(driver_src))


LIBNAME = libmpquad.a
MPQ_LIB = -L./lib -lmpquad

ALL_INCLUDES=$(GMP_INCLUDE) $(MPFR_INCLUDE) $(GMPFRXX_INCLUDE) $(MPQ_INCLUDE)
ALL_LIBS=$(MPQ_LIB) $(GMP_LIBS) $(MPFR_LIBS) $(GMPFRXX_LIBS) 

# Flags to turn on extra debugging and print routines.
#EXTRA_FLAGS=-g -DDEBUG
EXTRA_FLAGS=

all: 
	make -C gmpfrxx
	make ./lib/$(LIBNAME)
	make $(drivers)


# Static linking under linux
ifeq ($(findstring linux,$(hostos)),linux)
./lib/$(LIBNAME): $(objects)
	ar rv ./lib/$(LIBNAME) $^
endif

# Mac OS static linking
ifeq ($(findstring darwin,$(hostos)),darwin)
./lib/$(LIBNAME): $(objects)
	libtool -static -o ./lib/$(LIBNAME) $^
endif


%.o: %.C
	$(CXX) $(EXTRA_FLAGS) $(ALL_INCLUDES) -c $< -o $@ 

# Link target for driver programs.  
drivers/%: drivers/%.o
	$(CXX) $(EXTRA_FLAGS) $(ALL_INCLUDES) $< -o $@ $(ALL_LIBS)

echo:
	@echo $(objects)
#make -C gmpfrxx echo
#@echo $(GMP_LIBS)
#@echo $(GMP_INCLUDE)

#.PHONY : .depend
# Don't name the target .depend.  This way, assuming there will never be
# a file called "depend" in this directory, the action of the depend target
# will occur every time you type 'make depend'.  The make_dependencies script
# does not appear to recognize .cc files and creates a circular dependency
# with those.  The fix is probably simple enough, but I don't know enough
# perl.
depend: 
	@perl ./make_dependencies.pl -I./include "-So" $(srcfiles)  $(driver_src) > .depend
	@echo "Updated .depend"




clean:
	make -C gmpfrxx clean
	rm -rf *~ $(objects) ./lib/$(LIBNAME) $(drivers) drivers/*.o


include .depend
