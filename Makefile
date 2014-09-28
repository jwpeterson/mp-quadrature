# Include common Make variable definitions.
include Make.common

# source files
srcfiles 	:= $(wildcard src/*.C)

# header files
headerfiles 	:= $(wildcard include/*.h)

# object files
objects		:= $(patsubst %.C, %.o, $(srcfiles))

# src dependency files
src_depend	:= $(patsubst %.C, %.d, $(srcfiles))

# driver source codes
driver_src	:= $(wildcard drivers/*.C)

# driver executable files
drivers		:= $(patsubst %.C, %, $(driver_src))

# driver dependency files
drivers_depend	:= $(patsubst %.C, %.d, $(driver_src))


LIBNAME = libmpquad.a
MPQ_LIB = -L./lib -lmpquad

ALL_INCLUDES=$(GMP_INCLUDE) $(MPFR_INCLUDE) $(GMPFRXX_INCLUDE) $(MPQ_INCLUDE)
ALL_LIBS=$(MPQ_LIB) $(GMP_LIBS) $(MPFR_LIBS) $(GMPFRXX_LIBS) 

# Flags to turn on extra debugging and print routines.
#EXTRA_FLAGS=-g -DDEBUG
EXTRA_FLAGS=-Wall

all: 
	make -C gmpfrxx
	make ./lib/$(LIBNAME)
	make $(drivers)


# Static linking under linux
ifeq ($(findstring linux,$(hostos)),linux)
./lib/$(LIBNAME): $(objects)
	ar rv ./lib/$(LIBNAME) $^
endif

# Mac OS static linking.  Make sure to use the libtool in /usr/bin to avoid
# getting a GNU libtool that might otherwise be in your PATH.
ifeq ($(findstring darwin,$(hostos)),darwin)
./lib/$(LIBNAME): $(objects)
	mkdir -p lib
	/usr/bin/libtool -static -o ./lib/$(LIBNAME) $^
endif

# -MMD Like -MD except mention only user header files, not system header files.
# -MP  This option instructs CPP to add a phony target for each
#      dependency other than the main file, causing each to depend on
#      nothing.  These dummy rules work around errors make gives if you
#      remove header files without updating the Makefile to match.
# -MF  When used with -M or -MM, specifies a file to write the dependencies to.
# -MT  target = Change the target of the rule emitted by dependency
#      generation.  By default CPP takes the name of the main input file,
#      deletes any directory components and any file suffix such as .c,
#      and appends the platform's usual object suffix.  The result is the
%.o: %.C
	$(CXX) -MMD -MP $(EXTRA_FLAGS) $(ALL_INCLUDES) -c $< -o $@

# Link target for driver programs.  Because we currently build a
# static library, the executables all have an explicit dependence on
# the library.
drivers/%: drivers/%.o ./lib/$(LIBNAME)
	$(CXX) -MMD -MP $(EXTRA_FLAGS) $(ALL_INCLUDES) $< -o $@ $(ALL_LIBS)

echo:
	@echo $(src_depend) $(drivers_depend)
#@echo $(objects)
#make -C gmpfrxx echo
#@echo $(GMP_LIBS)
#@echo $(GMP_INCLUDE)

clean:
	make -C gmpfrxx clean
	rm -rf *~ $(objects) ./lib/$(LIBNAME) $(drivers) drivers/*.o $(src_depend) $(drivers_depend)

# Include dependency rules we generated for all the sources
-include $(src_depend) $(drivers_depend)
