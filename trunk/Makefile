# Include common Make variable definitions.
include Make.common

# source files
srcfiles 	:= $(wildcard src/*.C)

# object files
objects		:= $(patsubst %.C, %.o, $(srcfiles))

# driver source codes
driver_src	:= $(wildcard drivers/*.cc)

# driver executable files
drivers		:= $(patsubst %.cc, %, $(driver_src))


LIBNAME = libmpquad.a
MPQ_LIB = -L./lib -lmpquad

ALL_INCLUDES=$(GMP_INCLUDE) $(GMPFRXX_INCLUDE) $(MPQ_INCLUDE)
ALL_LIBS=$(MPQ_LIB) $(GMP_LIBS) $(GMPFRXX_LIBS) 

# Flags to turn on extra debugging and print routines.
#EXTRA_FLAGS=-g -DDEBUG
EXTRA_FLAGS=

all: 
	make -C gmpfrxx library
	make library
	make $(drivers)


# Static linking under linux
ifeq ($(findstring linux,$(hostos)),linux)
library: $(objects)
	ar rv ./lib/$(LIBNAME) $^
endif

# Mac OS static linking
ifeq ($(findstring darwin,$(hostos)),darwin)
library: $(objects)
	libtool -static -o ./lib/$(LIBNAME) $^
endif


%.o: %.C
	$(CXX) $(EXTRA_FLAGS) $(ALL_INCLUDES) -c $< -o $@ 

drivers/%: drivers/%.cc
	$(CXX) $(EXTRA_FLAGS) $(ALL_INCLUDES) $< -o $@ $(ALL_LIBS)

echo:
	@echo $(objects)
#make -C gmpfrxx echo
#@echo $(GMP_LIBS)
#@echo $(GMP_INCLUDE)



clean:
	make -C gmpfrxx clean
	rm -rf *~ $(objects) ./lib/$(LIBNAME) $(drivers)
