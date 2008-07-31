# Include common Make variable definitions.
include Make.common

# source files
srcfiles 	:= $(wildcard src/*.C)

# object files
objects		:= $(patsubst %.C, %.o, $(srcfiles))

JACOBI_LIBNAME = libjacobi.a
JACOBI_LIB = -L. -ljacobi

ALL_INCLUDES=$(GMP_INCLUDE) $(GMPFRXX_INCLUDE) $(MPQ_INCLUDE)
ALL_LIBS=$(JACOBI_LIB) $(GMP_LIBS) $(GMPFRXX_LIBS) 

# Flags to turn on extra debugging and print routines.
#EXTRA_FLAGS=-g -DDEBUG
EXTRA_FLAGS=

all: library
	make -C gmpfrxx library
	$(CXX) $(EXTRA_FLAGS) $(ALL_INCLUDES) jacobi_rule.cc -o jacobi_rule $(ALL_LIBS)

# Static linking under linux
ifeq ($(findstring linux,$(hostos)),linux)
library: $(objects)
	ar rv $(JACOBI_LIBNAME) $^
endif

# Mac OS static linking
ifeq ($(findstring darwin,$(hostos)),darwin)
library: $(objects)
	libtool -static -o $(JACOBI_LIBNAME) $^
endif


%.o: %.C
	$(CXX) $(EXTRA_FLAGS) $(ALL_INCLUDES) -c $< -o $@ 

echo:
	@echo $(objects)
#make -C gmpfrxx echo
#@echo $(GMP_LIBS)
#@echo $(GMP_INCLUDE)

print_gauss: print_gauss.cc
	$(CXX) $(ALL_INCLUDES) $< -o $@ $(ALL_LIBS)

conical_product: conical_product.cc
	$(CXX) $(ALL_INCLUDES) $< -o $@ $(ALL_LIBS)


clean:
	make -C gmpfrxx clean
	rm -rf *~ jacobi_rule $(objects) $(JACOBI_LIBNAME) print_gauss
