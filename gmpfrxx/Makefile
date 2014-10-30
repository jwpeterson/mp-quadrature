include ../Make.common

# Compile with debugging information and warnings turned on.
CFLAGS += -g -Wall

# An absolute PATH to the current directory.  This is needed for the -rpath stuff that libtool does.
GMPFRXX_DIR ?= $(shell pwd)

# The basename of the library object to build.
LIBNAME = libgmpfrxx

# The first target.  This is what get runs when you just type 'make'.
# It depends on the .lo object files
$(GMPFRXX_DIR)/$(LIBNAME).la: gmpfrxx.lo mpfr_mul_d.lo
	@echo "Linking library "$@"..."
	@../libtool --quiet --tag=CXX --mode=link $(CXX) -o $(GMPFRXX_DIR)/$(LIBNAME).la $< -rpath $(GMPFRXX_DIR)
	@../libtool --quiet --mode=install install -c $(LIBNAME).la $(GMPFRXX_DIR)

example : example.lo $(GMPFRXX_DIR)/$(LIBNAME).la
	@echo "Linking gmpfrxx example program..."
	@../libtool --quiet --tag=CXX --mode=link $(CXX) -o $@ $^ $(GMP_LIBS) $(MPFR_LIBS)

%.lo : %.cpp
	@echo "Compiling $<..."
	@../libtool --quiet --tag=CXX --mode=compile $(CXX) -MMD -MP -MT $@ -c $(CFLAGS) $(GMP_INCLUDE) $(MPFR_INCLUDE) $< -o $@

# There is a C file in this library as well
%.lo : %.c
	@echo "Compiling $<..."
	@../libtool --quiet --tag=CC --mode=compile $(CC) -MMD -MP -MT $@ -c $(CFLAGS) $(GMP_INCLUDE) $(MPFR_INCLUDE) $< -o $@

echo:
	@echo "GMP_INCLUDE=$(GMP_INCLUDE)"

.PHONY : clean

clean:
	@echo "Cleaning in gmpfrxx..."
	@rm -rf .libs *.lo example *.a *.so *.so.* *~ *.d *.dylib *.la *.o

# Include dependency files for all the C++ and C source files
-include example.d gmpfrxx.d mpfr_mul_d.d
