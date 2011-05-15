include ../Make.common

CFLAGS += -g -Wall

# # Linux
# INC    = -I/usr/include
# LIBS   = -L/usr/lib -lmpfr -lgmp -lgmpxx -lm


# Mac
#INC    = -I/sw/include
#LIBS   = -L/sw/lib -lmpfr -lgmp -lgmpxx -lm

# Above here should be included from Make.common


headers = $(wildcard *.h)


LIBNAME = libgmpfrxx.a
# Mac OS dynamic linking
#sharedlibrary: gmpfrxx.o mpfr_mul_d.o
#	$(CXX) -dynamiclib -Wl,-undefined,dynamic_lookup -o libgmpfrxx.so $^

# Linux
ifeq ($(findstring linux,$(hostos)),linux)
$(LIBNAME): gmpfrxx.o mpfr_mul_d.o
	ar rv $(LIBNAME) $^
endif

# Mac OS static linking
ifeq ($(findstring darwin,$(hostos)),darwin)
$(LIBNAME): gmpfrxx.o mpfr_mul_d.o
	libtool -static -o $(LIBNAME) $^
endif

echo:
	@echo "GMP_INCLUDE=$(GMP_INCLUDE)"

example : example.o gmpfrxx.o mpfr_mul_d.o
	$(CXX) -o $@ $(CFLAGS) $^ $(GMP_LIBS) $(MPFR_LIBS)

%.o : %.f
	$(FC) -c $(FFLAGS) $<

%.o : %.cc $(headers)
	$(CXX) -c $(CFLAGS) $(GMP_INCLUDE) $(MPFR_INCLUDE) $<

%.o : %.cpp $(headers)
	$(CXX) -c $(CFLAGS) $(GMP_INCLUDE) $(MPFR_INCLUDE) $<

%.o : %.c $(headers)
	$(CC) -c $(CFLAGS) $(GMP_INCLUDE) $(MPFR_INCLUDE) $<

files = README example.cpp gmpfrxx.h gmpfrxx.cpp mpfr_mul_d.h \
	mpfr_mul_d.c Makefile

.PHONY : always
gmpfrxx.zip : always
	cp -a $(files) gmpfrxx
	zip -r gmpfrxx gmpfrxx

.PHONY : clean
clean:
	rm -f *.o example *.a *.so *~
