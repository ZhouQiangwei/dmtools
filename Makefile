#ifndef CPP
#$(error CPP variable undefined)
#endif

RPATH = $(shell pwd)
SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
# Programs built by default; dmDMR is gated on the presence of GSL.
PROGS = dmtools bam2dm dmalign bam2motif

# Canonical bam2dm implementation. dmtools invokes the bam2dm binary produced
# from this source.
BAM2DM_SRC := bam2dm.cpp

CXX ?= g++
CC ?= gcc
AR = ar
RANLIB = ranlib
CFLAGS ?= -g -w -O3 -Wsign-compare
# Ensure C++ builds inherit the same baseline optimisation/warning set.
CXXFLAGS ?= $(CFLAGS)
# Detect an appropriate C++ standard flag that is available on a wide range of
# GCC releases (4.x through 8.x and newer). GCC 5+ understands -std=gnu++11,
# while very old releases only accept the earlier -std=gnu++0x spelling. We
# compute the compiler major version and pick the most suitable flag.
CXX_VERSION_MAJOR := $(shell $(CXX) -dumpversion | cut -d. -f1)
ifeq ($(shell expr $(CXX_VERSION_MAJOR) \>= 5),1)
    CXXSTD ?= -std=gnu++11
else
    CXXSTD ?= -std=gnu++0x
endif
CXXFLAGS += $(CXXSTD)
ASAN_FLAGS = -fsanitize=address -fno-omit-frame-pointer -O1 -g
#-fopenmp
#-Wall
## changed: = instaed of ?= above 4 lines
## -fopenmp for multi-threads
LIBS = -lm -lz
EXTRA_CFLAGS_PIC = -fpic
LDFLAGS_SUB = htslib/libhts.a -lz -lpthread -llzma -lbz2 -lcurl
LDLIBS =
INCLUDES = 
BIGWIG_CFLAGS ?=
BIGWIG_LIBS ?= -lBigWig

# Create a simple test-program to check if gcc can compile with curl
#tmpfile:=$(shell mktemp --suffix=.c)
# mktemp with -t is not portable across all platforms; provide a fallback to a
# predictable temporary directory when the first form fails. The fallback is
# only evaluated when the initial invocation fails to avoid emitting multiple
# paths.
tmpfile:=$(shell mktemp -t temphaoqiaoXXXXX.c 2>/dev/null || (TMPDIR=$${TMPDIR:-/tmp}; mktemp $$TMPDIR/temphaoqiaoXXXXX.c))
$(file >$(tmpfile),#include <curl/curl.h>)
$(file >>$(tmpfile),int main() { return 0; })
#HAVE_CURL:=$(shell $(CC) $(CFLAGS) $(EXTRA_CFLAGS_PIC) $(LIBS) -lcurl $(tmpfile) -o /dev/null >/dev/null 2>&1 && echo "YES")
#$(shell rm $(tmpfile))

# Check if CC (C compiler) supports libcurl
HAVE_CURL_CC := $(shell $(CC) $(CFLAGS) $(EXTRA_CFLAGS_PIC) $(LIBS) -lcurl $(tmpfile) -o /dev/null >/dev/null 2>&1 && echo "YES")

# Check if CXX (C++ compiler) supports libcurl
HAVE_CURL_CXX := $(shell $(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS_PIC) $(LIBS) -lcurl $(tmpfile) -o /dev/null >/dev/null 2>&1 && echo "YES")
$(shell rm $(tmpfile))

ifeq ($(and $(HAVE_CURL_CC), $(HAVE_CURL_CXX)),YES)
    # If both CC and CXX support libcurl, add the library to LIBS
    LIBS += -lcurl
else
   # and if not, disable CURL specific code compilation
   CFLAGS += -DNOCURL
   #$(info LIBS: $(CFLAGS))
endif

bigwig_tmp:=$(shell mktemp -t bigwigcheckXXXXX.c 2>/dev/null || (TMPDIR=$${TMPDIR:-/tmp}; mktemp $$TMPDIR/bigwigcheckXXXXX.c))
$(file >$(bigwig_tmp),#include <bigWig.h>)
$(file >>$(bigwig_tmp),int main(){ return bwInit(1024); })
HAVE_BIGWIG := $(shell $(CC) $(CFLAGS) $(BIGWIG_CFLAGS) $(bigwig_tmp) $(LIBS) $(BIGWIG_LIBS) -o /dev/null >/dev/null 2>&1 && echo "YES")
$(shell rm -f $(bigwig_tmp))

ifeq ($(HAVE_BIGWIG),YES)
    $(info Found libBigWig; enabling bigWig compatibility tests.)
else
    $(info libBigWig not found; skipping bigWig compatibility tests. Set BIGWIG_CFLAGS/BIGWIG_LIBS to override.)
    BIGWIG_LIBS :=
endif

#ifeq ($(HAVE_CURL),YES)
#	# If yes, add the library
#	LIBS += -lcurl
#else
#	# and if not, disable CURL specific code compilation
#	CFLAGS += -DNOCURL
#endif


prefix = /usr/local
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib

# Detect GSL for dmDMR; allow disabling with WITH_GSL=0.
ifeq ($(WITH_GSL),0)
    HAVE_GSL := NO
else
    gsl_tmp:=$(shell mktemp -t gslcheckXXXXX.c 2>/dev/null || (TMPDIR=$${TMPDIR:-/tmp}; mktemp $$TMPDIR/gslcheckXXXXX.c))
    $(file >$(gsl_tmp),#include <gsl/gsl_matrix_double.h>)
    $(file >>$(gsl_tmp),int main() { gsl_matrix *m = gsl_matrix_calloc(1,1); gsl_matrix_free(m); return 0; })
    HAVE_GSL := $(shell $(CXX) $(CXXFLAGS) -lgsl -lgslcblas -lm $(gsl_tmp) -o /dev/null >/dev/null 2>&1 && echo "YES")
    $(shell rm -f $(gsl_tmp))
endif

ifeq ($(HAVE_GSL),YES)
    PROGS += dmDMR
else
    $(info GSL not found; skipping dmDMR. Install libgsl-dev and rerun make WITH_GSL=1 to enable.)
endif

.PHONY: all clean lib test doc

.SUFFIXES: .c .o .pico

all: libs lib $(PROGS)

pycd: lib $(PROGS)

lib: lib-static lib-shared

lib-static: libBinaMeth.a

lib-shared: libBinaMeth.so

doc:
	doxygen

OBJS = io.o dmValues.o dmRead.o dmStats.o dmWrite.o

.c.o:
	$(CC) -I. $(CFLAGS) $(INCLUDES) -c -o $@ $<

.c.pico:
	$(CC) -I. $(CFLAGS) $(INCLUDES) $(EXTRA_CFLAGS_PIC) -c -o $@ $<

libBinaMeth.a: $(OBJS)
	-@rm -f $@
	$(AR) -rcs $@ $(OBJS)
	$(RANLIB) $@

libBinaMeth.so: $(OBJS:.o=.pico)
	$(CC) -shared -o $@ $(OBJS:.o=.pico) $(LDLIBS) $(LIBS)
#$(CC) -shared $(LDFLAGS) -o $@ $(OBJS:.o=.pico) $(LDLIBS) $(LIBS)

test/testLocal: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) $(BIGWIG_CFLAGS) test/testLocal.c libBinaMeth.a $(LIBS) $(BIGWIG_LIBS)

test/testRemoteManyContigs: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) $(BIGWIG_CFLAGS) test/testRemoteManyContigs.c libBinaMeth.a $(LIBS) $(BIGWIG_LIBS)

test/testRemote: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) $(BIGWIG_CFLAGS) test/testRemote.c libBinaMeth.a $(LIBS) $(BIGWIG_LIBS)

test/testWrite: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) $(BIGWIG_CFLAGS) test/testWrite.c libBinaMeth.a $(LIBS) $(BIGWIG_LIBS)

dmtools: libBinaMeth.so
	$(CC) -o $@ -I. -L. $(CFLAGS) dmtools.c dmSingleCell.c dmScShrinkage.c -lBinaMeth $(LIBS) -Wl,-rpath $(RPATH) -lpthread

#bam2dm: libBinaMeth.so
#	$(CXX) -o $@ -I. -L. $(CFLAGS) bam2dm.cpp -lBinaMeth -Wl,-rpath $(RPATH) htslib/libhts.a -llzma -lbz2 -lz

dmDMR:
	$(CXX) $(CXXFLAGS) -c -o regression.o regression.cpp -lgsl -lgslcblas -lm -lz
	$(CXX) $(CXXFLAGS) -o dmDMR dmDMR.cpp regression.o -I. -L. -lBinaMeth -Wl,-rpath $(RPATH) -lgsl -lgslcblas -lm -lz

dmalign:
	$(CXX) $(CXXFLAGS) -o genome2cg genome2cg.cpp
	$(CXX) $(CXXFLAGS) -o genomebinLen genomebinLen.cpp
	$(CXX) $(CXXFLAGS) dmalign.cpp -o dmalign -lz

htslib/libhts.a:
	$(MAKE) -C htslib libhts.a

bam2dm: libBinaMeth.a htslib/libhts.a
	$(CXX) $(CXXFLAGS) -no-pie $(BAM2DM_SRC) -o bam2dm -m64 -I. libBinaMeth.a -Wl,-rpath $(RPATH) htslib/libhts.a $(LDFLAGS_SUB)

bam2dm-asan: CXXFLAGS += $(ASAN_FLAGS)
bam2dm-asan: CFLAGS += $(ASAN_FLAGS)
bam2dm-asan: LDFLAGS_SUB += -fsanitize=address
bam2dm-asan: clean libBinaMeth.a htslib/libhts.a
	$(CXX) $(CXXFLAGS) -no-pie $(BAM2DM_SRC) -o bam2dm-asan -m64 -I. libBinaMeth.a -Wl,-rpath $(RPATH) htslib/libhts.a $(LDFLAGS_SUB)

bam2motif: libBinaMeth.a htslib/libhts.a
	$(CXX) $(CXXFLAGS) -no-pie bam2motif.cpp -o bam2motif -m64 -I. libBinaMeth.a -Wl,-rpath $(RPATH) htslib/libhts.a $(LDFLAGS_SUB)

test/exampleWrite: libBinaMeth.so
	$(CC) -o $@ -I. -L. $(CFLAGS) $(BIGWIG_CFLAGS) test/exampleWrite.c -lBinaMeth $(LIBS) $(BIGWIG_LIBS) -Wl,-rpath .

test/testIterator: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) $(BIGWIG_CFLAGS) test/testIterator.c libBinaMeth.a $(LIBS) $(BIGWIG_LIBS)

test/testBinOrder: htslib/libhts.a
	$(CC) -o $@ -I. $(CFLAGS) test/testBinOrder.c htslib/libhts.a $(LDFLAGS_SUB)

test/testCoordinates: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/testCoordinates.c libBinaMeth.a $(LIBS)

test/test_id_roundtrip: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/test_id_roundtrip.c libBinaMeth.a $(LIBS)

test/test_quantize_roundtrip: libBinaMeth.a dmtools
	$(CC) -o $@ -I. $(CFLAGS) test/test_quantize_roundtrip.c libBinaMeth.a $(LIBS)

install: dmtools bam2dm dmDMR dmalign bam2motif

BASE_TESTS = test/testBinOrder test/testCoordinates test/test_id_roundtrip test/test_quantize_roundtrip
BIGWIG_TESTS = test/testLocal test/testRemote test/testWrite test/testRemoteManyContigs test/exampleWrite test/testIterator

ifeq ($(HAVE_BIGWIG),YES)
    ENABLED_BIGWIG_TESTS := $(BIGWIG_TESTS)
else
    ENABLED_BIGWIG_TESTS :=
endif

test: $(BASE_TESTS) dmtools bam2dm-asan $(ENABLED_BIGWIG_TESTS)
ifeq ($(HAVE_BIGWIG),YES)
	./test/test.py
else
	@echo "Skipping bigWig compatibility tests (libBigWig not available)."
endif
	./test/testBinOrder.sh
	./test/testCoordinates
	./test/test_id_roundtrip
	./test/test_quantize_roundtrip
	./test/test_sc_shrinkage_toy.py
	./test/test_sc_shrinkage_integration.py
	./test/test_dmr_bb_toy.py

clean:
	rm -f *.o libBinaMeth.a libBinaMeth.so *.pico test/testLocal test/testRemote test/testWrite dmtools dmDMR bam2dm bam2motif dmalign test/exampleWrite test/testRemoteManyContigs test/testIterator genome2cg genomebinLen dmalign
	make clean -C htslib

cleandm:
	rm dmtools

install-env: libBinaMeth.a libBinaMeth.so
	install -d $(prefix)/lib $(prefix)/include
	install libBinaMeth.a $(prefix)/lib
	install libBinaMeth.so $(prefix)/lib
	install *.h $(prefix)/include

libs:
	chmod +x htslib/version.sh 2>/dev/null || true
	$(MAKE) -C htslib libhts.a
