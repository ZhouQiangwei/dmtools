#ifndef CPP
#$(error CPP variable undefined)
#endif

RPATH = $(shell pwd)
SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
PROGS = dmtools bam2dm dmDMR dmalign

CXX = g++
CC = gcc
AR = ar
RANLIB = ranlib
CFLAGS = -g -w -O3 -Wsign-compare
#-fopenmp
#-Wall
## changed: = instaed of ?= above 4 lines
## -fopenmp for multi-threads
LIBS = -lm -lz
EXTRA_CFLAGS_PIC = -fpic
LDFLAGS_SUB = htslib/libhts.a -lz -lpthread -llzma -lbz2 -lcurl
LDLIBS =
INCLUDES = 

# Create a simple test-program to check if gcc can compile with curl
#tmpfile:=$(shell mktemp --suffix=.c)
tmpfile:=$(shell mktemp -t temphaoqiaoXXXXX.c)
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
	$(CC) -o $@ -I. $(CFLAGS) test/testLocal.c libBinaMeth.a $(LIBS)

test/testRemoteManyContigs: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/testRemoteManyContigs.c libBinaMeth.a $(LIBS)

test/testRemote: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/testRemote.c libBinaMeth.a $(LIBS)

test/testWrite: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/testWrite.c libBinaMeth.a $(LIBS)

dmtools: libBinaMeth.so
	$(CC) -o $@ -I. -L. $(CFLAGS) dmtools.c -lBinaMeth $(LIBS) -Wl,-rpath $(RPATH) -lpthread

#bam2dm: libBinaMeth.so
#	$(CXX) -o $@ -I. -L. $(CFLAGS) bam2dm.cpp -lBinaMeth -Wl,-rpath $(RPATH) htslib/libhts.a -llzma -lbz2 -lz

dmDMR:
	$(CXX) $(CFLAGS) -c -o regression.o regression.cpp -lgsl -lgslcblas -lm -lz
	$(CXX) $(CFLAGS) -o dmDMR dmDMR.cpp regression.o -I. -L. -lBinaMeth -Wl,-rpath $(RPATH) -lgsl -lgslcblas -lm -lz

dmalign:
	$(CXX) $(CFLAGS) -o genome2cg genome2cg.cpp
	$(CXX) $(CFLAGS) -o genomebinLen genomebinLen.cpp
	$(CXX) $(CFLAGS) dmalign.cpp -o dmalign -lz

bam2dm:
	$(CXX) -std=c++11 $(CFLAGS) bam2dm.cpp -o bam2dm -m64 -I. -L. -lz -lBinaMeth -Wl,-rpath $(RPATH) $(LDFLAGS_SUB)

test/exampleWrite: libBinaMeth.so
	$(CC) -o $@ -I. -L. $(CFLAGS) test/exampleWrite.c -lBinaMeth $(LIBS) -Wl,-rpath .

test/testIterator: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/testIterator.c libBinaMeth.a $(LIBS)

install: dmtools bam2dm dmDMR dmalign

test: test/testLocal test/testRemote test/testWrite test/testLocal dmtools test/exampleWrite test/testRemoteManyContigs test/testIterator
	./test/test.py

clean:
	rm -f *.o libBinaMeth.a libBinaMeth.so *.pico test/testLocal test/testRemote test/testWrite dmtools dmDMR bam2dm dmalign test/exampleWrite test/testRemoteManyContigs test/testIterator genome2cg genomebinLen dmalign
	make clean -C htslib

cleandm:
	rm dmtools

install-env: libBinaMeth.a libBinaMeth.so
	install -d $(prefix)/lib $(prefix)/include
	install libBinaMeth.a $(prefix)/lib
	install libBinaMeth.so $(prefix)/lib
	install *.h $(prefix)/include

libs:
	make -C htslib
