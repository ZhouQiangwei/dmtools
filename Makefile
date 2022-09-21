#ifndef CPP
#$(error CPP variable undefined)
#endif

CPP = $(shell pwd)
SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
PROGS = bam2bm bmtools bmDMR

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
tmpfile:=$(shell mktemp --suffix=.c)
$(file >$(tmpfile),#include <curl/curl.h>)
$(file >>$(tmpfile),int main() { return 0; })
HAVE_CURL:=$(shell $(CC) $(CFLAGS) $(EXTRA_CFLAGS_PIC) $(LIBS) -lcurl $(tmpfile) -o /dev/null >/dev/null 2>&1 && echo "YES")
$(shell rm $(tmpfile))

ifeq ($(HAVE_CURL),YES)
	# If yes, add the library
	LIBS += -lcurl
else
	# and if not, disable CURL specific code compilation
	CFLAGS += -DNOCURL
endif


prefix = /usr/local
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib

.PHONY: all clean lib test doc

.SUFFIXES: .c .o .pico

all: libs lib $(PROGS)

lib: lib-static lib-shared

lib-static: libBigWig.a

lib-shared: libBigWig.so

doc:
	doxygen

OBJS = io.o bwValues.o bwRead.o bwStats.o bwWrite.o

.c.o:
	$(CC) -I. $(CFLAGS) $(INCLUDES) -c -o $@ $<

.c.pico:
	$(CC) -I. $(CFLAGS) $(INCLUDES) $(EXTRA_CFLAGS_PIC) -c -o $@ $<

libBigWig.a: $(OBJS)
	-@rm -f $@
	$(AR) -rcs $@ $(OBJS)
	$(RANLIB) $@

libBigWig.so: $(OBJS:.o=.pico)
	$(CC) -shared -o $@ $(OBJS:.o=.pico) $(LDLIBS) $(LIBS)
#$(CC) -shared $(LDFLAGS) -o $@ $(OBJS:.o=.pico) $(LDLIBS) $(LIBS)

test/testLocal: libBigWig.a
	$(CC) -o $@ -I. $(CFLAGS) test/testLocal.c libBigWig.a $(LIBS)

test/testRemoteManyContigs: libBigWig.a
	$(CC) -o $@ -I. $(CFLAGS) test/testRemoteManyContigs.c libBigWig.a $(LIBS)

test/testRemote: libBigWig.a
	$(CC) -o $@ -I. $(CFLAGS) test/testRemote.c libBigWig.a $(LIBS)

test/testWrite: libBigWig.a
	$(CC) -o $@ -I. $(CFLAGS) test/testWrite.c libBigWig.a $(LIBS)

bmtools: libBigWig.so
	$(CC) -o $@ -I. -L. $(CFLAGS) bmtools.c -lBigWig $(LIBS) -Wl,-rpath $(CPP) -lpthread

bmDMR:
	g++ $(CFLAGS) -c -o regression.o regression.cpp
	g++ $(CFLAGS) -o bmDMR bmDMR.cpp regression.o -I. -L. -lBigWig -Wl,-rpath $(CPP) -lgsl -lgslcblas -lm -lz

bam2bm:
	g++ $(CFLAGS) bam2bm.cpp -o bam2bm -m64 -I. -L. -lz -lBigWig -Wl,-rpath $(CPP) $(LDFLAGS_SUB)

test/exampleWrite: libBigWig.so
	$(CC) -o $@ -I. -L. $(CFLAGS) test/exampleWrite.c -lBigWig $(LIBS) -Wl,-rpath .

test/testBigBed: libBigWig.a
	$(CC) -o $@ -I. $(CFLAGS) test/testBigBed.c libBigWig.a $(LIBS)

test/testIterator: libBigWig.a
	$(CC) -o $@ -I. $(CFLAGS) test/testIterator.c libBigWig.a $(LIBS)

install: bam2bm bmtools bmDMR

test: test/testLocal test/testRemote test/testWrite test/testLocal bmtools test/exampleWrite test/testRemoteManyContigs test/testBigBed test/testIterator
	./test/test.py

clean:
	rm -f *.o libBigWig.a libBigWig.so *.pico test/testLocal test/testRemote test/testWrite bmtools bmDMR bam2bm test/exampleWrite test/testRemoteManyContigs test/testBigBed test/testIterator example_output.bw
	make clean -C htslib

install-old: libBigWig.a libBigWig.so
	install -d $(prefix)/lib $(prefix)/include
	install libBigWig.a $(prefix)/lib
	install libBigWig.so $(prefix)/lib
	install *.h $(prefix)/include

libs:
	make -C htslib
