git clone --recursive https://github.com/samtools/htslib.git
make -C htslib
make

### error, undefined SIZE_MAX
#htslib/htslib/kstring.h
#include <stdint.h>
#ifndef SIZE_MAX
#define SIZE_MAX (4294967295U)
#endif
