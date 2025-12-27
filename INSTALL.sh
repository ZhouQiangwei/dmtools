git clone --recursive https://github.com/samtools/htslib.git
make -C htslib
make

# Optional (for sc-export --to h5ad):
# python3 -m pip install anndata scipy pandas

### error, undefined SIZE_MAX
#htslib/htslib/kstring.h
#include <stdint.h>
#ifndef SIZE_MAX
#define SIZE_MAX (4294967295U)
#endif
