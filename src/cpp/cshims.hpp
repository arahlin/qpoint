#pragma once

// Single point of contact with the bundled C libraries. These stay as C and
// are linked, not rewritten: ERFA and the HEALPix C library are third-party,
// and sincos.c's polynomial approximations must stay bit-identical or every
// fast_math=True result diverges from the C build.
extern "C" {
#include "chealpix.h"
#include "erfa.h"
#include "fast_math.h"
}
