#ifndef _GF16_H_
#define _GF16_H_

#include <stdint.h>

// gf16 := gf2[x]/(x^4+x+1)
static inline uint8_t gf16_mul(uint8_t a, uint8_t b)
{
    uint8_t r8 = (a & 1) * b;
    r8 ^= (a & 2) * b;
    r8 ^= (a & 4) * b;
    r8 ^= (a & 8) * b;

    // reduction
    uint8_t r4 = r8 ^ (((r8 >> 4) & 5) * 3); // x^4 = x+1  , x^6 = x^3 + x^2
    r4 ^= (((r8 >> 5) & 1) * 6);             // x^5 = x^2 + x
    return (r4 & 0xf);
}

#endif // _GF16_H_
