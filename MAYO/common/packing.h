#ifndef PACKING_H
#define PACKING_H

#include <stdint.h>

void unpack_sk(uint8_t *E, uint8_t *P1a, uint8_t *oilspace, uint8_t *bilin, const uint8_t *sk);
void unpack_pk(uint8_t *E, uint8_t *P, const uint8_t *pk);
void unpack_sig(uint8_t *salt, uint8_t *sig, const uint8_t *sm);

#endif