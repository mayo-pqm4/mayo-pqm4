#ifndef API_H
#define API_H

#include "params.h"
#include <stddef.h>
#include <stdint.h>

#define CRYPTO_SECRETKEYBYTES (E_bytes + P1a_bytes + OilSpace_bytes + bilin_bytes)
#define CRYPTO_PUBLICKEYBYTES (E_bytes + P_bytes)
#define CRYPTO_BYTES (bytes(NN * KK) + salt_bytes)

int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

int crypto_sign(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk);

int crypto_sign_open(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);

#endif
