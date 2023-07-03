#include "MAYO.h"
#include <string.h>

void unpack_sk(uint8_t *E, uint8_t *P1a, uint8_t *oilspace, uint8_t *bilin, const uint8_t *sk)
{
	memcpy(E, sk, E_bytes);
	memcpy(P1a, sk + E_bytes, P1a_bytes);
	memcpy(oilspace, sk + E_bytes + P1a_bytes, OilSpace_bytes);
	memcpy(bilin, sk + E_bytes + P1a_bytes + OilSpace_bytes, bilin_bytes);
}

void unpack_pk(uint8_t *E, uint8_t *P, const uint8_t *pk)
{
	memcpy(E, pk, E_bytes);
	memcpy(P, pk + E_bytes, P_bytes);
}

void unpack_sig(uint8_t *salt, uint8_t *sig, const uint8_t *sm)
{
	memcpy(salt, sm, salt_bytes);
	memcpy(sig, sm + salt_bytes, preimage_bytes);
}