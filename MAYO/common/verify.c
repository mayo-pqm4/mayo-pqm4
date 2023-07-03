#include "MAYO.h"
#include "packing.h"
#include <string.h>

int verify(const uint8_t *msg, int msg_len, uint8_t *salt, uint8_t *preimage, uint8_t *E, uint8_t *P)
{

	uint8_t hash[M_bytes];
	uint8_t result[M_bytes] = {0};
	uint8_t xy[NN];
	uint8_t sol[M_bytes];
	uint8_t sol1[M_bytes];
	int countE = 0;

#if defined COMPACT_E_MATRICES
	uint8_t E_mat[E_bytes_expanded];
	shake256(E_mat, E_bytes_expanded, E, E_bytes);
#else
	uint8_t *E_mat = E;
#endif

	// compute hash(msg || salt)
	uint8_t conc[msg_len + salt_bytes];
	memcpy(conc, msg, msg_len);
	memcpy(conc + msg_len, salt, salt_bytes);
	shake256(hash, M_bytes, conc, msg_len + salt_bytes);

	for (int i = 0; i < KK; i++)
	{
		mayo_eval_leaktime(sol, P, preimage + i * bytes(NN), gf16mul_lut);

		mat_prod(sol1, E_mat + countE * M_bytes * MM, sol);
		for (int j = 0; j < M_bytes; j++)
		{
			result[j] ^= sol1[j];
		}
		countE += KK - i;
	}

	countE = 0;
	for (int i = 0; i < KK; i++)
	{
		countE++;
		for (int k = 0; k < bytes(NN); k++)
		{
			xy[k] = preimage[i * bytes(NN) + k];
		}
		for (int j = i + 1; j < KK; j++)
		{

			for (int k = 0; k < bytes(NN); k++)
			{
				xy[k + bytes(NN)] = preimage[j * bytes(NN) + k];
			}
			bilinear_eval_leaktime(sol, P, xy, gf16mul_lut);

			mat_prod(sol1, E_mat + countE * M_bytes * MM, sol);

			for (int k = 0; k < M_bytes; k++)
			{
				result[k] ^= sol1[k];
			}
			countE++;
		}
	}

	// Ignore last nibble if MM is odd
	if (MM % 2 != 0)
	{
		hash[M_bytes - 1] &= 0xF;
		result[M_bytes - 1] &= 0xF;
	}

	if (memcmp(result, hash, M_bytes) != 0)
	{
		return -1;
	}

	return 0;
}

int crypto_sign_open(uint8_t *m, size_t *mlen,
					 const uint8_t *sm, size_t smlen,
					 const uint8_t *pk)
{
	uint8_t E[E_bytes], salt[salt_bytes], sig[preimage_bytes], P[P_bytes];

	if (smlen < CRYPTO_BYTES)
		goto badsig;

	*mlen = smlen - CRYPTO_BYTES;

	unpack_pk(E, P, pk);
	unpack_sig(salt, sig, sm);

	if (verify(sm + CRYPTO_BYTES, *mlen, salt, sig, E, P) != 0)
	{
		goto badsig;
	}
	else
	{
		for (size_t i = 0; i < *mlen; ++i)
		{
			m[i] = sm[CRYPTO_BYTES + i];
		}
		return 0;
	}

badsig:
	*mlen = -1;
	for (size_t i = 0; i < smlen; ++i)
		m[i] = 0;

	return -1;
}
