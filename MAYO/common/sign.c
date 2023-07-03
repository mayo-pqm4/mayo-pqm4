#include "MAYO.h"
#include "packing.h"
#include <stdlib.h>
#include <string.h>

unsigned int randomint()
{
	uint8_t buffer[4];
	unsigned int out;

	randombytes(buffer, 4);

	out = (unsigned int)buffer[0] << 24 | (unsigned int)buffer[1] << 16 | (unsigned int)buffer[2] << 8 | (unsigned int)buffer[3] << 0;

	return out;
}

void shuffle(int *array, int size)
{
	int j;

	if (size > 1)
	{
		for (int i = size - 1; i > 0; i--)
		{
			j = randomint() % (i + 1);
			int t = array[j];
			array[j] = array[i];
			array[i] = t;
		}
	}
}

int cmpfunc(const void *a, const void *b)
{
	return (*(int *)a - *(int *)b);
}

void sample_system(uint8_t *E, uint8_t *P1a, uint8_t *bilin, uint8_t *vinegar, uint8_t *linear, uint8_t *cterms)
{
	uint8_t sol[M_bytes] = {0};
	uint8_t sol1[M_bytes] = {0};
	uint8_t solLin[M_bytes * OO] = {0};
	uint8_t solLin1[M_bytes * OO] = {0};
	uint8_t xy[2 * bytes(NN - OO)] = {0};
	int countE = 0;

#if defined COMPACT_E_MATRICES
	uint8_t E_mat[E_bytes_expanded];
	shake256(E_mat, E_bytes_expanded, E, E_bytes);
#else
	uint8_t *E_mat = E;
#endif

	for (int i = 0; i < KK; i++)
	{
		mayo_eval(sol, P1a, vinegar + i * bytes(NN)); // P(vi)

		mat_prod(sol1, E_mat + countE * M_bytes * MM, sol); // Eii*P(vi)

		for (int j = 0; j < M_bytes; j++)
		{ // sum_i
			cterms[j] ^= sol1[j];
		}
		mat_prod_batch(solLin, bilin, vinegar + i * bytes(NN)); // P'(vi, oi)

		mat_x_mat_prod(solLin1, E_mat + countE * M_bytes * MM, solLin); // Eii*P'(vi, oi)

		for (int j = 0; j < bytes(MM) * OO; j++)
		{ // sum_i;
			linear[i * bytes(MM) * OO + j] ^= solLin1[j];
		}

		countE += KK - i;
	}

	countE = 0;

	for (int i = 0; i < KK; i++)
	{
		countE++;
		for (int k = 0; k < bytes(NN - OO); k++)
		{
			xy[k] = vinegar[i * bytes(NN) + k]; // vin i
		}
		for (int j = i + 1; j < KK; j++)
		{
			for (int k = 0; k < bytes(NN - OO); k++)
			{
				xy[k + bytes(NN - OO)] = vinegar[j * bytes(NN) + k]; // vin j
			}
			bilinear_eval(sol, P1a, xy); // P'(vi, vj)

			mat_prod(sol1, E_mat + countE * M_bytes * MM, sol); // Eij*P'(vi, vj)

			for (int l = 0; l < M_bytes; l++)
			{ // sum_ij
				cterms[l] ^= sol1[l];
			}
			mat_prod_batch(solLin, bilin, vinegar + i * bytes(NN)); // P'(vi,oj)

			mat_x_mat_prod(solLin1, E_mat + countE * M_bytes * MM, solLin); // Eij*P'(vi,oj)
			for (int l = 0; l < bytes(MM) * OO; l++)
			{
				linear[j * bytes(MM) * OO + l] ^= solLin1[l];
			}
			mat_prod_batch(solLin, bilin, vinegar + j * bytes(NN)); // P'(vj,oi)

			mat_x_mat_prod(solLin1, E_mat + countE * M_bytes * MM, solLin); // Eij*P'(vj,oi)
			for (int l = 0; l < bytes(MM) * OO; l++)
			{
				linear[i * bytes(MM) * OO + l] ^= solLin1[l];
			}
			countE++;
		}
	}
}

int try_gauss(uint8_t *linear, uint8_t *cterms, uint8_t *squaresol, int *ind)
{

	uint8_t square[M_bytes * MM] = {0}; // obtained square system

	shuffle(ind, KK * OO);
	qsort(ind, (KK * OO) - MM, sizeof(int), cmpfunc);
	// ordered array of (KK*OO)-MM random indeces representing the columns to ignore in the linear system

	// construct the derived square system ignoring the selected columns
	int row_sq = 0;
	int index = 0;
	for (int row_lin = 0; row_lin < KK * OO; row_lin++)
	{
		if ((row_lin == ind[index]) & (index < (KK * OO) - MM))
		{
			index += 1;
		}
		else
		{
			for (int j = 0; j < M_bytes; j++)
			{
				square[row_sq * M_bytes + j] = linear[row_lin * M_bytes + j];
			}
			row_sq += 1;
		}
	}

	int succ;
	solve_linear(succ, squaresol, square, cterms); // try to solve the square system

	if (succ == 1)
	{
		return 1;
	}
	return 0;
}

void decode_oil(uint8_t *squaresol, int *ind, uint8_t *oil_space, uint8_t *preimage)
{
	uint8_t oilsol[bytes(OO) * KK] = {0};
	int sq = 0;
	int indexx = 0;
	int fl_lin = 0;
	int fl_sq = 0;
	int byte_sq = squaresol[0];

	for (int lin = 0; lin < KK * OO; lin++)
	{
		fl_lin = lin / 2;
		if ((ind[indexx] == lin) && (indexx < (KK * OO) - MM))
		{
			indexx += 1;
		}
		else
		{
			if ((lin % 2 == 0) && (sq % 2 == 0))
			{
				oilsol[fl_lin] |= (byte_sq & 0xF);
			}
			if ((lin % 2 == 1) && (sq % 2 == 1))
			{
				oilsol[fl_lin] |= (byte_sq & 0xF0);
			}
			if ((lin % 2 == 0) && (sq % 2 == 1))
			{
				oilsol[fl_lin] |= (byte_sq >> 4);
			}
			if ((lin % 2 == 1) && (sq % 2 == 0))
			{
				oilsol[fl_lin] |= (byte_sq << 4);
			}
			sq += 1;
			fl_sq = sq / 2;
			byte_sq = squaresol[fl_sq];
		}
	}

	mat_prod_multiple_oil(preimage, oil_space, oilsol);
}

void sign(const uint8_t *msg, int msg_len, uint8_t *E, uint8_t *P1a, uint8_t *oil_space, uint8_t *bilin, uint8_t *sig)
{
	uint8_t hash[M_bytes];
	uint8_t salt[salt_bytes];
	uint8_t preimage[preimage_bytes];
	uint8_t vinegar[bytes(NN) * KK] = {0};
	uint8_t cterms[M_bytes];
	uint8_t linear[M_bytes * KK * OO];
	uint8_t squaresol[M_bytes]; // solution to the square system
	int s = 0;

	int ind[KK * OO];
	for (int i = 0; i < KK * OO; i++)
	{
		ind[i] = i;
	}

	while (s == 0)
	{
		memset(cterms, 0, M_bytes);
		memset(linear, 0, bytes(MM) * KK * OO);
		memset(squaresol, 0, M_bytes);

		randombytes(salt, salt_bytes);

		// compute hash(msg || salt)
		uint8_t conc[msg_len + salt_bytes];
		memcpy(conc, msg, msg_len);
		memcpy(conc + msg_len, salt, salt_bytes);
		shake256(hash, M_bytes, conc, msg_len + salt_bytes);

		// sample vinegar
		for (int i = 0; i < KK; i++)
		{
			randombytes(vinegar + i * bytes(NN), bytes(NN - OO));
			vinegar[i * bytes(NN) + bytes(NN - OO) - 1] &= 0xF;
		}

		sample_system(E, P1a, bilin, vinegar, linear, cterms);

		for (int j = 0; j < M_bytes; j++)
		{
			cterms[j] ^= hash[j];
		}

		s = try_gauss(linear, cterms, squaresol, ind);
	}

	// reconstruct the solution to the original linear system by adding zeroes in the ignored indeces
	decode_oil(squaresol, ind, oil_space, preimage);

	memcpy(sig, salt, salt_bytes);
	for (int i = 0; i < preimage_bytes; i++)
	{
		sig[i + salt_bytes] = preimage[i] ^ vinegar[i];
	}
}

int crypto_sign(uint8_t *sm, size_t *smlen,
				const uint8_t *m, size_t mlen,
				const uint8_t *sk)
{
	uint8_t E[E_bytes], P1a[P1a_bytes], oilspace[OilSpace_bytes], bilin[bilin_bytes];

	for (size_t i = 0; i < mlen; ++i)
		sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];

	unpack_sk(E, P1a, oilspace, bilin, sk);
	sign(m, mlen, E, P1a, oilspace, bilin, sm);
	*smlen = mlen + CRYPTO_BYTES;
	return 0;
}
