#include "MAYO.h"
#include "gf16.h"
#include <string.h>

/*************************************************
 * Name:        randomize_matrix
 *
 * Description: create a random matrix with uint8_t entries in [0,...,15] from a seed
 *
 * Arguments:   - int nRows: number of rows
 *              - int nCols: number of cols
 *              - uint8_t matrix[nRows][nCols]: matrix to be randomized
 *              - uint8_t *seed: seed to be expanded
 *              - size_t seed_len: seed length
 **************************************************/
void randomize_matrix(int nRows, int nCols, uint8_t matrix[nRows][nCols], uint8_t *seed, size_t seed_len)
{
    shake256incctx ctx;
    shake256_inc_init(&ctx);
    shake256_inc_absorb(&ctx, seed, seed_len);
    shake256_inc_finalize(&ctx);

    for (int i = 0; i < nRows; i++)
    {
        shake256_inc_squeeze(matrix[i], nCols, &ctx);
        for (int j = 0; j < nCols; j++)
        {
            matrix[i][j] &= 0xF;
        }
    }

    shake256_inc_ctx_release(&ctx);
}

/*************************************************
 * Name:        create_UTrandom_matrix
 *
 * Description: create an upper triangular random matrix with uint8_t entries in [0,...,15] from a seed
 *
 * Arguments:   - int nRows: number of rows
 *              - int nCols: number of cols
 *              - uint8_t matrix[nRows][nCols]: matrix to be randomized
 *              - uint8_t *seed: seed to be expanded
 *              - size_t seed_len: seed length
 **************************************************/
void randomize_ut_matrix(int nRows, int nCols, uint8_t matrix[nRows][nCols], uint8_t *seed, size_t seed_len)
{
    shake256incctx ctx;
    shake256_inc_init(&ctx);
    shake256_inc_absorb(&ctx, seed, seed_len);
    shake256_inc_finalize(&ctx);

    for (int i = 0; i < nRows; i++)
    {
        shake256_inc_squeeze(matrix[i], nCols, &ctx);
        for (int j = 0; j < nCols; j++)
        {
            if (i <= j)
            {
                matrix[i][j] &= 0xF;
            }
            else
            {
                matrix[i][j] = 0;
            }
        }
    }

    shake256_inc_ctx_release(&ctx);
}

/*************************************************
 * Name:        mat_mul_mat
 *
 * Description: multiply two matrices
 *
 * Arguments:   - int nRowsA: number of rows of A
 *              - int nColsA: number of cols of A
 *              - int nColsB: number of cols of B
 *              - uint8_t A[nRowsA][nColsA]: first matrix to multiply
 *              - uint8_t B[nColsA][nColsB]: second matrix to multiply
 *              - uint8_t result[nRowsA][nColsB]: result of the multiplication
 **************************************************/
void mat_mul_mat(int nRowsA, int nColsA, int nColsB, uint8_t A[nRowsA][nColsA], uint8_t B[nColsA][nColsB], uint8_t result[nRowsA][nColsB])
{
    for (int i = 0; i < nRowsA; i++)
    {
        for (int j = 0; j < nColsB; j++)
        {
            result[i][j] = 0;
            for (int l = 0; l < nColsA; l++)
            {
                result[i][j] ^= gf16_mul(A[i][l], B[l][j]);
            }
        }
    }
}

/*************************************************
 * Name:        mat_mul_tran
 *
 * Description: multiply the first matrix with the transpose of the second matrix
 *
 * Arguments:   - int nRowsA: number of rows of A
 *              - int nColsA: number of cols of A
 *              - int nColsB: number of cols of B
 *              - uint8_t A[nRowsA][nColsA]: first matrix to multiply
 *              - uint8_t B[nRowsB][nColsA]: second matrix to multiply
 *              - uint8_t result[nRowsA][nRowsB]: result of the multiplication
 **************************************************/
void mat_mul_tran(int nRowsA, int nColsA, int nRowsB, uint8_t A[nRowsA][nColsA], uint8_t B[nRowsB][nColsA], uint8_t result[nRowsA][nRowsB])
{
    for (int i = 0; i < nRowsA; i++)
    {
        for (int j = 0; j < nRowsB; j++)
        {
            result[i][j] = 0;
            for (int l = 0; l < nColsA; l++)
            {
                result[i][j] ^= gf16_mul(A[i][l], B[j][l]);
            }
        }
    }
}

/*************************************************
 * Name:        mat_add_tran_mul_tran
 *
 * Description: multiply the sum of the first matrix with its transpose by the transpose of the second matrix
 *
 * Arguments:   - int nRowsA: number of rows of A
 *              - int nColsA: number of cols of A
 *              - int nColsB: number of cols of B
 *              - uint8_t A[nRowsA][nColsA]: first matrix to multiply
 *              - uint8_t B[nRowsB][nColsA]: second matrix to multiply
 *              - uint8_t result[nRowsA][nRowsB]: result of the multiplication
 **************************************************/
void mat_add_tran_mul_tran(int nRowsA, int nColsA, int nRowsB, uint8_t A[nRowsA][nColsA], uint8_t B[nRowsB][nColsA], uint8_t result[nRowsA][nRowsB])
{
    for (int i = 0; i < nRowsA; i++)
    {
        for (int j = 0; j < nRowsB; j++)
        {
            result[i][j] = 0;
            for (int l = 0; l < nColsA; l++)
            {
                result[i][j] = gf16_mul(A[i][l] ^ A[l][i], B[j][l]) ^ result[i][j];
            }
        }
    }
}

/*************************************************
 * Name:        upper_mat
 *
 * Description: make a target matrix upper triangular
 *
 * Arguments:   - int nRows: number of rows
 *              - int nCols: number of columns
 *              - uint8_t matrix[nRows][nCols]: matrix to convert to UT
 **************************************************/
void upper_mat(int nRows, int nCols, uint8_t matrix[nRows][nCols])
{
    for (int i = 0; i < nRows; i++)
    {
        for (int j = i + 1; j < nCols; j++)
        {
            matrix[i][j] = (matrix[i][j] ^ matrix[j][i]);
            matrix[j][i] = 0;
        }
    }
}

/*************************************************
 * Name:        gf16matrix_to_bytes
 *
 * Description: convert a matrix with entries in GF16 to an array of bytes
 *
 * Arguments:   - int nRows: number of rows
 *              - int nCols: number of columns
 *              - uint8_t matrix[nRows][nCols]: input matrix
 *              - uint8_t *result: output array
 **************************************************/
void gf16matrix_to_bytes(int nRows, int nCols, uint8_t matrix[nRows][nCols], uint8_t *result)
{
    int index = 0;
    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < nCols - (nCols % 2); j += 2)
        {
            result[index] = matrix[i][j] + 16 * matrix[i][j + 1];
            index++;
        }

        if ((nCols % 2) == 1)
        {
            result[index] = matrix[i][nCols - 1];
            index++;
        }
    }
}

/*************************************************
 * Name:        ut_gf16matrix_list_to_bytes
 *
 * Description: convert a list of square UT matrices with entries in GF16 to an array of bytes
 *
 * Arguments:   - int mDim: dimensions of each matrix
 *              - int mLen: number of matrices
 *              - uint8_t matrix_list[mLen][mDim][mDim]: input matrices list
 *              - uint8_t *result: output array
 **************************************************/
void ut_gf16matrix_list_to_bytes(int mDim, int mLen, uint8_t matrix_list[mLen][mDim][mDim], uint8_t *result)
{
    int index = 0;

    for (int i = 0; i < mDim; ++i)
    {
        for (int j = i; j < mDim; ++j)
        {
            for (int k = 0; k < mLen - (mLen % 2); k += 2)
            {
                result[index] = matrix_list[k][i][j] + 16 * matrix_list[k + 1][i][j];
                index++;
            }

            if ((mLen % 2) == 1)
            {
                result[index] = matrix_list[mLen - 1][i][j];
                index++;
            }
        }
    }
}

/*************************************************
 * Name:        gf16matrix_list_to_bytes
 *
 * Description: convert a list of matrices with entries in GF16 to an array of bytes
 *
 * Arguments:   - int nRows: number of rows of each matrix
 *              - int nCols: number of cols of each matrix
 *              - int mLen: number of matrices
 *              - uint8_t matrix_list[mLen][nRows][nCols]: input matrices list
 *              - uint8_t *result: output array
 **************************************************/
void gf16matrix_list_to_bytes(int nRows, int nCols, int mLen, uint8_t matrix_list[mLen][nRows][nCols], uint8_t *result)
{
    int index = 0;

    for (int i = 0; i < nRows; ++i)
    {
        for (int j = 0; j < nCols; ++j)
        {
            for (int k = 0; k < mLen - (mLen % 2); k += 2)
            {
                result[index] = matrix_list[k][i][j] + 16 * matrix_list[k + 1][i][j];
                index++;
            }

            if ((mLen % 2) == 1)
            {
                result[index] = matrix_list[mLen - 1][i][j];
                index++;
            }
        }
    }
}

/*************************************************
 * Name:        compute_P2
 *
 * Description: compute a list of matrices given by P2[i] = -O*P1a[i]*O^t - P1b[i]*O^t
 *
 * Arguments:   - uint8_t O_matrix[OO][NN-OO]: oil space matrix
 *              - uint8_t P1a[MM][NN-OO][NN-OO]: list of public P1a matrices
 *              - uint8_t P1b[MM][NN-OO][OO]: list of public P1b matrices
 *              - uint8_t P2[MM][OO][OO]: list of resulting P2 matrices
 **************************************************/
void compute_P2(uint8_t O_matrix[OO][NN - OO], uint8_t P1a[MM][NN - OO][NN - OO], uint8_t P1b[MM][NN - OO][OO], uint8_t P2[MM][OO][OO])
{
    uint8_t tmpOne[OO][NN - OO];
    uint8_t tmpTwo[OO][OO];
    uint8_t tmpThree[OO][OO];

    for (int i = 0; i < MM; ++i)
    {
        mat_mul_mat(OO, NN - OO, NN - OO, O_matrix, P1a[i], tmpOne);
        mat_mul_tran(OO, NN - OO, OO, tmpOne, O_matrix, tmpTwo);
        mat_mul_mat(OO, NN - OO, OO, O_matrix, P1b[i], tmpThree);
        for (int j = 0; j < OO; ++j)
        {
            for (int k = 0; k < OO; ++k)
            {
                (P2[i])[j][k] = tmpTwo[j][k] ^ tmpThree[j][k];
            }
        }
    }

    for (int i = 0; i < MM; i++)
    {
        upper_mat(OO, OO, P2[i]);
    }
}

/*************************************************
 * Name:        compute_bilin
 *
 * Description: compute a list of matrices given by bilin[i] = (P1a[i] + P1a[i]^t)*O^t + P1b[i]
 *
 * Arguments:   - uint8_t O_matrix[OO][NN-OO]: oil space matrix
 *              - uint8_t P1a[MM][NN-OO][NN-OO]: list of public P1a matrices
 *              - uint8_t P1b[MM][NN-OO][OO]: list of public P1b matrices
 *              - uint8_t bilin[MM][NN-OO][OO]: list of resulting bilin matrices
 **************************************************/
void compute_bilin(uint8_t O_matrix[OO][NN - OO], uint8_t P1a[MM][NN - OO][NN - OO], uint8_t P1b[MM][NN - OO][OO], uint8_t bilin[MM][NN - OO][OO])
{
    uint8_t prod[NN - OO][OO];

    for (int i = 0; i < MM; ++i)
    {
        mat_add_tran_mul_tran(NN - OO, NN - OO, OO, P1a[i], O_matrix, prod);

        for (int j = 0; j < NN - OO; ++j)
        {
            for (int k = 0; k < OO; ++k)
            {
                (bilin[i])[j][k] = (prod[j][k] ^ (P1b[i])[j][k]);
            }
        }
    }
}

/*************************************************
 * Name:        assemble_P_to_bytes
 *
 * Description: assemble a list of block matrices of the form
 *              [[P1a[i], P1b[i]], [0, P2[i]]] and convert to an
 *              array of bytes
 *
 * Arguments:   - uint8_t P1a[MM][NN-OO][NN-OO]: list of public P1a matrices
 *              - uint8_t P1b[MM][NN-OO][OO]: list of public P1b matrices
 *              - uint8_t P2[MM][OO][OO]: list of public P2 matrices
 *              - uint8_t *result: output array
 **************************************************/
void assemble_P_to_bytes(uint8_t P1a[MM][NN - OO][NN - OO], uint8_t P1b[MM][NN - OO][OO], uint8_t P2[MM][OO][OO], uint8_t *result)
{
    int index = 0;

    for (int i = 0; i < NN; ++i)
    {
        for (int j = i; j < NN; ++j)
        {
            for (int k = 0; k < MM - (MM % 2); k += 2)
            {
                if (i < NN - OO && j < NN - OO)
                {
                    result[index] = P1a[k][i][j] + 16 * P1a[k + 1][i][j];
                }
                else if (i < NN - OO)
                {
                    result[index] = P1b[k][i][j - (NN - OO)] + 16 * P1b[k + 1][i][j - (NN - OO)];
                }
                else if (j < NN - OO)
                {
                    result[index] = 0;
                }
                else
                {
                    result[index] = P2[k][i - (NN - OO)][j - (NN - OO)] + 16 * P2[k + 1][i - (NN - OO)][j - (NN - OO)];
                }

                index++;
            }

            if ((MM % 2) == 1)
            {
                if (i < NN - OO && j < NN - OO)
                {
                    result[index] = P1a[MM - 1][i][j];
                }
                else if (i < NN - OO)
                {
                    result[index] = P1b[MM - 1][i][j - (NN - OO)];
                }
                else if (j < NN - OO)
                {
                    result[index] = 0;
                }
                else
                {
                    result[index] = P2[MM - 1][i - (NN - OO)][j - (NN - OO)];
                }

                index++;
            }
        }
    }
}

/*************************************************
 * Name:        derive_sk
 *
 * Description: derive secret key from a the public key seed and the
 *              oil space seed
 *
 * Arguments:   - uint8_t *sk: secret key
 *              - uint8_t *pk_seed: public key seed
 *              - uint8_t *oil_seed: oil space seed
 **************************************************/
void derive_sk(uint8_t *sk, uint8_t *pk_seed, uint8_t *oil_seed)
{
    uint8_t P_seed[2 * MM * pk_seed_bytes];
    uint8_t O_matrix[OO][NN - OO];
    uint8_t P1a[MM][NN - OO][NN - OO];
    uint8_t bilin[MM][NN - OO][OO];

    shake256(P_seed, 2 * MM * pk_seed_bytes, pk_seed, pk_seed_bytes);

    {
        randomize_matrix(OO, NN - OO, O_matrix, oil_seed, oil_seed_bytes);

        for (int i = 0; i < MM; i++)
        {
            randomize_ut_matrix(NN - OO, NN - OO, P1a[i], P_seed + i * pk_seed_bytes, pk_seed_bytes);
        }

        uint8_t P1b[MM][NN - OO][OO];
        for (int i = 0; i < MM; i++)
        {
            randomize_matrix(NN - OO, OO, P1b[i], P_seed + (MM + i) * pk_seed_bytes, pk_seed_bytes);
        }

        /* bilin = (P1a + P1a^T)*O^T + P1b */
        compute_bilin(O_matrix, P1a, P1b, bilin);
    }

    shake256(sk, E_bytes, P_seed, pk_seed_bytes);
    ut_gf16matrix_list_to_bytes(NN - OO, MM, P1a, sk + E_bytes);
    gf16matrix_to_bytes(OO, NN - OO, O_matrix, sk + E_bytes + P1a_bytes);
    gf16matrix_list_to_bytes(NN - OO, OO, MM, bilin, sk + E_bytes + P1a_bytes + OilSpace_bytes);
}

/*************************************************
 * Name:        derive_pk
 *
 * Description: derive public key from a the public key seed and the
 *              oil space seed
 *
 * Arguments:   - uint8_t *pk: public key
 *              - uint8_t *pk_seed: public key seed
 *              - uint8_t *oil_seed: oil space seed
 **************************************************/
void derive_pk(uint8_t *pk, uint8_t *pk_seed, uint8_t *oil_seed)
{
    uint8_t P_seed[2 * MM * pk_seed_bytes];
    uint8_t P2[MM][OO][OO];

    shake256(P_seed, 2 * MM * pk_seed_bytes, pk_seed, pk_seed_bytes);

    uint8_t O_matrix[OO][NN - OO];
    randomize_matrix(OO, NN - OO, O_matrix, oil_seed, oil_seed_bytes);

    uint8_t P1a[MM][NN - OO][NN - OO] = {0};
    for (int i = 0; i < MM; i++)
    {
        randomize_ut_matrix(NN - OO, NN - OO, P1a[i], P_seed + i * pk_seed_bytes, pk_seed_bytes);
    }

    uint8_t P1b[MM][NN - OO][OO];
    for (int i = 0; i < MM; i++)
    {
        randomize_matrix(NN - OO, OO, P1b[i], P_seed + (MM + i) * pk_seed_bytes, pk_seed_bytes);
    }

    compute_P2(O_matrix, P1a, P1b, P2);

    shake256(pk, E_bytes, P_seed, pk_seed_bytes);
    assemble_P_to_bytes(P1a, P1b, P2, pk + E_bytes);
}

/*************************************************
 * Name:        crypto_sign_keypair
 *
 * Description: Generates public and private key.
 *
 * Arguments:   - uint8_t *pk: pointer to output public key
 *              - uint8_t *sk: pointer to output private key
 *
 * Returns 0 (success)
 **************************************************/
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk)
{
    uint8_t sk_seed[sk_seed_bytes];
    uint8_t exp_seed[pk_seed_bytes + oil_seed_bytes];
    uint8_t *pk_seed;
    uint8_t *oil_seed;

    randombytes(sk_seed, sk_seed_bytes);

    shake256(exp_seed, pk_seed_bytes + oil_seed_bytes, sk_seed, sk_seed_bytes);
    pk_seed = exp_seed;
    oil_seed = exp_seed + pk_seed_bytes;

    derive_sk(sk, pk_seed, oil_seed);
    derive_pk(pk, pk_seed, oil_seed);

    return 0;
}
