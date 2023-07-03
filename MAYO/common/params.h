#ifndef PARAMS_H
#define PARAMS_H

#include "config.h"

#define tri_terms(n_var) ((n_var) * ((n_var) + 1) / 2)
#define bytes(n_var) (((n_var) + ((n_var) % 2)) / 2)

#if defined MAYO_ORIG_I
#define MM 67
#define NN 66
#define OO 5
#define KK 14
#elif defined MAYO_NEW_I
#define MM 64
#define NN 66
#define OO 7
#define KK 10
#elif defined MAYO_NEW_I_K14
#define MM 64
#define NN 66
#define OO 7
#define KK 14
#endif

#define sk_seed_bytes 32
#define pk_seed_bytes 32
#define oil_seed_bytes 32

#define M_bytes (bytes(MM))

#define OilSpace_bytes (bytes(NN - OO) * OO)
#define P_bytes (tri_terms(NN) * M_bytes)
#define P1a_bytes (tri_terms(NN - OO) * M_bytes)
#define bilin_bytes ((NN - OO) * OO * M_bytes)

#if defined COMPACT_E_MATRICES
#define E_bytes 32
#define E_bytes_expanded (M_bytes * MM * tri_terms(KK))
#else
#define E_bytes (M_bytes * MM * tri_terms(KK))
#endif

#define preimage_bytes (bytes(NN) * KK)
#define salt_bytes 32

#endif
