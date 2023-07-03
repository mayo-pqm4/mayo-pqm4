// #ifndef MAYO_ASM_H
// #define MAYO_ASM_H
#include <stddef.h>
#include <stdint.h>

extern const uint8_t gf16mul_lut[256];

// sign and verify
extern void mat_prod_32_64_asm(uint8_t *result, const uint8_t *rmat, uint8_t *vec);

#define mat_prod(result, rmat, vec) mat_prod_32_64_asm(result, rmat, vec)

// sign
extern void mayo_eval_59_32_asm(uint8_t *result, const uint8_t *tmat, uint8_t *vec);
extern void bilinear_eval_59_32(uint8_t *result, const uint8_t *tmat, uint8_t *vec);
extern void mat_prod_batch_7_59_32(uint8_t *result, const uint8_t *trimat, uint8_t *vec);
extern void mat_x_mat_prod_32_64_7(uint8_t *result, const uint8_t *rmat, uint8_t *rmat2);
extern void mat_prod_multiple_oil_7_7x59_14(uint8_t *result, const uint8_t *rmat, uint8_t *vec);
extern unsigned int solve_linear_64(uint8_t *sol, uint8_t *inp_mat, uint8_t *c_terms);

#define mayo_eval(result, tmat, vec) mayo_eval_59_32_asm(result, tmat, vec)
#define bilinear_eval(result, tmat, vec) bilinear_eval_59_32(result, tmat, vec)
#define mat_prod_batch(result, trimat, vec) mat_prod_batch_7_59_32(result, trimat, vec)
#define mat_x_mat_prod(result, rmat, rmat2) mat_x_mat_prod_32_64_7(result, rmat, rmat2)
#define mat_prod_multiple_oil(result, rmat, vec) mat_prod_multiple_oil_7_7x59_14(result, rmat, vec)
#define solve_linear(solved, sol, inp_mat, c_terms) (solved = solve_linear_64(sol, inp_mat, c_terms))

// verify
extern void mayo_eval_66_32_leaktime_asm(uint8_t *result, const uint8_t *trimat, const uint8_t *x, const uint8_t *lut);
extern void bilinear_eval_66_32_leaktime_asm(uint8_t *y, const uint8_t *trimat, const uint8_t *xx, const uint8_t *lut);

#define mayo_eval_leaktime(result, trimat, x, lut) mayo_eval_66_32_leaktime_asm(result, trimat, x, lut)
#define bilinear_eval_leaktime(y, trimat, xx, lut) bilinear_eval_66_32_leaktime_asm(y, trimat, xx, lut)