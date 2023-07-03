.syntax unified
.cpu cortex-m4
.thumb

#include "bitslice.inc"
#include "madd_bitsliced.inc"

//extern void mat_prod_16_32_asm(uint32_t *c, uint32_t *a, uint32_t *b);
.global mat_prod_32_64_asm
.type mat_prod_32_64_asm, %function
.align 2
mat_prod_32_64_asm:
    push {r4-r11, r14}
    result .req r0                          
    rmat .req r1                            
    vec .req r2                             
    mat0    .req r0
    mat1    .req r2
    mat2    .req r10
    mat3    .req r12
    buf_vec     .req r5
    tmp0    .req r6
    tmp1    .req r7
    tmp2    .req r8
    tmp3    .req r9
    accu0    .req r3
    accu1    .req r11
    accu2    .req r4
    accu3    .req r14

    one      .req s0
    ctr1     .req s1
    ctr2     .req s2

    push.w {result}
    push.w {vec}
    vmov one, #0.5
    mov.w accu0, #0
    mov.w accu1, #0
    mov.w accu2, #0
    mov.w accu3, #0
    vmov ctr1, #4                               
    
    1:
        pop.w {vec}
        ldr.w buf_vec, [vec], #4                
        push.w {vec}
        .set kk, 0
        .rept 8                                 
            ldr.w tmp1, [rmat, #4]             
            ldr.w tmp2, [rmat, #8]
            ldr.w tmp3, [rmat, #12]
            ldr.w tmp0, [rmat], #32
            // bitslice on the fly
            bitslice mat0, mat1, mat2, mat3, tmp0, tmp1, tmp2, tmp3
            madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, buf_vec, tmp0, tmp1, tmp2, tmp3
            .if kk != 7
            lsr.w buf_vec, buf_vec, #4
            .endif
            .set kk, kk+1
        .endr
    vsub.f32 ctr1, ctr1, one
    vcmp.f32 ctr1, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 1b

    pop.w {vec}
    pop.w {result}
    // un-bitslice on the fly
    bitslice tmp0, tmp1, tmp2, tmp3, accu0, accu1, accu2, accu3
    str.w tmp0, [result]
    str.w tmp1, [result, #4]
    str.w tmp2, [result, #8]
    str.w tmp3, [result, #12]

    sub vec, #32                
    sub rmat, #32*64-16         
                                
    push.w {result}
    push.w {vec}

    mov.w accu0, #0
    mov.w accu1, #0
    mov.w accu2, #0
    mov.w accu3, #0
    
    vmov ctr1, #4
    1:
        pop.w {vec}
        ldr.w buf_vec, [vec], #4
        push.w {vec}
        .set kk, 0
        .rept 8                                 
            ldr.w tmp1, [rmat, #4]             
            ldr.w tmp2, [rmat, #8]
            ldr.w tmp3, [rmat, #12]
            ldr.w tmp0, [rmat], #32
            // bitslice on the fly
            bitslice mat0, mat1, mat2, mat3, tmp0, tmp1, tmp2, tmp3
            madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, buf_vec, tmp0, tmp1, tmp2, tmp3
            .if kk != 7
            lsr.w buf_vec, buf_vec, #4
            .endif
            .set kk, kk+1
        .endr
    vsub.f32 ctr1, ctr1, one
    vcmp.f32 ctr1, #0.0
    vmrs apsr_nzcv, FPSCR
    bhi.w 1b

    pop.w {vec}
    pop.w {result}
    // un-bitslice on the fly
    bitslice tmp0, tmp1, tmp2, tmp3, accu0, accu1, accu2, accu3
    str.w tmp0, [result, #16]
    str.w tmp1, [result, #20]
    str.w tmp2, [result, #24]
    str.w tmp3, [result, #28]

    pop.w {r4-r11, pc}