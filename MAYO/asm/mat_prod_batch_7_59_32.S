.syntax unified
.cpu cortex-m4
.thumb

#include "bitslice.inc"
#include "madd_bitsliced.inc"

//extern void mat_prod_batch_7_59_32(uint32_t *c, uint32_t *a, uint8_t *b);
.global mat_prod_batch_7_59_32
.type mat_prod_batch_7_59_32, %function
.align 2
mat_prod_batch_7_59_32:
    push {r4-r11, r14}
    c_ptr .req r0
    a_ptr .req r1
    b_ptr .req r2
    mat0    .req r0
    mat1    .req r2
    mat2    .req r10
    mat3    .req r12
    b_32     .req r5
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

    push.w {c_ptr}
    push.w {b_ptr}

    vmov one, #0.5
    vmov ctr1, #4                     
    vsub.f32 ctr1, ctr1, one        
    1:
        mov.w accu0, #0
        mov.w accu1, #0
        mov.w accu2, #0
        mov.w accu3, #0
        vmov ctr2, #4                     
        vsub.f32 ctr2, ctr2, one        
        2:
            pop.w {b_ptr}
            ldr.w b_32, [b_ptr], #4     
            push.w {b_ptr}
            .set kk, 0
            .rept 8
                ldr.w tmp0, [a_ptr]
                ldr.w tmp1, [a_ptr, #4]
                ldr.w tmp2, [a_ptr, #8]
                ldr.w tmp3, [a_ptr, #12]
                add a_ptr, #32*7
                // bitslice on the fly
                bitslice mat0, mat1, mat2, mat3, tmp0, tmp1, tmp2, tmp3

                madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
                .if kk != 7
                lsr.w b_32, b_32, #4
                .endif
                .set kk, kk+1
            .endr
        vsub.f32 ctr2, ctr2, one
        vcmp.f32 ctr2, #0.0
        vmrs apsr_nzcv, FPSCR
        bhi.w 2b

        pop.w {b_ptr}
        ldr.w b_32, [b_ptr]
        push.w {b_ptr}
        .set kk, 0
        .rept 3
            ldr.w tmp0, [a_ptr]
            ldr.w tmp1, [a_ptr, #4]
            ldr.w tmp2, [a_ptr, #8]
            ldr.w tmp3, [a_ptr, #12]

            add a_ptr, #32*7
            // bitslice on the fly
            bitslice mat0, mat1, mat2, mat3, tmp0, tmp1, tmp2, tmp3
            madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
            .if kk != 2
                lsr.w b_32, b_32, #4
            .endif
            .set kk, kk+1
            .endr

        pop.w {b_ptr}
        pop.w {c_ptr}
        // un-bitslice on the fly
        bitslice tmp0, tmp1, tmp2, tmp3, accu0, accu1, accu2, accu3

        str.w tmp1, [c_ptr, #4]
        str.w tmp2, [c_ptr, #8]
        str.w tmp3, [c_ptr, #12]
        str.w tmp0, [c_ptr], #16
        push.w {c_ptr}
        mov r3, #59*7*32-16
        sub a_ptr, a_ptr, r3
        sub b_ptr, b_ptr, #28
        push.w {b_ptr}

        mov.w accu0, #0
        mov.w accu1, #0
        mov.w accu2, #0
        mov.w accu3, #0
        vmov ctr2, #4                     
        vsub.f32 ctr2, ctr2, one        
        2:
            pop.w {b_ptr}
            ldr.w b_32, [b_ptr], #4     
            push.w {b_ptr}
            .set kk, 0
            .rept 8
                ldr.w tmp0, [a_ptr]
                ldr.w tmp1, [a_ptr, #4]
                ldr.w tmp2, [a_ptr, #8]
                ldr.w tmp3, [a_ptr, #12]
                add a_ptr, #32*7
                // bitslice on the fly
                bitslice mat0, mat1, mat2, mat3, tmp0, tmp1, tmp2, tmp3

                madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
                .if kk != 7
                lsr.w b_32, b_32, #4
                .endif
                .set kk, kk+1
            .endr
        vsub.f32 ctr2, ctr2, one
        vcmp.f32 ctr2, #0.0
        vmrs apsr_nzcv, FPSCR
        bhi.w 2b

        pop.w {b_ptr}
        ldr.w b_32, [b_ptr]
        push.w {b_ptr}
        .set kk, 0
        .rept 3
            ldr.w tmp0, [a_ptr]
            ldr.w tmp1, [a_ptr, #4]
            ldr.w tmp2, [a_ptr, #8]
            ldr.w tmp3, [a_ptr, #12]
            add a_ptr, #32*7
            // bitslice on the fly
            bitslice mat0, mat1, mat2, mat3, tmp0, tmp1, tmp2, tmp3
            madd_bitsliced accu0, accu1, accu2, accu3, mat0, mat1, mat2, mat3, b_32, tmp0, tmp1, tmp2, tmp3
            .if kk != 2
                lsr.w b_32, b_32, #4
            .endif
            .set kk, kk+1
            .endr

        pop.w {b_ptr}
        pop.w {c_ptr}
        // un-bitslice on the fly
        bitslice tmp0, tmp1, tmp2, tmp3, accu0, accu1, accu2, accu3

        str.w tmp1, [c_ptr, #4]
        str.w tmp2, [c_ptr, #8]
        str.w tmp3, [c_ptr, #12]
        str.w tmp0, [c_ptr], #16
        push.w {c_ptr}
        mov r3, #59*7*32-16 
        sub a_ptr, a_ptr, r3
        sub b_ptr, b_ptr, #28
        push.w {b_ptr}

        vsub.f32 ctr1, ctr1, one
        vcmp.f32 ctr1, #0.0
        vmrs apsr_nzcv, FPSCR
        bhi.w 1b

    add sp, #8
    pop.w {r4-r11, pc}
