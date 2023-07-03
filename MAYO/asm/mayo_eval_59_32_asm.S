.syntax unified
.cpu cortex-m4
.thumb

#include "bitslice.inc"
#include "madd_bitsliced.inc"

.global mayo_eval_59_32_asm
.type mayo_eval_59_32_asm, %function
.align 2

mayo_eval_59_32_asm:
    push.w {r4-r11, r14}

    one      .req s0
    ctr1     .req s1
    ctr2     .req s2

    sub.w sp, #92								
    # initialize y with 0
    mov.w r12, #0
    str.w r12, [sp, #60]
    str.w r12, [sp, #64]
    str.w r12, [sp, #68]
    str.w r12, [sp, #72]
	str.w r12, [sp, #76]
    str.w r12, [sp, #80]
    str.w r12, [sp, #84]
    str.w r12, [sp, #88]

    # re-organize x
    .set j, 0
	.rept 7										
		ldr.w r3, [r2, #4*j]
		and.w r4, r3, #0xF0F0F0F0
		lsr.w r4, r4, #4
		and.w r3, r3, #0x0F0F0F0F

		.set i, 0
		.rept 4
			strb.w r3, [sp, #2*i+8*j]
			strb.w r4, [sp, #2*i+8*j+1]			
			lsr.w r3, r3, #8
			lsr.w r4, r4, #8
			.set i, i+1
		.endr
		.set j, j+1
    .endr
	
	ldrh.w r3, [r2, #28]				
    and.w r4, r3, #0xF0F0F0F0
    lsr.w r4, r4, #4
    and.w r3, r3, #0x0F0F0F0F
    .set i, 0
    .rept 2
		strb.w r3, [sp, #28*2 + i*2]
		strb.w r4, [sp, #28*2 + i*2+1]			
		lsr.w r3, r3, #8
		lsr.w r4, r4, #8
		.set i, i+1
    .endr
	
    mov.w r2, sp

	vmov s6, r1						

    vmov.w one, #0.5
    vmov.w ctr1, #30							
	vsub.f32 ctr1, ctr1, one

    vmov.w s3, r0
    2:
		vmov.w s4, r2 							
		mov.w r4, #0
		mov.w r5, #0
		mov.w r6, #0
		mov.w r7, #0
		vmov.f32 ctr2, ctr1
		1:
			ldr.w r9,  [r1, #4]
			ldr.w r10, [r1, #8]
			ldr.w r11, [r1, #12]
			ldr.w r8,  [r1], #32
			ldrb.w r0, [r2], #1

			vmov s5, r2
			bitslice r12, r14, r3, r2, r8, r9, r10, r11
			madd_bitsliced r4, r5, r6, r7, r12, r14, r3, r2, r0, r8, r9, r10, r11
			vmov r2, s5

			vsub.f32 ctr2, ctr2, one
			vcmp.f32 ctr2, #0.0
			vmrs apsr_nzcv, FPSCR
			bhi.w 1b

		vmov r2, s4					
		ldr.w r8,  [sp, #60]		
		ldr.w r9,  [sp, #64]
		ldr.w r10, [sp, #68]
		ldr.w r11, [sp, #72]
		ldrb.w r3, [r2], #1
		vmov s5, r2
		madd_bitsliced r8, r9, r10, r11, r4, r5, r6, r7, r3, r0, r2, r12, r14
		vmov r2, s5

		str.w r8,  [sp, #60]
		str.w r9,  [sp, #64]
		str.w r10, [sp, #68]
		str.w r11, [sp, #72]
		vsub.f32 ctr1, ctr1, one
		vcmp.f32 ctr1, #0.0
		vmrs apsr_nzcv, FPSCR
		bhi.w 2b

    vmov r0, s3
    # un-bitslice
    bitslice r5, r2, r3, r4, r8, r9, r10, r11
    str.w r5, [r0]
    str.w r2, [r0, #4]
    str.w r3, [r0, #8]
    str.w r4, [r0, #12]
	
	vmov r1, s6						
	add.w r1, r1, #16
	
	vmov.w one, #0.5
    vmov.w ctr1, #30							
	vsub.f32 ctr1, ctr1, one

	mov.w r2, sp
	
    vmov.w s3, r0
    2:
		vmov.w s4, r2
		mov.w r4, #0
		mov.w r5, #0
		mov.w r6, #0
		mov.w r7, #0
		vmov.f32 ctr2, ctr1
		1:
			ldr.w r9,  [r1, #4]
			ldr.w r10, [r1, #8]
			ldr.w r11, [r1, #12]
			ldr.w r8,  [r1], #32
			ldrb.w r0, [r2], #1

			vmov s5, r2
			bitslice r12, r14, r3, r2, r8, r9, r10, r11
			madd_bitsliced r4, r5, r6, r7, r12, r14, r3, r2, r0, r8, r9, r10, r11
			vmov r2, s5

			vsub.f32 ctr2, ctr2, one
			vcmp.f32 ctr2, #0.0
			vmrs apsr_nzcv, FPSCR
			bhi.w 1b

		vmov r2, s4
		ldr.w r8,  [sp, #76]
		ldr.w r9,  [sp, #80]
		ldr.w r10, [sp, #84]
		ldr.w r11, [sp, #88]
		ldrb.w r3, [r2], #1
		vmov s5, r2
		madd_bitsliced r8, r9, r10, r11, r4, r5, r6, r7, r3, r0, r2, r12, r14
		vmov r2, s5

		str.w r8,  [sp, #76]
		str.w r9,  [sp, #80]
		str.w r10, [sp, #84]
		str.w r11, [sp, #88]
		vsub.f32 ctr1, ctr1, one
		vcmp.f32 ctr1, #0.0
		vmrs apsr_nzcv, FPSCR
		bhi.w 2b

    vmov r0, s3
    # un-bitslice
    bitslice r5, r2, r3, r4, r8, r9, r10, r11
    str.w r5, [r0,#16]
    str.w r2, [r0,#20]
    str.w r3, [r0,#24]
    str.w r4, [r0,#28]

    add.w sp, #92
    pop.w {r4-r11, pc}
