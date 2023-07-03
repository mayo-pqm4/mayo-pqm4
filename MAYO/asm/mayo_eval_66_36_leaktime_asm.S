.syntax unified
.cpu cortex-m4
.thumb

#include "bitslice.inc"
#include "madd_bitsliced.inc"

.global mayo_eval_66_36_leaktime_asm
.type mayo_eval_66_36_leaktime_asm, %function
.align 2
mayo_eval_66_36_leaktime_asm:
	push.w {r4-r11, r14}				
	ctr1   .req r0						
	trimat .req r1						
	ctr2   .req r2						
	xi	   .req r3						

	xj     .req r5						
	lutgf    .req r0
	mat0     .req r6
	mat1     .req r7
	mat2     .req r8
	mat3     .req r9
	buf0     .req r10
	buf1     .req r11
	buf2     .req r12
	buf3     .req r14
	tmp    	 .req sp					
	xptr     .req r4					
	
	sub.w sp, sp, #540+36+68			
										
	add.w xptr, sp, #540+36				
	vmov.w s0, r0						
	vmov.w s2, r3						

	// init tmp to zero
	mov.w r0, #0						
	.set i, 0							
	.rept 67							
		strd.w r0, r0, [sp, #36+i*8]	
		.set i, i+1						
	.endr
	str.w r0, [sp, #36+536]

	.set i, 0							
	.rept 8
										
		ldr.w r0, [r2, #i*4]				
		and.w r3, r0, #0xF					
		and.w r5, r0, #0xF0					
		add.w r3, r3, r5, lsl#4				
		and.w r5, r0, #0xF00				
		add.w r3, r3, r5, lsl#8			
		and.w r5, r0, #0xF000				
		add.w r3, r3, r5, lsl#12			
		str.w r3, [xptr, #i*8]				
		ubfx.w r3, r0, #16, #4				
		and.w r5, r0, #0xF00000				
		add.w r3, r3, r5, lsr#12			
		and.w r5, r0, #0xF000000
		add.w r3, r3, r5, lsr#8
		and.w r5, r0, #0xF0000000
		add.w r3, r3, r5, lsr#4
		str.w r3, [xptr, #i*8+4]			
		.set i, i+1
	.endr								
								
	ldrb.w r0, [r2, #32]				
	and.w r3, r0, #0xF
	and.w r5, r0, #0xF0
	add.w r3, r3, r5, lsl#4				
	str.w r3, [xptr, #64]				

	mov.w ctr1, #66	
	vmov.w s1, xptr						
	outer: 	
		vmov.w xptr, s1					
		ldrb.w xi, [xptr], #1			
		vmov.w s1, xptr					

		// if 0, do nothing				
		cmp.w xi, #0		
		beq.w skip_outer

		mov.w ctr2, ctr1				
		vmov.w s3, ctr1					
		vmov.w lutgf, s2				
		sub.w xptr, xptr, #1			
		inner: 
			ldrb.w xj, [xptr], #1		
			// if 0, do nothing
			cmp.w xj, #0
			beq.w skip_inner			
			orr.w xj, xi, xj, lsl#4		
			ldrb.w xj, [lutgf, xj]	

			// compute address of buffer
			lsl.w xj, xj, #2        	
            add.w xj, xj, xj, lsl#3     
            add.w xj, sp, xj           
			
			ldr.w buf0, [xj, #0]        
			ldr.w buf1, [xj, #4]
			ldr.w buf2, [xj, #8]
			ldr.w buf3, [xj, #12]
			ldr.w mat1, [trimat, #4]
			ldr.w mat2, [trimat, #8]
			ldr.w mat3, [trimat, #12]
			ldr.w mat0, [trimat], #16
			eor.w buf0, buf0, mat0		
			str.w buf0, [xj, #0]		
			eor.w buf1, buf1, mat1
			str.w buf1, [xj, #4]
			eor.w buf2, buf2, mat2
			str.w buf2, [xj, #8]
			eor.w buf3, buf3, mat3
			str.w buf3, [xj, #12]

			ldr.w buf0, [xj, #16]		
			ldr.w buf1, [xj, #20]
			ldr.w buf2, [xj, #24]
			ldr.w buf3, [xj, #28]
			ldr.w mat1, [trimat, #4]
			ldr.w mat2, [trimat, #8]
			ldr.w mat3, [trimat, #12]
			ldr.w mat0, [trimat], #16
			eor.w buf0, buf0, mat0
			str.w buf0, [xj, #16]
			eor.w buf1, buf1, mat1
			str.w buf1, [xj, #20]
			eor.w buf2, buf2, mat2
			str.w buf2, [xj, #24]
			eor.w buf3, buf3, mat3
			str.w buf3, [xj, #28]			

            ldrh.w buf0, [xj, #32]          
            ldrh.w mat0, [trimat], #2       
            eor.w buf0, buf0, mat0 
            strh.w buf0, [xj, #32]          

			subs.w ctr2, ctr2, #1
			bne.w inner
			cont_inner:

		vmov.w ctr1, s3
		subs.w ctr1, ctr1, #1			
		bne.w outer
	cont_outer:
	// do the actual multiplication
	add.w r1, sp, #36					
    ldr.w r7, [r1, #4]					
    ldr.w r8, [r1, #8]
    ldr.w r9, [r1, #12]
    ldr.w r6, [r1], #36                 
    bitslice r2, r3, r4, r5, r6, r7, r8, r9     
    mov.w r0, #2								
    1:                                      
        ldr.w r7, [r1, #4]
        ldr.w r8, [r1, #8]
        ldr.w r9, [r1, #12]
        ldr.w r6, [r1], #36

        bitslice r10, r11, r12, r14, r6, r7, r8, r9
        madd_bitsliced r2, r3, r4, r5, r10, r11, r12, r14, r0, r6, r7, r8, r9

        add.w r0, r0, #1				
        cmp.w r0, #16					
        bne.w 1b

    bitslice r6, r7, r8, r9, r2, r3, r4, r5		
	vmov.w r0, s0						
    str.w r6, [r0]                      
    str.w r7, [r0, #4]
    str.w r8, [r0, #8]
    str.w r9, [r0, #12]
    
    sub.w r1, r1, #36*15-16				
    ldr.w r7, [r1, #4]
    ldr.w r8, [r1, #8]
    ldr.w r9, [r1, #12]
    ldr.w r6, [r1], #36
    bitslice r2, r3, r4, r5, r6, r7, r8, r9
    mov.w r0, #2
    1:
        ldr.w r7, [r1, #4]
        ldr.w r8, [r1, #8]
        ldr.w r9, [r1, #12]
        ldr.w r6, [r1], #36
        bitslice r10, r11, r12, r14, r6, r7, r8, r9
        madd_bitsliced r2, r3, r4, r5, r10, r11, r12, r14, r0, r6, r7, r8, r9

        add.w r0, r0, #1
        cmp.w r0, #16
        bne.w 1b

    bitslice r6, r7, r8, r9, r2, r3, r4, r5
    vmov.w r0, s0
    str.w r6, [r0, #16]                 
    str.w r7, [r0, #20]
    str.w r8, [r0, #24]
    str.w r9, [r0, #28]

    sub.w r1, r1, #36*15-16             
    ldrh.w r6, [r1], #36               	
    bitslice_single r2, r3, r4, r5, r6    
    mov.w r0, #2								
    1:                                      
        ldrh.w r6, [r1], #36                

        bitslice_single r10, r11, r12, r14, r6
        madd_bitsliced r2, r3, r4, r5, r10, r11, r12, r14, r0, r6, r7, r8, r9

        add.w r0, r0, #1				
        cmp.w r0, #16					
        bne.w 1b

    unbitslice_single r6, r2, r3, r4, r5		
	vmov.w r0, s0						        
    strh.w r6, [r0, #32]                        

	add.w sp, sp, #540+36+68
	pop.w {r4-r11, pc}

	skip_inner:
		add.w trimat, trimat, #34               
		subs.w ctr2, ctr2, #1
		bne.w inner
		b cont_inner
	skip_outer:                                 
        add.w trimat, trimat,  ctr1, lsl#1      
        add.w trimat, trimat,  ctr1, lsl#5
		subs.w ctr1, ctr1, #1
		bne.w outer
		b cont_outer
		