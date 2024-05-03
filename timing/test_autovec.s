	.arch armv8-a
	.file	"test_autovec.c"
// GNU C17 (Debian 12.2.0-14) version 12.2.0 (aarch64-linux-gnu)
//	compiled by GNU C version 12.2.0, GMP version 6.2.1, MPFR version 4.1.1-p1, MPC version 1.3.1, isl version isl-0.25-GMP

// warning: MPFR header version 4.1.1-p1 differs from library version 4.2.0.
// GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
// options passed: -mlittle-endian -mabi=lp64 -O2 -fasynchronous-unwind-tables
	.text
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align	3
.LC0:
	.string	"%f "
	.section	.text.startup,"ax",@progbits
	.align	2
	.p2align 4,,11
	.global	main
	.type	main, %function
main:
.LFB11:
	.cfi_startproc
	mov	x12, 8240	//,
	sub	sp, sp, x12	//,,
	.cfi_def_cfa_offset 8240
// test_autovec.c:7:     double a[1024] = {0};
	add	x3, sp, 48	// tmp109,,
	mov	x2, 8192	//,
	mov	x0, x3	//, tmp109
	mov	w1, 0	//,
// test_autovec.c:4: int main() {
	stp	x29, x30, [sp]	//,,
	.cfi_offset 29, -8240
	.cfi_offset 30, -8232
	mov	x29, sp	//,
	stp	x19, x20, [sp, 16]	//,,
	str	x21, [sp, 32]	//,
	.cfi_offset 19, -8224
	.cfi_offset 20, -8216
	.cfi_offset 21, -8208
// test_autovec.c:7:     double a[1024] = {0};
	bl	memset		//
	mov	x3, x0	// tmp109,
	mov	x19, x0	// ivtmp.25, tmp109
	adrp	x0, .LC1	// tmp125,
	add	x20, x3, 1024	// _54, tmp109,
	movi	v3.4s, 0x4	// tmp112
	ldr	q2, [x0, #:lo12:.LC1]	// vect_vec_iv_.12,
	mov	x0, x3	// ivtmp.35, tmp109
	.p2align 3,,7
.L2:
	mov	v0.16b, v2.16b	// vect_vec_iv_.12, vect_vec_iv_.12
	add	v2.4s, v2.4s, v3.4s	// vect_vec_iv_.12, vect_vec_iv_.12, tmp112
// test_autovec.c:11: 	a[i] = i;
	sxtl	v1.2d, v0.2s	// vect__1.13, vect_vec_iv_.12
	sxtl2	v0.2d, v0.4s	// vect__1.13, vect_vec_iv_.12
	scvtf	v1.2d, v1.2d	// tmp114, vect__1.13
	scvtf	v0.2d, v0.2d	// tmp116, vect__1.13
	stp	q1, q0, [x0]	// tmp114, tmp116, MEM <vector(2) double> [(double *)_52]
	add	x0, x0, 32	// ivtmp.35, ivtmp.35,
	cmp	x0, x20	// ivtmp.35, _54
	bne	.L2		//,
	mov	x0, x3	// ivtmp.29, tmp109
	.p2align 3,,7
.L3:
// test_autovec.c:15: 	a[i] *= 2;
	ldr	q0, [x0]	// MEM <vector(2) double> [(double *)_18], MEM <vector(2) double> [(double *)_18]
	fadd	v0.2d, v0.2d, v0.2d	// vect__3.9, MEM <vector(2) double> [(double *)_18], MEM <vector(2) double> [(double *)_18]
	str	q0, [x0], 16	// vect__3.9, MEM <vector(2) double> [(double *)_18]
	cmp	x0, x20	// ivtmp.29, _54
	bne	.L3		//,
	adrp	x21, .LC0	// tmp123,
// test_autovec.c:19: 	printf("%f ", a[i]);
	add	x21, x21, :lo12:.LC0	// tmp120, tmp123,
	.p2align 3,,7
.L4:
// test_autovec.c:19: 	printf("%f ", a[i]);
	ldr	d0, [x19], 8	//, MEM[(double *)_20]
	mov	x0, x21	//, tmp120
	bl	printf		//
// test_autovec.c:18:     for(int i = 0; i < LEN; i++) {
	cmp	x19, x20	// ivtmp.25, _54
	bne	.L4		//,
// test_autovec.c:21:     printf("\n");
	mov	w0, 10	//,
	bl	putchar		//
// test_autovec.c:22: }
	ldp	x29, x30, [sp]	//,,
	mov	w0, 0	//,
	ldp	x19, x20, [sp, 16]	//,,
	mov	x12, 8240	//,
	ldr	x21, [sp, 32]	//,
	add	sp, sp, x12	//,,
	.cfi_restore 29
	.cfi_restore 30
	.cfi_restore 21
	.cfi_restore 19
	.cfi_restore 20
	.cfi_def_cfa_offset 0
	ret	
	.cfi_endproc
.LFE11:
	.size	main, .-main
	.section	.rodata.cst16,"aM",@progbits,16
	.align	4
.LC1:
	.word	0
	.word	1
	.word	2
	.word	3
	.ident	"GCC: (Debian 12.2.0-14) 12.2.0"
	.section	.note.GNU-stack,"",@progbits
