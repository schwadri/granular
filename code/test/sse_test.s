	.section	__TEXT,__text,regular,pure_instructions
	.globl	_main
	.align	4, 0x90
_main:
Leh_func_begin1:
	pushq	%rbp
Ltmp0:
	movq	%rsp, %rbp
Ltmp1:
	subq	$128, %rsp
Ltmp2:
	movl	$1000000, -44(%rbp)
	movl	-44(%rbp), %eax
	imull	$4, %eax, %eax
	movslq	%eax, %rax
	movq	%rax, -40(%rbp)
	movq	-40(%rbp), %rax
	movabsq	$4, %rcx
	imulq	%rcx, %rax
	movq	%rax, %rdi
	movq	%rcx, -112(%rbp)
	callq	__Znam
	movq	%rax, -56(%rbp)
	movl	-44(%rbp), %eax
	imull	$4, %eax, %eax
	movslq	%eax, %rax
	movq	%rax, -32(%rbp)
	movq	-32(%rbp), %rax
	movq	-112(%rbp), %rcx
	imulq	%rcx, %rax
	movq	%rax, %rdi
	callq	__Znam
	movq	%rax, -64(%rbp)
	leaq	-80(%rbp), %rax
	movq	%rax, %rdi
	callq	__ZN4util8time_posC1Ev
	movq	-72(%rbp), %rax
	movq	-80(%rbp), %rcx
	movq	__ZSt4cout@GOTPCREL(%rip), %rdx
	leaq	(%rdx), %rdx
	movq	%rdx, %rdi
	movq	%rcx, %rsi
	movq	%rax, -120(%rbp)
	callq	__ZNSolsEx
	leaq	L_.str(%rip), %rcx
	movq	%rax, %rdi
	movq	%rcx, %rsi
	callq	__ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc
	movq	%rax, %rdi
	movq	-120(%rbp), %rsi
	callq	__ZNSolsEx
	leaq	L_.str1(%rip), %rcx
	movq	%rax, %rdi
	movq	%rcx, %rsi
	callq	__ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc
	movq	__ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_@GOTPCREL(%rip), %rcx
	leaq	(%rcx), %rcx
	movq	%rax, %rdi
	movq	%rcx, %rsi
	callq	__ZNSolsEPFRSoS_E
	movl	$0, -100(%rbp)
	jmp	LBB1_2
LBB1_1:
	movq	-56(%rbp), %rax
	movl	-100(%rbp), %ecx
	movslq	%ecx, %rcx
	movabsq	$4, %rdx
	imulq	%rdx, %rcx
	addq	%rcx, %rax
	movq	-56(%rbp), %rcx
	movl	-100(%rbp), %edx
	movslq	%edx, %rdx
	movabsq	$4, %rsi
	imulq	%rsi, %rdx
	addq	%rdx, %rcx
	movaps	(%rcx), %xmm0
	movq	-64(%rbp), %rcx
	movslq	-100(%rbp), %rdx
	movaps	(%rcx,%rdx,4), %xmm1
	addps	%xmm1, %xmm0
	movaps	%xmm0, (%rax)
	movl	-100(%rbp), %eax
	leal	4(%rax), %eax
	movl	%eax, -100(%rbp)
LBB1_2:
	movl	-100(%rbp), %eax
	movl	-44(%rbp), %ecx
	cmpl	%ecx, %eax
	jl	LBB1_1
	leaq	-24(%rbp), %rax
	movq	%rax, %rdi
	callq	__ZN4util8time_posC1Ev
	leaq	-24(%rbp), %rax
	leaq	-80(%rbp), %rcx
	movq	%rax, %rdi
	movq	%rcx, %rsi
	callq	__ZN4utilmiERKNS_8time_posES2_
	leaq	-96(%rbp), %rcx
	movq	%rax, (%rcx)
	movq	%rdx, -88(%rbp)
	movq	-96(%rbp), %rax
	movq	__ZSt4cout@GOTPCREL(%rip), %rcx
	movq	%rcx, %rdi
	movq	%rax, %rsi
	movq	%rdx, -128(%rbp)
	callq	__ZNSolsEx
	leaq	L_.str(%rip), %rcx
	movq	%rax, %rdi
	movq	%rcx, %rsi
	callq	__ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc
	movq	%rax, %rdi
	movq	-128(%rbp), %rsi
	callq	__ZNSolsEx
	leaq	L_.str1(%rip), %rcx
	movq	%rax, %rdi
	movq	%rcx, %rsi
	callq	__ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc
	movq	__ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_@GOTPCREL(%rip), %rcx
	movq	%rax, %rdi
	movq	%rcx, %rsi
	callq	__ZNSolsEPFRSoS_E
	movl	$0, -104(%rbp)
	jmp	LBB1_5
LBB1_4:
	movl	-104(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -104(%rbp)
LBB1_5:
	movl	-104(%rbp), %eax
	movl	-44(%rbp), %ecx
	cmpl	%ecx, %eax
	jl	LBB1_4
	movq	-56(%rbp), %rax
	cmpq	$0, %rax
	je	LBB1_8
	movq	-56(%rbp), %rax
	movq	%rax, %rdi
	callq	__ZdaPv
LBB1_8:
	movq	-64(%rbp), %rax
	cmpq	$0, %rax
	je	LBB1_10
	movq	-64(%rbp), %rax
	movq	%rax, %rdi
	callq	__ZdaPv
LBB1_10:
	movl	$0, -8(%rbp)
	movl	-8(%rbp), %eax
	movl	%eax, -4(%rbp)
	movl	-4(%rbp), %eax
	addq	$128, %rsp
	popq	%rbp
	ret
Leh_func_end1:

	.section	__TEXT,__StaticInit,regular,pure_instructions
	.align	4, 0x90
__GLOBAL__I_main:
Leh_func_begin2:
	pushq	%rbp
Ltmp3:
	movq	%rsp, %rbp
Ltmp4:
	movl	$1, %eax
	movl	$65535, %ecx
	movl	%eax, %edi
	movl	%ecx, %esi
	callq	__Z41__static_initialization_and_destruction_0ii
	popq	%rbp
	ret
Leh_func_end2:

	.section	__TEXT,__textcoal_nt,coalesced,pure_instructions
	.globl	__ZN4utilmiERKNS_8time_posES2_
.weak_definition __ZN4utilmiERKNS_8time_posES2_
	.align	4, 0x90
__ZN4utilmiERKNS_8time_posES2_:
Leh_func_begin3:
	pushq	%rbp
Ltmp5:
	movq	%rsp, %rbp
Ltmp6:
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	-16(%rbp), %rcx
	movq	(%rcx), %rcx
	subq	%rcx, %rax
	movq	%rax, -64(%rbp)
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	-16(%rbp), %rcx
	movq	8(%rcx), %rcx
	subq	%rcx, %rax
	movq	%rax, -56(%rbp)
	movq	-56(%rbp), %rax
	cmpq	$0, %rax
	jge	LBB3_2
	movq	-64(%rbp), %rax
	movabsq	$1, %rcx
	subq	%rcx, %rax
	movq	%rax, -64(%rbp)
	movq	-56(%rbp), %rax
	movabsq	$1000000, %rcx
	addq	%rcx, %rax
	movq	%rax, -56(%rbp)
LBB3_2:
	movq	-64(%rbp), %rax
	movq	%rax, -48(%rbp)
	movq	-56(%rbp), %rax
	movq	%rax, -40(%rbp)
	movq	-48(%rbp), %rax
	movq	%rax, -32(%rbp)
	movq	-40(%rbp), %rax
	movq	%rax, -24(%rbp)
	movq	-32(%rbp), %rax
	movq	-24(%rbp), %rcx
	movq	%rcx, %rdx
	popq	%rbp
	ret
Leh_func_end3:

	.section	__TEXT,__StaticInit,regular,pure_instructions
	.align	4, 0x90
__Z41__static_initialization_and_destruction_0ii:
Leh_func_begin4:
	pushq	%rbp
Ltmp7:
	movq	%rsp, %rbp
Ltmp8:
	subq	$16, %rsp
Ltmp9:
	movl	%edi, -4(%rbp)
	movl	%esi, -8(%rbp)
	movl	-4(%rbp), %eax
	cmpl	$1, %eax
	jne	LBB4_3
	movl	-8(%rbp), %eax
	cmpl	$65535, %eax
	jne	LBB4_3
	leaq	__ZStL8__ioinit(%rip), %rax
	movq	%rax, %rdi
	callq	__ZNSt8ios_base4InitC1Ev
	leaq	___tcf_0(%rip), %rax
	movabsq	$0, %rcx
	movq	___dso_handle@GOTPCREL(%rip), %rdx
	leaq	(%rdx), %rdx
	movq	%rax, %rdi
	movq	%rcx, %rsi
	callq	___cxa_atexit
LBB4_3:
	addq	$16, %rsp
	popq	%rbp
	ret
Leh_func_end4:

	.section	__TEXT,__text,regular,pure_instructions
	.align	4, 0x90
___tcf_0:
Leh_func_begin5:
	pushq	%rbp
Ltmp10:
	movq	%rsp, %rbp
Ltmp11:
	subq	$16, %rsp
Ltmp12:
	movq	%rdi, -8(%rbp)
	leaq	__ZStL8__ioinit(%rip), %rax
	movq	%rax, %rdi
	callq	__ZNSt8ios_base4InitD1Ev
	addq	$16, %rsp
	popq	%rbp
	ret
Leh_func_end5:

	.section	__TEXT,__textcoal_nt,coalesced,pure_instructions
	.globl	__ZN4util8time_posC1Ev
.weak_definition __ZN4util8time_posC1Ev
	.align	1, 0x90
__ZN4util8time_posC1Ev:
Leh_func_begin6:
	pushq	%rbp
Ltmp13:
	movq	%rsp, %rbp
Ltmp14:
	subq	$32, %rsp
Ltmp15:
	movq	%rdi, -8(%rbp)
	leaq	-24(%rbp), %rax
	movabsq	$0, %rcx
	movq	%rax, %rdi
	movq	%rcx, %rsi
	callq	_gettimeofday
	movq	-24(%rbp), %rax
	movq	-8(%rbp), %rcx
	movq	%rax, (%rcx)
	movl	-16(%rbp), %eax
	movslq	%eax, %rax
	movq	-8(%rbp), %rcx
	movq	%rax, 8(%rcx)
	addq	$32, %rsp
	popq	%rbp
	ret
Leh_func_end6:

.zerofill __DATA,__bss,__ZStL8__ioinit,1,3
	.section	__TEXT,__cstring,cstring_literals
L_.str:
	.asciz	 " s "

L_.str1:
	.asciz	 " us "

	.section	__DATA,__mod_init_func,mod_init_funcs
	.align	3
	.quad	__GLOBAL__I_main
	.section	__TEXT,__eh_frame,coalesced,no_toc+strip_static_syms+live_support
EH_frame0:
Lsection_eh_frame:
Leh_frame_common:
Lset0 = Leh_frame_common_end-Leh_frame_common_begin
	.long	Lset0
Leh_frame_common_begin:
	.long	0
	.byte	1
	.asciz	 "zR"
	.byte	1
	.byte	120
	.byte	16
	.byte	1
	.byte	16
	.byte	12
	.byte	7
	.byte	8
	.byte	144
	.byte	1
	.align	3
Leh_frame_common_end:
	.globl	_main.eh
_main.eh:
Lset1 = Leh_frame_end1-Leh_frame_begin1
	.long	Lset1
Leh_frame_begin1:
Lset2 = Leh_frame_begin1-Leh_frame_common
	.long	Lset2
Ltmp16:
	.quad	Leh_func_begin1-Ltmp16
Lset3 = Leh_func_end1-Leh_func_begin1
	.quad	Lset3
	.byte	0
	.byte	4
Lset4 = Ltmp0-Leh_func_begin1
	.long	Lset4
	.byte	14
	.byte	16
	.byte	134
	.byte	2
	.byte	4
Lset5 = Ltmp1-Ltmp0
	.long	Lset5
	.byte	13
	.byte	6
	.align	3
Leh_frame_end1:

__GLOBAL__I_main.eh:
Lset6 = Leh_frame_end2-Leh_frame_begin2
	.long	Lset6
Leh_frame_begin2:
Lset7 = Leh_frame_begin2-Leh_frame_common
	.long	Lset7
Ltmp17:
	.quad	Leh_func_begin2-Ltmp17
Lset8 = Leh_func_end2-Leh_func_begin2
	.quad	Lset8
	.byte	0
	.byte	4
Lset9 = Ltmp3-Leh_func_begin2
	.long	Lset9
	.byte	14
	.byte	16
	.byte	134
	.byte	2
	.byte	4
Lset10 = Ltmp4-Ltmp3
	.long	Lset10
	.byte	13
	.byte	6
	.align	3
Leh_frame_end2:

	.globl	__ZN4utilmiERKNS_8time_posES2_.eh
.weak_definition __ZN4utilmiERKNS_8time_posES2_.eh
__ZN4utilmiERKNS_8time_posES2_.eh:
Lset11 = Leh_frame_end3-Leh_frame_begin3
	.long	Lset11
Leh_frame_begin3:
Lset12 = Leh_frame_begin3-Leh_frame_common
	.long	Lset12
Ltmp18:
	.quad	Leh_func_begin3-Ltmp18
Lset13 = Leh_func_end3-Leh_func_begin3
	.quad	Lset13
	.byte	0
	.byte	4
Lset14 = Ltmp5-Leh_func_begin3
	.long	Lset14
	.byte	14
	.byte	16
	.byte	134
	.byte	2
	.byte	4
Lset15 = Ltmp6-Ltmp5
	.long	Lset15
	.byte	13
	.byte	6
	.align	3
Leh_frame_end3:

__Z41__static_initialization_and_destruction_0ii.eh:
Lset16 = Leh_frame_end4-Leh_frame_begin4
	.long	Lset16
Leh_frame_begin4:
Lset17 = Leh_frame_begin4-Leh_frame_common
	.long	Lset17
Ltmp19:
	.quad	Leh_func_begin4-Ltmp19
Lset18 = Leh_func_end4-Leh_func_begin4
	.quad	Lset18
	.byte	0
	.byte	4
Lset19 = Ltmp7-Leh_func_begin4
	.long	Lset19
	.byte	14
	.byte	16
	.byte	134
	.byte	2
	.byte	4
Lset20 = Ltmp8-Ltmp7
	.long	Lset20
	.byte	13
	.byte	6
	.align	3
Leh_frame_end4:

___tcf_0.eh:
Lset21 = Leh_frame_end5-Leh_frame_begin5
	.long	Lset21
Leh_frame_begin5:
Lset22 = Leh_frame_begin5-Leh_frame_common
	.long	Lset22
Ltmp20:
	.quad	Leh_func_begin5-Ltmp20
Lset23 = Leh_func_end5-Leh_func_begin5
	.quad	Lset23
	.byte	0
	.byte	4
Lset24 = Ltmp10-Leh_func_begin5
	.long	Lset24
	.byte	14
	.byte	16
	.byte	134
	.byte	2
	.byte	4
Lset25 = Ltmp11-Ltmp10
	.long	Lset25
	.byte	13
	.byte	6
	.align	3
Leh_frame_end5:

	.globl	__ZN4util8time_posC1Ev.eh
.weak_definition __ZN4util8time_posC1Ev.eh
__ZN4util8time_posC1Ev.eh:
Lset26 = Leh_frame_end6-Leh_frame_begin6
	.long	Lset26
Leh_frame_begin6:
Lset27 = Leh_frame_begin6-Leh_frame_common
	.long	Lset27
Ltmp21:
	.quad	Leh_func_begin6-Ltmp21
Lset28 = Leh_func_end6-Leh_func_begin6
	.quad	Lset28
	.byte	0
	.byte	4
Lset29 = Ltmp13-Leh_func_begin6
	.long	Lset29
	.byte	14
	.byte	16
	.byte	134
	.byte	2
	.byte	4
Lset30 = Ltmp14-Ltmp13
	.long	Lset30
	.byte	13
	.byte	6
	.align	3
Leh_frame_end6:


.subsections_via_symbols
