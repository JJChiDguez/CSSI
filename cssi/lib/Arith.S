//##########################################
//# Arith.S: Assembly code for integer and
//#          field arithmetic
//# Created by: José Eduardo Ochoa Jiménez
//##########################################

// Registers that are used for parameter passing:
#define reg_p1  %rdi
#define reg_p2  %rsi
#define reg_p3  %rdx

prime:
         .quad  0xAC0E7A06FFFFFFFF
         .quad  0x12
mu:
         .quad  0xAC0E7A0700000001
         .quad  0x96F0AD1DFAEEAC43

//**************************************************************
//  Shift-right of 2-word operands
//  Operation: c [reg_p1] = a [reg_p2] div 2
//**************************************************************
.global SHR
.p2align 4
SHR:

 movq  (reg_p2), %rax
 movq 8(reg_p2), %rcx

 shrq  $1, %rcx
 rcrq  $1, %rax

 movq      %rax,   (reg_p1)
 movq      %rcx,  8(reg_p1)

 ret

//**************************************************************
//  Shift-left of 2-word operands
//  Operation: c [reg_p1] = a [reg_p2] * 2
//**************************************************************
.global SHL
.p2align 4
SHL:

 movq  (reg_p2), %rax
 movq 8(reg_p2), %rcx
 movq        $0, %r8

 shlq  $1, %rax
 rclq  $1, %rcx
 rclq  $1, %r8

 movq      %rax,   (reg_p1)
 movq      %rcx,  8(reg_p1)
 movq       %r8, 16(reg_p1)

 ret

//**************************************************************
//  Shift-left of 4-word operands
//  Operation: c [reg_p1] = a [reg_p2] * 2
//**************************************************************
.global SHL_D
.p2align 4
SHL_D:

 movq   (reg_p2), %rax
 movq  8(reg_p2), %rcx
 movq 16(reg_p2), %r8
 movq 24(reg_p2), %r9
 movq         $0, %r10

 shlq  $1, %rax
 rclq  $1, %rcx
 rclq  $1, %r8
 rclq  $1, %r9
 rclq  $1, %r10

 movq      %rax,   (reg_p1)
 movq      %rcx,  8(reg_p1)
 movq       %r8, 16(reg_p1)
 movq       %r9, 24(reg_p1)
 movq      %r10, 32(reg_p1)

 ret

//**************************************************************
//  Integer Addition of 2-word operands
//  Operation: c [reg_p1] = a [reg_p2] + b [reg_p3]
//**************************************************************
.global ADD
.p2align 4
ADD:

 movq  (reg_p2), %rax
 movq 8(reg_p2), %rcx
 movq        $0, %r8

 addq  (reg_p3), %rax
 adcq 8(reg_p3), %rcx
 adcq        $0, %r8

 movq      %rax,   (reg_p1)
 movq      %rcx,  8(reg_p1)
 movq      %r8 , 16(reg_p1)

 ret

//**************************************************************
//  Integer Addition of 4-word operands
//  Operation: c [reg_p1] = a [reg_p2] + b [reg_p3]
//**************************************************************
.global ADD_D
.p2align 4
ADD_D:

 movq   (reg_p2), %rax
 movq  8(reg_p2), %rcx
 movq 16(reg_p2), %r8
 movq 24(reg_p2), %r9
 movq        $0, %r10

 addq   (reg_p3), %rax
 adcq  8(reg_p3), %rcx
 adcq 16(reg_p3), %r8
 adcq 24(reg_p3), %r9
 adcq        $0, %r10

 movq      %rax,   (reg_p1)
 movq      %rcx,  8(reg_p1)
 movq      %r8 , 16(reg_p1)
 movq      %r9 , 24(reg_p1)
 movq      %r10, 32(reg_p1)

 ret


//**************************************************************
//  Integer Subtraction of 2-word operands
//  Operation: c [reg_p1] = a [reg_p2] - b [reg_p3]
//**************************************************************
.global SUB
.p2align 4
SUB:

 movq  (reg_p2), %rax
 movq 8(reg_p2), %rcx
 movq        $0, %r8

 subq  (reg_p3), %rax
 sbbq 8(reg_p3), %rcx
 sbbq        $0, %r8

 movq      %rax,   (reg_p1)
 movq      %rcx,  8(reg_p1)
 movq      %r8 , 16(reg_p1)

 ret

//**************************************************************
//  Integer Multiplication of 2-word operands
//  Operation: c [reg_p1] = a [reg_p2] * b [reg_p3]
//**************************************************************
.global MUL
.p2align 4
MUL:

 push %r12
 push %r13
 push %r14
 
 movq    reg_p3, %rcx

 // Integer multiplication
 // c = reg_p2*reg_p3
 movq  (reg_p2), %r8
 movq    (%rcx), %rax
 mulq  %r8
 movq      %rax, %r9
 movq      %rdx, %r10
 
 movq   8(%rcx), %rax
 mulq  %r8
 movq      %rax, %r11
 movq      %rdx, %r12

 movq 8(reg_p2), %r8
 movq    (%rcx), %rax
 mulq  %r8
 movq      %rax, %r13
 movq      %rdx, %r14
 
 movq   8(%rcx), %rax
 mulq  %r8

 addq %r11, %r10
 adcq %r12, %rax
 adcq   $0, %rdx
 addq %r13, %r10
 adcq %r14, %rax
 adcq   $0, %rdx

 movq  %r9,  (reg_p1)
 movq %r10, 8(reg_p1)
 movq %rax, 16(reg_p1)
 movq %rdx, 24(reg_p1)

 pop %r14
 pop %r13
 pop %r12
 
 ret

//**************************************************************
//  Integer Squaring of 2-word operands
//  Operation: c [reg_p1] = a [reg_p2]^2
//**************************************************************
.global SQR
.p2align 4
SQR:

 push %r13
 push %r14

 // Integer squaring
 // c = reg_p2*reg_p3
 movq  (reg_p2), %rax
 mulq  %rax
 movq      %rax, %r9
 movq      %rdx, %r10
 
 movq 8(reg_p2), %r8
 movq  (reg_p2), %rax
 mulq  %r8
 movq      %rax, %r13
 movq      %rdx, %r14
 
 movq  8(reg_p2), %rax
 mulq  %rax

 addq %r13, %r10
 adcq %r14, %rax
 adcq   $0, %rdx
 addq %r13, %r10
 adcq %r14, %rax
 adcq   $0, %rdx

 movq  %r9,  (reg_p1)
 movq %r10, 8(reg_p1)
 movq %rax, 16(reg_p1)
 movq %rdx, 24(reg_p1)

 pop %r14
 pop %r13
 
 ret

//**************************************************************
//  REDC: Montgomery modular reduction
//  Operation: c [reg_p1] = a [reg_p2] mod p
//**************************************************************
.global REDC
.p2align 4
REDC:

 push %r12
 push %r13
 push %r14
 push %r15
 
 // Montgomery reduction
 // m = mu*c mod 2^(64*2)
 movq  (reg_p2), %r9
 movq  (mu), %rax
 mulq  %r9
 movq      %rax, %r11
 movq      %rdx, %r12
 
 movq 8(reg_p2), %r10
 movq (mu), %rax
 mulq  %r10
 movq      %rax, %r13
 
 movq (mu+8), %rax
 mulq  %r9

 addq %r13, %r12
 addq %rax, %r12

 // u = m*p
 movq (prime), %rax
 mulq  %r11
 movq      %rax, %r8
 movq      %rdx, %r9

 movq (prime), %rax
 mulq  %r12
 movq      %rax, %r10
 movq      %rdx, %r13

 movq (prime+8), %rax
 mulq  %r11
 movq      %rax, %r14
 movq      %rdx, %r15

 movq (prime+8), %rax
 mulq  %r12

 addq %r10, %r9
 adcq %r13, %rax
 adcq   $0, %rdx
 addq %r14, %r9
 adcq %r15, %rax
 adcq   $0, %rdx

 // r = (c + u) div 2^(2*64) 
 movq $0, %r14
 addq   (reg_p2), %r8
 adcq  8(reg_p2), %r9
 adcq 16(reg_p2), %rax
 adcq 24(reg_p2), %rdx
 adcq   $0, %r14

 subq      (prime), %rax
 sbbq    (prime+8), %rdx
 sbbq           $0, %r14

 movq      %r14, %r11
 andq   (prime), %r14
 andq (prime+8), %r11

 addq      %r14, %rax
 adcq      %r11, %rdx

 movq      %rax,  (reg_p1)
 movq      %rdx, 8(reg_p1)

 pop %r15
 pop %r14
 pop %r13
 pop %r12

 ret



//**************************************************************
//  Addition in Fp
//  Operation: c [reg_p1] = (a [reg_p2] + b [reg_p3]) mod prime
//**************************************************************
.global ADD_Fp
.p2align 4
ADD_Fp:

 movq  (reg_p2), %rax
 movq 8(reg_p2), %rcx
 movq        $0, %r8

 addq  (reg_p3), %rax
 adcq 8(reg_p3), %rcx
 adcq        $0, %r8

 subq   (prime), %rax
 sbbq (prime+8), %rcx
 sbbq        $0, %r8
 
 movq       %r8, %r9
 andq   (prime), %r8
 andq (prime+8), %r9

 addq       %r8, %rax
 adcq       %r9, %rcx

 movq      %rax,  (reg_p1)
 movq      %rcx, 8(reg_p1)

 ret

//**************************************************************
//  Subtraction in Fp
//  Operation: c [reg_p1] = (a [reg_p2] - b [reg_p3]) mod prime
//**************************************************************
.global SUB_Fp
.p2align 4
SUB_Fp:

 movq  (reg_p2), %rax
 movq 8(reg_p2), %rcx
 movq        $0, %r8

 subq  (reg_p3), %rax
 sbbq 8(reg_p3), %rcx
 sbbq        $0, %r8

 movq       %r8, %r9
 andq   (prime), %r8
 andq (prime+8), %r9

 addq       %r8, %rax
 adcq       %r9, %rcx

 movq      %rax,  (reg_p1)
 movq      %rcx, 8(reg_p1)

 ret

//**************************************************************
//  Multiplication in Fp
//  Operation: c [reg_p1] = (a [reg_p2] * b [reg_p3]) mod prime
//**************************************************************
.global MUL_Fp
.p2align 4
MUL_Fp:

 push %r12
 push %r13
 push %r14
 push %r15
 push %rbx

 movq    reg_p3, %rcx

 // Integer multiplication
 // c = reg_p2*reg_p3
 movq  (reg_p2), %r8
 movq    (%rcx), %rax
 mulq  %r8
 movq      %rax, %r9
 movq      %rdx, %r10
 
 movq   8(%rcx), %rax
 mulq  %r8
 movq      %rax, %r11
 movq      %rdx, %r12

 movq 8(reg_p2), %r8
 movq    (%rcx), %rax
 mulq  %r8
 movq      %rax, %r13
 movq      %rdx, %r14
 
 movq   8(%rcx), %rax
 mulq  %r8
 movq      %rax, %r15
 movq      %rdx, %rbx

 addq %r11, %r10
 adcq %r12, %r15
 adcq   $0, %rbx
 addq %r13, %r10
 adcq %r14, %r15
 adcq   $0, %rbx

 // Montgomery reduction
 // m = mu*c mod 2^(64*2)
 movq (mu), %rax
 mulq  %r9
 movq      %rax, %r11
 movq      %rdx, %r12
 
 movq (mu), %rax
 mulq  %r10
 movq      %rax, %r13
 
 movq (mu+8), %rax
 mulq  %r9

 addq %r13, %r12
 addq %rax, %r12

 // u = m*p
 movq (prime), %rax
 mulq  %r11
 movq      %rax, (reg_p1)
 movq      %rdx, %r8

 movq (prime+8), %rax
 mulq  %r12
 movq  %rax, %r13
 movq  %rdx, %r14

 movq (prime), %rax
 mulq  %r12

 addq      %rax, %r8  
 adcq      %rdx, %r13
 adcq        $0, %r14
 
 movq (prime+8), %rax
 mulq  %r11

 addq      %rax, %r8  
 adcq      %rdx, %r13
 adcq        $0, %r14

 // r = (c + u) div 2^(2*64)
 movq   $0, %rax
 addq  %r9, (reg_p1)
 adcq %r10, %r8
 adcq %r15, %r13
 adcq %rbx, %r14
 adcq   $0, %rax

 subq      (prime), %r13
 sbbq    (prime+8), %r14
 sbbq           $0, %rax

 movq      %rax, %r11
 andq   (prime), %rax
 andq (prime+8), %r11

 addq      %rax, %r13
 adcq      %r11, %r14

 movq      %r13,  (reg_p1)
 movq      %r14, 8(reg_p1)

 pop %rbx
 pop %r15
 pop %r14
 pop %r13
 pop %r12

 ret

//**************************************************************
//  Squaring in Fp
//  Operation: c [reg_p1] = a [reg_p2]^2 mod prime
//**************************************************************
.global SQR_Fp
.p2align 4
SQR_Fp:

 push %r12
 push %r13
 push %r14
 push %r15
 push %rbx

 // Integer squaring
 // c = reg_p2*reg_p3
 movq  (reg_p2), %rax
 mulq  %rax
 movq      %rax, %r9
 movq      %rdx, %r10
 
 movq 8(reg_p2), %r8
 movq  (reg_p2), %rax
 mulq  %r8
 movq      %rax, %r13
 movq      %rdx, %r14
 
 movq  8(reg_p2), %rax
 mulq  %rax
 movq      %rax, %r15
 movq      %rdx, %rbx

 addq %r13, %r10
 adcq %r14, %r15
 adcq   $0, %rbx
 addq %r13, %r10
 adcq %r14, %r15
 adcq   $0, %rbx

 // Montgomery reduction
 // m = mu*c mod 2^(64*2)
 movq (mu), %rax
 mulq  %r9
 movq      %rax, %r11
 movq      %rdx, %r12
 
 movq (mu), %rax
 mulq  %r10
 movq      %rax, %r13
 
 movq (mu+8), %rax
 mulq  %r9

 addq %r13, %r12
 addq %rax, %r12

 // u = m*p
 movq (prime), %rax
 mulq  %r11
 movq      %rax, (reg_p1)
 movq      %rdx, %r8

 movq (prime+8), %rax
 mulq  %r12
 movq  %rax, %r13
 movq  %rdx, %r14

 movq (prime), %rax
 mulq  %r12

 addq      %rax, %r8  
 adcq      %rdx, %r13
 adcq        $0, %r14
 
 movq (prime+8), %rax
 mulq  %r11

 addq      %rax, %r8  
 adcq      %rdx, %r13
 adcq        $0, %r14

 // r = (c + u) div 2^(2*64)
 movq   $0, %rax
 addq  %r9, (reg_p1)
 adcq %r10, %r8
 adcq %r15, %r13
 adcq %rbx, %r14
 adcq   $0, %rax

 subq      (prime), %r13
 sbbq    (prime+8), %r14
 sbbq           $0, %rax

 movq      %rax, %r11
 andq   (prime), %rax
 andq (prime+8), %r11

 addq      %rax, %r13
 adcq      %r11, %r14

 movq      %r13,  (reg_p1)
 movq      %r14, 8(reg_p1)

 pop %rbx
 pop %r15
 pop %r14
 pop %r13
 pop %r12
 
 ret

