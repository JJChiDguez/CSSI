#ifndef _ARITH_H_
#define _ARITH_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <openssl/md5.h>

// 128-bit prime number
#include "setup_FF.h"

/*Integer Arithmetic functions for 4-word operands (c is the result)*/
extern void ADD_D(uint64_t *c, uint64_t *a, uint64_t *b);
extern void SHL_D(uint64_t *c, uint64_t *a);

/*Integer Arithmetic functions for 2-word operands (c is the result)*/
extern void SHR(uint64_t *c, uint64_t *a);
extern void SHL(uint64_t *c, uint64_t *a);
extern void ADD(uint64_t *c, uint64_t *a, uint64_t *b);
extern void SUB(uint64_t *c, uint64_t *a, uint64_t *b);
extern void MUL(uint64_t *c, uint64_t *a, uint64_t *b);
extern void SQR(uint64_t *c, uint64_t *a);
extern void REDC(uint64_t *c, uint64_t *a);

/*Finite Field Arithmetic functions for 2-word operands (c is the result)*/
extern void ADD_Fp(uint64_t *c, uint64_t *a, uint64_t *b);
extern void SUB_Fp(uint64_t *c, uint64_t *a, uint64_t *b);
extern void MUL_Fp(uint64_t *c, uint64_t *a, uint64_t *b);
extern void SQR_Fp(uint64_t *c, uint64_t *a);
void INV_Fp(uint64_t *b, uint64_t *a);

/*returns 0 if x==y, 1 if x>y and -1 fi x<y*/
static inline int compara(const uint64_t *x, const uint64_t *y, int NUM)
{
	int i;
	for (i=NUM-1; i >= 0; i--) 
    {
		if (x[i] != y[i]) 
            return x[i] > y[i] ? 1 : -1; 
	}
	return 0;
}

/*Generates a 2-words number a < p*/
static inline void random_NUM(uint64_t *a)
{
    int c;

    do
    {
        a[0] = rand();
        a[0] = a[0]<<32 | rand();
    
        a[1] = rand();
        a[1] = a[1]<<32 | rand();

        c = compara(a, p, 2);
    }
    while(c > 0);
}

/*Prints a number x of lenght NUM using the name c */
/* TIPO=1 prints the words of x with a space between them */
static inline void print_NUM(uint64_t *x, int NUM, int TIPO, char *c)
{
    int i;
    printf("%s := 0x", c);
    for(i=NUM-1; i > -1; i--){
        if(TIPO == 1)
            printf("%.16" PRIX64 " ", x[i]);
        else
            printf("%.16" PRIX64 "", x[i]);
    }
    printf(";\n");
}


/*Finite Field quadratic Arithmetic functions for 4-word operands (c is the result)*/

static inline void print_num(uint64_t *x, int NUM)
{
    int i;
    printf("0x");
    for(i=NUM-1; i > -1; i--)
            printf("%.16" PRIX64 "", x[i]);
}

static inline void print_hex(char output_string[WORD_N * 2 * 16 + 1], uint64_t *x, int NUM)
{
    int i;
    unsigned char v_string[16 + 1];
    output_string[0] = '\0';
    for(i=NUM-1; i > -1; i--)
    {
       v_string[0] = '\0';
		 sprintf(v_string, "%.16" PRIX64 "", x[i]);
		 v_string[strlen(v_string)] = '\0';
		 strcat(output_string, v_string);
		 output_string[strlen(output_string)] = '\0';
    }
}

static inline void print_ELE(uint64_t (*x)[2], int NUM,  char *c)
{
	printf("%s := ", c);
	print_num( x[0], NUM);
	printf(" + ");
	print_num( x[1], NUM);
	printf(" * i;\n");
}

int compare_in_Fp2(uint64_t local_aux_00[2][WORD_N], uint64_t local_aux_01[2][WORD_N]);

void ADD_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N], uint64_t b[2][WORD_N]);
void SUB_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N], uint64_t b[2][WORD_N]);
void MUL_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N], uint64_t b[2][WORD_N]);
void SQR_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N]);
void INV_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N]);

void   NEG_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N]);
void   TIMES_i(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N]);
void CONJUGATE(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N]);

void  EXP_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N], uint64_t *n);
void SQRT_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N]);

void MOD_MAPSTO_MONT(uint64_t a[2][WORD_N]);
void MONT_MAPSTO_MOD(uint64_t a[2][WORD_N]);


//Measuring the perfomance
static uint64_t get_cycles()
{
   uint32_t lo, hi;
   asm volatile("rdtsc":"=a"(lo),"=d"(hi));
   return ((uint64_t)hi<<32) | lo;
};

#endif
