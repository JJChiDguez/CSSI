#include "FieldArith.h"

// ---------------------------------------------------------------
// compare_in_Fp2() function :::::::::::::::::::::::::::::::::::::
// inputs:
//        local_aux_00 : an element of Fp2,
//        local_aux_01 : an element of Fp2
// output:
//         0  if local_aux_00 == local_aux_01,
//        -1  if local_aux_00 < local_aux_01, or
//         1  if local_aux_00 > local_aux_01,
int compare_in_Fp2(uint64_t local_aux_00[2][WORD_N], uint64_t local_aux_01[2][WORD_N])
{
	int local_1st = compara(local_aux_00[0], local_aux_01[0], WORD_N);
	int local_2nd = compara(local_aux_00[1], local_aux_01[1], WORD_N);

	if(local_1st != 0)
		return local_1st;
	else
		return local_2nd;
};


// ---------------------------------------------------------------
// INV_Fp() function :::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp,
// output:
//        b : a^-1
void INV_Fp(uint64_t *b, uint64_t *a)
{
    uint64_t x1[3]={0}, x2[3]={0}, x[5]={0};
    uint64_t u[3], v[3];
    int k=0, aux;

    memcpy(u, a, sizeof(uint64_t)*2);
    memcpy(v, p, sizeof(uint64_t)*2);
    x1[0] = 1;
    

    while(v[0] != 0 || v[1] != 0)
    {
        if(!(v[0]&0x1))
        {
            SHR(v, v);
            SHL(x1, x1);
        }
        else if(!(u[0]&0x1))
        {
            SHR(u, u);
            SHL(x2, x2);
        }
        else if(compara(v, u, 2) >= 0)
        {
            SUB(v, v, u);
            SHR(v, v);
            ADD(x2, x2, x1);
            SHL(x1, x1);
            
            if(x1[2] == 1)
                SUB(x1, x1, p);
        }
        else
        {
            SUB(u, u, v);
            SHR(u, u);
            ADD(x1, x1, x2);
            SHL(x2, x2);
        }
        k++;
    }

    if (k < 128)
    {
        MUL_Fp(x1, x1, base_mont_2);
        k = k + 128;
    }
    MUL_Fp(x1, x1, base_mont_2);
    
    aux = (256-k) / 64;
    k = (256-k)%64;
    
    memcpy(&x[aux], x1, sizeof(uint64_t)*2);

    while(k>0)
    {
        SHL_D(x, x);
        k--;
    }
    
    REDC(b, x);
   
}



// ---------------------------------------------------------------
// ADD_Fp2() function ::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
//        b : an element of Fp2,
// output:
//        c : a + b
void ADD_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N], uint64_t b[2][WORD_N])
{
	ADD_Fp(c[0], a[0], b[0]);
	ADD_Fp(c[1], a[1], b[1]);
}	// 2 ADDS in Fp

// ---------------------------------------------------------------
// SUB_Fp2() function ::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
//        b : an element of Fp2,
// output:
//        c : a - b
void SUB_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N], uint64_t b[2][WORD_N])
{
	SUB_Fp(c[0], a[0], b[0]);
	SUB_Fp(c[1], a[1], b[1]);
}	// 2 ADDS (SUBS) in Fp

// ---------------------------------------------------------------
// MUL_Fp2() function ::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
//        b : an element of Fp2,
// output:
//        c : a * b
void MUL_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N], uint64_t b[2][WORD_N])
{
    uint64_t z0[WORD_N], z1[WORD_N], z2[WORD_N], z3[WORD_N];

    ADD_Fp(z0, a[0], a[1]);	// a[0] + a[1]
    ADD_Fp(z1, b[0], b[1]);	// b[0] + b[1]
    
    MUL_Fp(c[1], z0, z1);	// (a[0] + a[1]) * (b[0] + b[1])
    MUL_Fp(z2, a[0], b[0]);	// a[0] * b[0]
    MUL_Fp(z3, a[1], b[1]);	// a[1] * b[1]
    
    SUB_Fp(c[0], z2, z3);	//  a[0] * b[0] -  a[1] * b[1]
    
    SUB_Fp(c[1], c[1], z2);	//  (a[0] + a[1]) * (b[0] + b[1]) - a[0] * b[0]
    SUB_Fp(c[1], c[1], z3); //  (a[0] + a[1]) * (b[0] + b[1]) - a[0] * b[0] - a[1] * b[1] = a[1] * b[0] + a[0] * b[1]
}	// 3 MULS + 5 ADDS in Fp

// ---------------------------------------------------------------
// SQR_Fp2() function ::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
// output:
//        c : a^2
void SQR_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N])
{
    uint64_t z0[WORD_N], z1[WORD_N], z2[WORD_N], z3[WORD_N];
    
    ADD_Fp(z0, a[0], a[0]);	// 2 * a[0]
    ADD_Fp(z1, a[0], a[1]);	// a[0] + a[1]
    SUB_Fp(z2, a[0], a[1]);	// a[0] - a[1]
    
    MUL_Fp(c[0], z1, z2);	// (a[0] + a[1]) * (a[0] - a[1]) = a[0]^2 - a[1]^2
    MUL_Fp(c[1], z0, a[1]);	// 2 * a[0] * a[1]
}	// 2 MULS + 3 ADDS in Fp

// ---------------------------------------------------------------
// INV_Fp2() function ::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
// output:
//        c : a^-1
void INV_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N])
{
    uint64_t N0[WORD_N], N1[WORD_N], S[WORD_N], S1[WORD_N], S2[WORD_N], AUX[WORD_N];
    
    SQR_Fp(N0, a[0]);	// a[0] ^ 2
    SQR_Fp(N1, a[1]);	// a[1] ^ 2
    
    ADD_Fp(S, N0, N1);	// a[0] ^ 2 + a[1] ^ 2 = Norm(a[0] + a[1] * i)
    
    INV_Fp(S1, S);		// 1 / (a[0] ^ 2 + a[1] ^ 2)
    
    SUB_Fp(S2, zeroM, a[1]);	// -a[1]
    
    MUL_Fp(AUX, S1, a[0]);	//  a[0] / (a[0] ^ 2 + a[1] ^ 2)
    memcpy(c[0], AUX, sizeof(uint64_t) * WORD_N);
    MUL_Fp(c[1], S1, S2);   // -a[1] / (a[0] ^ 2 + a[1] ^ 2)
}	// 1 INV + 2 MULS + 2 SQR + 2 ADDS in Fp

// ---------------------------------------------------------------
// EXP_Fp2() function ::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
//        n : an integer less than p - 1
// output:
//        c : a^n
void EXP_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N], uint64_t *n)
{
    uint64_t y[2][WORD_N], x[2][WORD_N], z[2][WORD_N];
    memcpy(y[0], unoM, sizeof(uint64_t) * WORD_N);
    memcpy(y[1], zeroM, sizeof(uint64_t) * WORD_N); 
    memcpy(x, a, sizeof(uint64_t) * WORD_N * 2);
   
    uint64_t v[WORD_N];
    memcpy(v, n, sizeof(uint64_t) * WORD_N);

    while( compara(v, uno, WORD_N) != 0)	// because v > 0: we have v != 1 => v > 1
    {
        if (!(v[0]&0x0000000000000001)){
            SQR_Fp2(z, x);
            memcpy(x, z, sizeof(uint64_t) * WORD_N * 2);
        }
        else{
            MUL_Fp2(z, y, x);
            memcpy(y, z, sizeof(uint64_t) * WORD_N * 2);
            SQR_Fp2(z, x);
            memcpy(x, z, sizeof(uint64_t) * WORD_N * 2);
        }
        SHR(v, v);
    }
    
    MUL_Fp2(c, y, x);
}

// ---------------------------------------------------------------
// NEG_Fp2() function ::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
// output:
//        c : -a
void NEG_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N])
{
    SUB_Fp(c[0], zeroM, a[0]);
    SUB_Fp(c[1], zeroM, a[1]);
}

// ---------------------------------------------------------------
// CONJUGATE() function ::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
// output:
//        c : the conjugate of a, i.e., if a = a_0 + a_1 * i, then
//            c = a_0 - a_1 * i
void CONJUGATE(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N])
{
    memcpy(c[0], a[0], sizeof(uint64_t) * WORD_N);
    SUB_Fp(c[1], zeroM, a[1]);
}

// ---------------------------------------------------------------
// TIMES_i() function ::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
// output:
//        c : i * a
void TIMES_i(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N])
{
    memcpy(c[1], a[0], sizeof(uint64_t) * WORD_N);
    SUB_Fp(c[0], zeroM, a[1]);
}

// ---------------------------------------------------------------
// MOD_MAPSTO_MONT() function ::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
// output:
//        c : a * R
void MOD_MAPSTO_MONT(uint64_t a[2][WORD_N])
{
   MUL_Fp(a[0], a[0], base_mont_2);
   MUL_Fp(a[1], a[1], base_mont_2);
}

// ---------------------------------------------------------------
// MONT_MAPSTO_MOD() function ::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
// output:
//        c : a * (R^-1)
void MONT_MAPSTO_MOD(uint64_t a[2][WORD_N])
{
   MUL_Fp(a[0], a[0], uno);
   MUL_Fp(a[1], a[1], uno);
}

// ---------------------------------------------------------------
// SQRT_Fp2() function :::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        a : an element of Fp2,
// output:
//        c : a^(1/2)
void SQRT_Fp2(uint64_t c[2][WORD_N], uint64_t a[2][WORD_N])
{
    uint64_t a1[2][WORD_N],  alpha[2][WORD_N], a0[2][WORD_N], x0[2][WORD_N];
    EXP_Fp2(a1, a, constant_sqrt1);	//    a1 <- a^( [p-3]/4)
    SQR_Fp2(a0, a1);				//    a0 <- a1^2
    MUL_Fp2(alpha, a0, a);			// alpha <- a1^2 * a
    MUL_Fp2(x0, a1, a);			//    x0 <- a1 * a
    
    if( (compara(alpha[1], zeroM, WORD_N) == 0) && (compara(alpha[0], negM, WORD_N) == 0) )
        TIMES_i(c, x0);
    else{
        ADD_Fp(alpha[0], alpha[0], unoM); 	// alpha <- alpha + 1
        EXP_Fp2(a0, alpha, constant_sqrt2);	//    a0 <- (alpha + 1)^( [p - 1] / 2)
        MUL_Fp2(c, a0, x0);				//     c <- a0 * x0
    }
}
