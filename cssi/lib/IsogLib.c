#include "IsogLib.h"

// ---------------------------------------------------------------
// get_isomorphism() function ::::::::::::::::::::::::::::::::::::
// inputs:
//        E0 : a supersingular elliptic curve defined over Fp2,
//        E1 : a supersingular elliptic curve defined over Fp2,
// output:
//        I : a 1-isogeny (i.e., isomorphism) between E0 and E1
void get_isomorphism(Isomorphism *I, Curve E0, Curve E1){
	// GENERAL FORM: E0 and E1 have their parameters different from zero
	uint64_t  AUX0[2][WORD_N], AUX1[2][WORD_N], AUX2[2][WORD_N];

	MUL_Fp2(AUX0, E0.A, E1.B);
	MUL_Fp2(AUX1, E0.B, E1.A);
	INV_Fp2(AUX2, AUX0);
	MUL_Fp2(AUX0, AUX2, AUX1);
	SQRT_Fp2(I->Z, AUX0);
}	// 1*SQRT + 1*INV + 3*MUL

// ---------------------------------------------------------------
// eval_isomorphism() function :::::::::::::::::::::::::::::::::::
// inputs:
//        I : an isomorphism at E,
//        P : a point of the supersingular elliptic curve E,
// output:
//        Q : I(P)
void eval_isomorphism(Point *Q, Isomorphism I, Point P){
	memcpy( Q->X, P.X, sizeof(uint64_t) * WORD_N * 2);
	memcpy( Q->Y, P.Y, sizeof(uint64_t) * WORD_N * 2);
	MUL_Fp2(Q->Z, P.Z, I.Z);
}	// 1*MUL


// ---------------------------------------------------------------
// get_2_isog() function :::::::::::::::::::::::::::::::::::::::::
// inputs:
//        E0 : a supersingular elliptic curve defined over Fp2,
//         K : an order-2 point,
// output:
//        f : a degree-2 isogeny from E0 to E0/<K>,
//       EA : the supersingular elliptic curve E0/<K>, which is 
//            2-isogenous to E0
void get_2_isog(Curve *EA, Isogeny *f, Point K, Curve E0){
	memcpy(f->Kx, K.X, sizeof(uint64_t) * WORD_N * 2);	//  K.X
	SQR_Fp2(f->KzKz, K.Z);										// (K.Z)^2
	MUL_Fp2(f->KzKzKz, f->KzKz, K.Z);						// (K.Z)^3

	uint64_t AUX[2][WORD_N], FLAG[2][WORD_N], KXKX[2][WORD_N];
	
	SQR_Fp2(AUX, f->KzKz);										// (K.Z)^4
	MUL_Fp2(FLAG, AUX, E0.A);									// (K.Z)^4 * E0.A

	SQR_Fp2(KXKX, K.X);											// (K.X)^2
	ADD_Fp2(AUX, KXKX, KXKX);									// 2 * (K.X)^2
	ADD_Fp2(AUX, AUX, KXKX);									// 3 * (K.X)^2

	ADD_Fp2(f->Gxza, AUX, FLAG);								// 3 * (K.X)^2 + (K.Z)^4 * E0.A

	ADD_Fp2(AUX, f->Gxza, f->Gxza);							// 2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A)
	ADD_Fp2(AUX, AUX, AUX);										// 4 * (3 * (K.X)^2 + (K.Z)^4 * E0.A)
	ADD_Fp2(EA->A, AUX, f->Gxza);								// 5 * (3 * (K.X)^2 + (K.Z)^4 * E0.A)
	SUB_Fp2(EA->A, FLAG, EA->A);								// (K.Z)^4 * E0.A - 5 * (3 * (K.X)^2 + (K.Z)^4 * E0.A)

	ADD_Fp2(AUX, AUX, AUX);										// 8 * (3 * (K.X)^2 + (K.Z)^4 * E0.A)
	SUB_Fp2(AUX, AUX, f->Gxza);								// 7 * (3 * (K.X)^2 + (K.Z)^4 * E0.A)
	MUL_Fp2(EA->B, AUX, K.X);									// 7 * K.X * (3 * (K.X)^2 + (K.Z)^4 * E0.A)
	SQR_Fp2(AUX, f->KzKzKz);									// (K.Z)^6
	MUL_Fp2(FLAG, AUX, E0.B);									// (K.Z)^6 * E0.B
	SUB_Fp2(EA->B, FLAG, EA->B);								// (K.Z)^6 * E0.B - 7 * K.X * (3 * (K.X)^2 + (K.Z)^4 * E0.A)
}	// 4*MUL + 4*SQR + 7*ADD + 3*SUB

// ---------------------------------------------------------------
// eval_2_isog() function :::::::::::::::::::::::::::::::::::
// inputs:
//        f : a degree-2 isogeny at E,
//        P0 : a point of the supersingular elliptic curve E,
// output:
//        PA : f(P)
void eval_2_isog(Point *PA, Isogeny f, Point P0){
	uint64_t PZPZ[2][WORD_N], PZPZPZPZ[2][WORD_N], 
				PZPZKX[2][WORD_N], KZKZPX[2][WORD_N], 
				SUB_AUX[2][WORD_N], KZKZKZ[2][WORD_N], 
				AUX_1[2][WORD_N], AUX_2[2][WORD_N], AUX_3[2][WORD_N];

	SQR_Fp2(PZPZ, P0.Z);							// (P0.Z)^2
	SQR_Fp2(PZPZPZPZ, PZPZ);					// (P0.Z)^4
	MUL_Fp2(PZPZKX, PZPZ, f.Kx);				// (P0.Z)^2 * f.Kx
	MUL_Fp2(KZKZPX, f.KzKz, P0.X);			// (f.Kz)^2 * P0.X
	SUB_Fp2(SUB_AUX, KZKZPX, PZPZKX);		// (f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx
	MUL_Fp2(AUX_1, KZKZPX, SUB_AUX);			// [(f.Kz)^2 * P0.X] * [(f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx]
	MUL_Fp2(AUX_2, f.Gxza, PZPZPZPZ);		// f.Gxza * (P0.Z)^4
	ADD_Fp2(AUX_1, AUX_1, AUX_2);				// [(f.Kz)^2 * P0.X] * [(f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx] + f.Gxza * (P0.Z)^4

	MUL_Fp2(PA->X, SUB_AUX, AUX_1);			// { (f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx } * 
														// { [(f.Kz)^2 * P0.X] * [(f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx] + f.Gxza * (P0.Z)^4 }

	SQR_Fp2(AUX_1, SUB_AUX);					// [(f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx]^2
	SUB_Fp2(AUX_1, AUX_1, AUX_2);				// [(f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx]^2 - f.Gxza * (P0.Z)^4
	MUL_Fp2(AUX_3, AUX_1, SUB_AUX);			//	{ (f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx } *
														// { [(f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx]^2 - f.Gxza * (P0.Z)^4 }
	MUL_Fp2(AUX_2, AUX_3, f.KzKzKz);			//	(f.Kz)^3 * { (f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx } *
														// { [(f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx]^2 - f.Gxza * (P0.Z)^4 }	
	MUL_Fp2(PA->Y, P0.Y, AUX_2);				//	P0.Y * (f.Kz)^3 * { (f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx } *
														// { [(f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx]^2 - f.Gxza * (P0.Z)^4 }	

	MUL_Fp2(PA->Z, P0.Z, SUB_AUX);			// P0.Z * { (f.Kz)^2 * P0.X - (P0.Z)^2 * f.Kx }
}	// 9*MUL + 3*SQR + 1*ADD + 2*SUB

// ---------------------------------------------------------------
// get_3_isog() function :::::::::::::::::::::::::::::::::::::::::
// inputs:
//        E0 : a supersingular elliptic curve defined over Fp2,
//         K : an order-3 point,
// output:
//        f : a degree-3 isogeny from E0 to E0/<K>
//       EA : the supersingular elliptic curve E0/<K>, which is 
//            3-isogenous to E0
void get_3_isog(Curve *EA, Isogeny *f, Point K, Curve E0){
	memcpy(f->Kx, K.X, sizeof(uint64_t) * WORD_N * 2);	//  K.X
	memcpy(f->Ky, K.Y, sizeof(uint64_t) * WORD_N * 2);	//  K.Y
	SQR_Fp2(f->KyKy4, K.Y);										// (K.Y)^2
	ADD_Fp2(f->KyKy4, f->KyKy4, f->KyKy4);					// 2 * (K.Y)^2
	ADD_Fp2(f->KyKy4, f->KyKy4, f->KyKy4);					// 4 * (K.Y)^2
	SQR_Fp2(f->KzKz, K.Z);										// (K.Z)^2
	MUL_Fp2(f->KzKzKz, f->KzKz, K.Z);						// (K.Z)^3

	uint64_t AUX[2][WORD_N], KXKX[2][WORD_N], FLAG[2][WORD_N];
	
	SQR_Fp2(AUX, f->KzKz);										// (K.Z)^4
	MUL_Fp2(FLAG, AUX, E0.A);									// (K.Z)^4 * E0.A

	SQR_Fp2(KXKX, K.X);											// (K.X)^2
	ADD_Fp2(AUX, KXKX, KXKX);									// 2 * (K.X)^2
	ADD_Fp2(AUX, AUX, KXKX);									// 3 * (K.X)^2

	ADD_Fp2(f->Gxza, AUX, FLAG);								// 3 * (K.X)^2 + (K.Z)^4 * E0.A
	ADD_Fp2(f->Gxza, f->Gxza, f->Gxza);						// 2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A)

	ADD_Fp2(AUX, f->Gxza, f->Gxza);							// 2 * [2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A) ]
	ADD_Fp2(AUX, AUX, AUX);										// 4 * [2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A) ]
	ADD_Fp2(EA->A, AUX, f->Gxza);								// 5 * [2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A) ]
	SUB_Fp2(EA->A, FLAG, EA->A);								// E0.A * (K.Z)^4 - 5 * [2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A) ]

	MUL_Fp2(AUX, K.X, f->Gxza);								// K.X * [2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A) ]
	ADD_Fp2(AUX, f->KyKy4, AUX);									// 4*(K.Y)^2 +  K.X * [2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A) ]
	ADD_Fp2(EA->B, AUX, AUX);									// 2 * { 4*(K.Y)^2 +  K.X * [2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A) ] }
	ADD_Fp2(EA->B, EA->B, EA->B);								// 4 * { 4*(K.Y)^2 +  K.X * [2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A) ] }
	ADD_Fp2(EA->B, EA->B, EA->B);								// 8 * { 4*(K.Y)^2 +  K.X * [2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A) ] }
	SUB_Fp2(EA->B, EA->B, AUX);								// 7 * { 4*(K.Y)^2 +  K.X * [2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A) ] }
	SQR_Fp2(AUX, f->KzKzKz);									// (K.Z)^6
	MUL_Fp2(FLAG, AUX, E0.B);									// (K.Z)^6 * E0.B
	SUB_Fp2(EA->B, FLAG, EA->B);								// (K.Z)^6 * E0.B - 7 * { 4*(K.Y)^2 +  K.X * [2 * (3 * (K.X)^2 + (K.Z)^4 * E0.A) ] }
}	// 4*MUL + 5*SQR + 13*ADD + 3*SUB

// ---------------------------------------------------------------
// eval_3_isog() function :::::::::::::::::::::::::::::::::::
// inputs:
//        f : a degree-3 isogeny at E,
//        P : a point of the supersingular elliptic curve E,
// output:
//        Q : f(P)
void eval_3_isog(Point *Q, Isogeny f, Point P){
	uint64_t PZPZ[2][WORD_N], PZPZPZ[2][WORD_N], PZPZPZPZ[2][WORD_N], PZPZPZPZPZPZ[2][WORD_N],
				PZPZKX[2][WORD_N], KZKZPX[2][WORD_N], 
				PZPZPZKY[2][WORD_N], KZKZKZPY[2][WORD_N], 
				SUB1_AUX[2][WORD_N], SUB2_AUX[2][WORD_N],
				KZKZKZ[2][WORD_N], 
				AUX1[2][WORD_N], AUX2[2][WORD_N], AUX3[2][WORD_N], 
				AUX4[2][WORD_N], AUX5[2][WORD_N], AUX6[2][WORD_N];
				
	SQR_Fp2(PZPZ, P.Z);								// (P.Z)^2
	MUL_Fp2(PZPZPZ, PZPZ, P.Z);					// (P.Z)^3
	SQR_Fp2(PZPZPZPZ, PZPZ);						// (P.Z)^4
	SQR_Fp2(PZPZPZPZPZPZ, PZPZPZ);				// (P.Z)^6

	MUL_Fp2(PZPZKX, PZPZ, f.Kx);					// (P.Z)^2 * f.Kx
	MUL_Fp2(KZKZPX, f.KzKz, P.X);					// (f.Kz)^2 * P.X
	SUB_Fp2(SUB1_AUX, KZKZPX, PZPZKX);			// (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx

	MUL_Fp2(Q->Z, P.Z, SUB1_AUX);					// P.Z * [ (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx ]

	MUL_Fp2(PZPZPZKY, PZPZPZ, f.Ky);				// (P.Z)^3 * f.Ky
	MUL_Fp2(KZKZKZPY, f.KzKzKz, P.Y);			// (f.Kz)^3 * P.Y
	SUB_Fp2(SUB2_AUX, KZKZKZPY, PZPZPZKY);		// (f.Kz)^3 * P.Y - (P.Z)^3 * f.Ky

	MUL_Fp2(AUX1, PZPZPZPZ, SUB1_AUX);			// (P.Z)^4 * [ (f.Kz)^2 * P.X - (P.Z)^2 * f.Kx ]
	MUL_Fp2(AUX6, f.Gxza, AUX1);					// f.Gxza * (P.Z)^4 * [ (f.Kz)^2 * P.X - (P.Z)^2 * f.Kx ]
	MUL_Fp2(AUX3, f.KyKy4, PZPZPZPZPZPZ);		// f.KyKy4 * (P.Z)^6
	ADD_Fp2(AUX2, AUX6, AUX3);						// f.Gxza * (P.Z)^4 * [ (f.Kz)^2 * P.X - (P.Z)^2 * f.Kx ] + f.KyKy4 * (P.Z)^6
	SQR_Fp2(AUX4, SUB1_AUX);						// [ (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx ]^2
	MUL_Fp2(AUX5, f.KzKz, AUX4);					// f.KzKz * { [ (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx ]^2 }
	MUL_Fp2(Q->X, P.X, AUX5);						// P.X * f.KzKz * { [ (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx ]^2 }
	ADD_Fp2(Q->X, Q->X, AUX2);						// P.X * f.KzKz * { [ (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx ]^2 } + 
															// f.Gxza * (P.Z)^4 * [ (f.Kz)^2 * P.X - (P.Z)^2 * f.Kx ] + f.KyKy4 * (P.Z)^6

	MUL_Fp2(AUX5, SUB1_AUX, AUX4);				// [ (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx ]^3
	MUL_Fp2(AUX4, f.KzKzKz, AUX5);				// f.KzKzKz * [ (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx ]^3
	MUL_Fp2(Q->Y, P.Y, AUX4);						// P.Y * f.KzKzKz * [ (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx ]^3	
	MUL_Fp2(AUX5, AUX6, SUB2_AUX);				// f.Gxza * (P.Z)^4 * [ (f.Kz)^2 * P.X - (P.Z)^2 * f.Kx ] * 
															// [ (f.Kz)^3 * P.Y - (P.Z)^3 * f.Ky ]
	SUB_Fp2(Q->Y, Q->Y, AUX5);						// P.Y * f.KzKzKz * [ (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx ]^3	-
															// f.Gxza * (P.Z)^4 * [ (f.Kz)^2 * P.X - (P.Z)^2 * f.Kx ] * 
															// [ (f.Kz)^3 * P.Y - (P.Z)^3 * f.Ky ]
	MUL_Fp2(AUX5, P.Y, f.KzKzKz);					// P.Y * f.KzKzKz
	MUL_Fp2(AUX4, AUX3, AUX5);						// P.Y * f.KzKzKz * f.KyKy4 * (P.Z)^6
	ADD_Fp2(AUX4, AUX4, AUX4);						// 2 * [P.Y * f.KzKzKz * f.KyKy4 * (P.Z)^6 ]
	SUB_Fp2(Q->Y, Q->Y, AUX4);						// P.Y * f.KzKzKz * [ (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx ]^3	-
															// f.Gxza * (P.Z)^4 * [ (f.Kz)^2 * P.X - (P.Z)^2 * f.Kx ] * 
															// [ (f.Kz)^3 * P.Y - (P.Z)^3 * f.Ky ] - 
															// 2 * [ P.Y * f.KzKzKz * f.KyKy4 * (P.Z)^6 ]
	MUL_Fp2(AUX4, f.Ky, PZPZPZ);					// f.Ky * (P.Z)^3
	MUL_Fp2(AUX5, AUX6, AUX4);						// f.Gxza * f.Ky * (P.Z)^7 * [ (f.Kz)^2 * P.X - (P.Z)^2 * f.Kx ]
	SUB_Fp2(Q->Y, Q->Y, AUX5);						// P.Y * f.KzKzKz * [ (f.Kz)^2 * P.X - (P0.Z)^2 * f.Kx ]^3	-
															// f.Gxza * (P.Z)^4 * [ (f.Kz)^2 * P.X - (P.Z)^2 * f.Kx ] * 
															// [ (f.Kz)^3 * P.Y - (P.Z)^3 * f.Ky ] - 
															// 2 * [ P.Y * f.KzKzKz * f.KyKy4 * (P.Z)^6 ] + 
															// f.Gxza * f.Ky * (P.Z)^7 * [ (f.Kz)^2 * P.X - (P.Z)^2 * f.Kx ]
}	// 19*MUL + 4*SQR + 3*ADD + 4*SUB

// ---------------------------------------------------------------
// spliting_scalars_Floor() function :::::::::::::::::::::::::::::
// inputs:
//   SPLITS : a sequence of points P_i such that if P_0 has order d^t, then 
//            P_i has order d^floor(t/2^i) for each i>= 1. Here the 
//            last element doesn't has order d.
//     ORDS : the sequence of floor(t/2^i) for each i>= 0,
//   lenght : the current element position to be added,
//  ocal_SM : the morphism [d],
//        a : the coefficient a of a supersingular elliptic curve E
// output:
//   SPLITS : a sequence of points P_i such that if P_0 has order d^t, then 
//            P_i has order d^floor(t/2^i) for each i>= 0. Here the 
//            last element has order d.
void  spliting_scalars_Floor(Point SPLITS[Log2_E], uint64_t ORDS[Log2_E], uint64_t *lenght, void (*local_SM)(), uint64_t a[2][WORD_N]){
	uint64_t aux_e = ORDS[*lenght], acc_e = aux_e;
	while(aux_e >= 2){
		local_SM(&(SPLITS[(*lenght) + 1]), SPLITS[(*lenght)],  aux_e / 2, a);
		acc_e -= (aux_e / 2);
		aux_e = aux_e - (aux_e / 2);
		ORDS[(*lenght) + 1] = acc_e;
		(*lenght) += 1;
	};
}

// ---------------------------------------------------------------
// spliting_scalars_Ceil() function ::::::::::::::::::::::::::::::
// inputs:
//   SPLITS : a sequence of points P_i such that if P_0 has order d^t, then 
//            P_i has order d^ceil(t/2^i) for each i>= 1. Here the 
//            last element doesn't has order d.
//     ORDS : the sequence of ceil(t/2^i) for each i>= 0,
//   lenght : the current element position to be added,
// local_SM : the morphism [d],
//        a : the coefficient a of a supersingular elliptic curve E
// output:
//   SPLITS : a sequence of points P_i such that if P_0 has order d^t, then 
//            P_i has order d^ceil(t/2^i) for each i>= 0. Here the 
//            last element has order d.
void spliting_scalars_Ceil(Point SPLITS[Log2_E], uint64_t ORDS[Log2_E], uint64_t *lenght, void (*local_SM)(), uint64_t a[2][WORD_N]){
	uint64_t aux_e = ORDS[*lenght], acc_e = aux_e;
	while(aux_e >= 2){
		local_SM(&(SPLITS[(*lenght) + 1]), SPLITS[(*lenght)],  aux_e - (aux_e / 2), a);
		acc_e -= (aux_e - (aux_e / 2));
		aux_e = aux_e / 2;
		ORDS[(*lenght) + 1] = acc_e;
		(*lenght) += 1;
	};
}

// ---------------------------------------------------------------
// get_isogenous_curve() function ::::::::::::::::::::::::::::::::
// inputs:
//          K : an order-(d^e) point,
//   local_SM : the morphism [d],
//  local_get : a degree-d isogeny computation function,
// local_eval : a degree-d isogeny evaluation function
//         E0 : a supersingular elliptic curve defined over Fp2,
//      level : the integer e
// output:
//        E_k : the supersingular elliptic curve E0/<R>, which is 
//              (d^e)-isogenous to E0
void get_isogenous_curve(Curve *E_k, Point R, void (*local_SM)(), void (*local_get)(), void (*local_eval)(), Curve E0, uint64_t level){	
	// setup
	memcpy( &(E_k->A[0]), E0.A[0], sizeof(uint64_t) * WORD_N);
	memcpy( &(E_k->A[1]), E0.A[1], sizeof(uint64_t) * WORD_N);
	memcpy( &(E_k->B[0]), E0.B[0], sizeof(uint64_t) * WORD_N);
	memcpy( &(E_k->B[1]), E0.B[1], sizeof(uint64_t) * WORD_N);

	Point local_SPLITS[Log2_E];
	Isogeny phi_k;
	uint64_t local_ORDS[Log2_E], local_i, local_k, local_j;

	// initialization
	Point_Assign(&local_SPLITS[0], R);
	local_ORDS[0] = level;
	local_i = 0;
	local_k = 0;

	spliting_scalars_Floor(local_SPLITS, local_ORDS, &local_i, local_SM, E_k->A);

	while( local_k < (level - (level/2) ) ){
		if(local_ORDS[local_i] != 1)
			spliting_scalars_Floor(local_SPLITS, local_ORDS, &local_i, local_SM, E_k->A);
		else{
			local_get(E_k, &phi_k, local_SPLITS[local_i], *E_k);
			for(local_j = 0; local_j < local_i; local_j++){
				local_eval(&local_SPLITS[local_j], phi_k, local_SPLITS[local_j]);
				local_ORDS[local_j] -= 1;
			}
			local_k += 1;
			local_i -= 1;
		}
	}

	while( local_k < level ){
		if(local_ORDS[local_i] != 1)
			spliting_scalars_Ceil(local_SPLITS, local_ORDS, &local_i, local_SM, E_k->A);
		else{
			local_get(E_k, &phi_k, local_SPLITS[local_i], *E_k);
			for(local_j = 0; local_j < local_i; local_j++){
				local_eval(&local_SPLITS[local_j], phi_k, local_SPLITS[local_j]);
				local_ORDS[local_j] -= 1;
			}
			local_k += 1;
			local_i -= 1;
		}
	}	
}

// ---------------------------------------------------------------
// get_isogenous_curve() function ::::::::::::::::::::::::::::::::
// inputs:
//         W0 : a point in the supersingular elliptic curve,
//         W1 : a point in the supersingular elliptic curve,
//          K : an order-(d^e) point,
//   local_SM : the morphism [d],
//  local_get : a degree-d isogeny computation function,
// local_eval : a degree-d isogeny evaluation function
//         E0 : a supersingular elliptic curve defined over Fp2,
//      level : the integer e
// output:
//         W0 : the image of the initial point W0 [under the degree-(d^e) isogeny],
//         W1 : the image of the initial point W1 [under the degree-(d^e) isogeny],
void eval_isogeny(Point *W0, Point *W1, Point R, void (*local_SM)(), void (*local_get)(), void (*local_eval)(), Curve E0, uint64_t level){
	// setup
	Curve E_k;
	memcpy( &(E_k.A[0]), E0.A[0], sizeof(uint64_t) * WORD_N);
	memcpy( &(E_k.A[1]), E0.A[1], sizeof(uint64_t) * WORD_N);
	memcpy( &(E_k.B[0]), E0.B[0], sizeof(uint64_t) * WORD_N);
	memcpy( &(E_k.B[1]), E0.B[1], sizeof(uint64_t) * WORD_N);

	Point local_SPLITS[Log2_E];
	Isogeny phi_k;
	uint64_t local_ORDS[Log2_E], local_i, local_k, local_j;

	// initialization
	Point_Assign(&local_SPLITS[0], R);
	local_ORDS[0] = level;
	local_i = 0;
	local_k = 0;

	spliting_scalars_Floor(local_SPLITS, local_ORDS, &local_i, local_SM, E_k.A);

	while( local_k < (level - (level/2) ) ){
		if(local_ORDS[local_i] != 1)
			spliting_scalars_Floor(local_SPLITS, local_ORDS, &local_i, local_SM, E_k.A);
		else{
			local_get(&E_k, &phi_k, local_SPLITS[local_i], E_k);
			for(local_j = 0; local_j < local_i; local_j++){
				local_eval(&local_SPLITS[local_j], phi_k, local_SPLITS[local_j]);
				local_ORDS[local_j] -= 1;
			}

			local_eval(W0, phi_k, *W0);
			local_eval(W1, phi_k, *W1);

			local_k += 1;
			local_i -= 1;
		}
	}

	while( local_k < level ){
		if(local_ORDS[local_i] != 1)
			spliting_scalars_Ceil(local_SPLITS, local_ORDS, &local_i, local_SM, E_k.A);
		else{
			local_get(&E_k, &phi_k, local_SPLITS[local_i], E_k);
			for(local_j = 0; local_j < local_i; local_j++){
				local_eval(&local_SPLITS[local_j], phi_k, local_SPLITS[local_j]);
				local_ORDS[local_j] -= 1;
			}

			local_eval(W0, phi_k, *W0);
			local_eval(W1, phi_k, *W1);

			local_k += 1;
			local_i -= 1;
		}
	}		
}
