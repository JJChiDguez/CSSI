#include "CurveArith.h"

// ---------------------------------------------------------------
// jInvariant() function :::::::::::::::::::::::::::::::::::::::::
// inputs:
//        E : a supersingular elliptic curve defined over Fp2,
// output:
//        J : j(E)
void jInvariant(uint64_t J[2][WORD_N], Curve E){
	//Computes the jInvariant of the curve E 
	uint64_t U[2][WORD_N], D[2][WORD_N], r[2][WORD_N], s[2][WORD_N], t[2][WORD_N], u[2][WORD_N];
	uint64_t aux_temp[2][WORD_N];

	SQR_Fp2(U, E.A);// U := a^2;
	MUL_Fp2(aux_temp, U, E.A);// U := a^3;
	memcpy(U, aux_temp, sizeof(uint64_t) * WORD_N * 2);
	ADD_Fp2(t, U, U);// t := 2a^3;
	ADD_Fp2(U, t, t);// U := 4a^3

	SQR_Fp2(u, E.B);//u := b^2;
	ADD_Fp2(t, u, u);//t := 2b^2;
	ADD_Fp2(s, t, t);//s := 4b^2;
	ADD_Fp2(s, s, s);//s := 8b^2;
	ADD_Fp2(r, s, s);//r := 16b^2;
	ADD_Fp2(s, r, s);//s := 24b^2;
	ADD_Fp2(t, s, t);//t := 26b^2;
	ADD_Fp2(t, t, u);//t := 27b^2;
	ADD_Fp2(D, t, U);//D := 27b^2 + 4a^3

	INV_Fp2(D, D);//D := D^-1;

	ADD_Fp2(t, U, U);//t := 2U;
	ADD_Fp2(t, t, t);//t := 4U;
	ADD_Fp2(t, t, t);//t := 8U;
	ADD_Fp2(t, t, t);//t := 16U;
	ADD_Fp2(t, t, t);//t := 32U;
	ADD_Fp2(t, t, t);//t := 64U;
	ADD_Fp2(s, t, t);//s := 128U;
	ADD_Fp2(s, t, s);//s := 192U;
	ADD_Fp2(t, s, s);//t := 384U;
	ADD_Fp2(t, t, t);//t := 768U;
	ADD_Fp2(t, t, t);//t := 1536U;
	ADD_Fp2(U, t, s);//U := 1728U;

	MUL_Fp2(J, U, D);
}

// ---------------------------------------------------------------
// Point_Assign() function :::::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
// output:
//        Q : a copy of the point P
void Point_Assign(Point *Q, Point P){
	//assign point P to point Q 
	memcpy(Q->X, P.X, sizeof(uint64_t) * WORD_N * 2);
	memcpy(Q->Y, P.Y, sizeof(uint64_t) * WORD_N * 2);
	memcpy(Q->Z, P.Z, sizeof(uint64_t) * WORD_N * 2);
}

// ---------------------------------------------------------------
// is_equal_point() function :::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
//        Q : a point of a supersingular elliptic curve,
// output:
//        1 if P == Q, or
//        0 if P != Q
int is_equal_point(Point P, Point Q){
	//Return 1 if points are equals
	uint64_t X_P[2][WORD_N], X_Q[2][WORD_N], 
				Y_P[2][WORD_N], Y_Q[2][WORD_N],
				Z2_P[2][WORD_N], Z2_Q[2][WORD_N], 
				Z3_P[2][WORD_N], Z3_Q[2][WORD_N];

	SQR_Fp2(Z2_P, P.Z); MUL_Fp2(Z3_P, Z2_P, P.Z);	// P.Z^2 and P.Z^3
	SQR_Fp2(Z2_Q, Q.Z); MUL_Fp2(Z3_Q, Z2_Q, Q.Z);	// Q.Z^2 and Q.Z^3

	MUL_Fp2(X_P, Z2_Q, P.X); MUL_Fp2(Y_P, Z3_Q, P.Y);
	MUL_Fp2(X_Q, Z2_P, Q.X); MUL_Fp2(Y_Q, Z3_P, Q.Y);

	if( 	(compara(X_P[0], X_Q[0], WORD_N) == 0) && (compara(X_P[1], X_Q[1], WORD_N) == 0) &&
			(compara(Y_P[0], Y_Q[0], WORD_N) == 0) && (compara(Y_P[1], Y_Q[1], WORD_N) == 0) ){
		return 1;
	}
	else
		return 0;
}

// ---------------------------------------------------------------
// is_Zero_point() function ::::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
// output:
//        1 if P == Infinity point, or
//        0 if P != Infinity point
int is_Zero_point(Point P)
{
	//return 1 if P is equal to the infinity point
	if( (compara(P.Z[0], zeroM, WORD_N) == 0) && (compara(P.Z[1], zeroM, WORD_N) == 0) ){
		return 1;	
	}
	else
		return 0;
}

// ---------------------------------------------------------------
// Pt_DBL() function :::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
//        a : the coefficient a of a supersingular elliptic curve E
// output:
//       P2 : [2]P
void  Pt_DBL(Point *P2, Point P0, uint64_t a[2][WORD_N])
{    
	if(is_Zero_point(P0)){
		Point_Assign(P2, P0);
	}
	else{
	uint64_t XX[2][WORD_N], 
				YY[2][WORD_N], YYYY[2][WORD_N], 
				ZZ[2][WORD_N], 
				S[2][WORD_N], M[2][WORD_N], T[2][WORD_N], AUX[2][WORD_N];

		SQR_Fp2(XX, P0.X);		//   XX:=X1^2;
		SQR_Fp2(YY, P0.Y);		//   YY:=Y1^2;
		SQR_Fp2(YYYY, YY);		// YYYY:=YY^2;
		SQR_Fp2(ZZ, P0.Z);		//   ZZ:=Z1^2;
		ADD_Fp2(AUX, P0.X, YY);
		SQR_Fp2(S, AUX);
		SUB_Fp2(S, S, XX);
		SUB_Fp2(S, S, YYYY);
		ADD_Fp2(S, S, S);			//    S:=2*((X1+YY)^2-XX-YYYY);
		SQR_Fp2(AUX, ZZ);
		MUL_Fp2(M, a, AUX);
		ADD_Fp2(AUX, XX, XX);
		ADD_Fp2(AUX, AUX, XX);
		ADD_Fp2(M, AUX, M);		//    M:=3*XX+a*ZZ^2;
		SQR_Fp2(T, M);
		SUB_Fp2(T, T, S);
		SUB_Fp2(P2->X, T, S);	//   X3:=M^2-2*S;
		SUB_Fp2(AUX, S, P2->X);
		MUL_Fp2(P2->Y, M, AUX);
		ADD_Fp2(AUX, YYYY, YYYY);
		ADD_Fp2(AUX, AUX, AUX);
		ADD_Fp2(AUX, AUX, AUX);
		SUB_Fp2(P2->Y, P2->Y, AUX);	//   Y3:=M*(S-X3)-8*YYYY;
		ADD_Fp2(AUX, P0.Y, P0.Z);
		SQR_Fp2(P2->Z, AUX);
		SUB_Fp2(P2->Z, P2->Z, YY);
		SUB_Fp2(P2->Z, P2->Z, ZZ);	//   Z3:=(Y1+Z1)^2-YY-ZZ;
	}
}// 1M + 8S + 1D + 10add + 2times2 + 1times3 + 1times8: 

// ---------------------------------------------------------------
// Pt_ADD() function :::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
//        Q : a point of a supersingular elliptic curve,
//        a : the coefficient a of a supersingular elliptic curve E
// output:
//       R : P + Q
void  Pt_ADD(Point *R, Point P, Point Q, uint64_t a[2][WORD_N])
{
	if (is_Zero_point(P)){
		Point_Assign(R, Q);
	}
	else if(is_Zero_point(Q)){
		Point_Assign(R, P);
	}
	else if(is_equal_point(P,Q)){
		Pt_DBL(R, P, a);
	}
	else{
		// Given points P and Q Store in R the sum P+Q
		uint64_t Z2Z2Z2[2][WORD_N], Z2Z2[2][WORD_N], 
					Z1Z1[2][WORD_N], 
					AUX[2][WORD_N],
					U1[2][WORD_N], U2[2][WORD_N],
					S1[2][WORD_N], S2[2][WORD_N],
					H[2][WORD_N], I[2][WORD_N], J[2][WORD_N],
					r[2][WORD_N], V[2][WORD_N];


		SQR_Fp2(Z2Z2, Q.Z);		//   Z2Z2:=Z2^2;
		MUL_Fp2(Z2Z2Z2, Q.Z, Z2Z2);	// Z2Z2Z2:=Z2*Z2Z2;

		SQR_Fp2(Z1Z1, P.Z);		//   Z1Z1:=Z1^2;
		MUL_Fp2(U1, P.X, Z2Z2);	//     U1:=X1*Z2Z2;
		MUL_Fp2(U2, Q.X, Z1Z1);	//     U2:=X2*Z1Z1;
		MUL_Fp2(S1, P.Y, Z2Z2Z2);	//     S1:=Y1*Z2Z2Z2;
		MUL_Fp2(AUX, P.Z, Z1Z1);
		MUL_Fp2(S2, Q.Y, AUX);		//     S2:=Y2*Z1*Z1Z1;
		SUB_Fp2(H, U2, U1);		//      H:=U2-U1;
		ADD_Fp2(AUX, H, H);
		SQR_Fp2(I, AUX);		//      I:=(2*H)^2;
		MUL_Fp2(J, H, I);		//      J:=H*I;
		SUB_Fp2(r, S2, S1);
		ADD_Fp2(r, r, r);		//      r:=2*(S2-S1);
		MUL_Fp2(V, U1, I);		//      V:=U1*I;

		SQR_Fp2(R->X, r);
		SUB_Fp2(R->X, R->X, J);
		ADD_Fp2(AUX, V, V);
		SUB_Fp2(R->X, R->X, AUX);	//     X3:=r^2-J-2*V;

		SUB_Fp2(AUX, V, R->X);
		MUL_Fp2(R->Y, r, AUX);
		MUL_Fp2(AUX, S1, J);
		ADD_Fp2(AUX, AUX, AUX);
		SUB_Fp2(R->Y, R->Y, AUX);	//     Y3:=r*(V-X3)-2*S1*J;

		ADD_Fp2(AUX, P.Z, Q.Z);
		SQR_Fp2(R->Z, AUX);
		SUB_Fp2(AUX, R->Z, Z1Z1);
		SUB_Fp2(AUX, AUX, Z2Z2);
		MUL_Fp2(R->Z, AUX, H);	//     Z3:=((Z1+Z2)^2-Z1Z1-Z2Z2)*H;
	}
}// 10M + 4S + 9add + 4times2 after 1M + 1S for caching: 


// ---------------------------------------------------------------
// Pt_Neg() function :::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
// output:
//       Pm : -P
void Pt_Neg(Point *Pm, Point P){
	//Return -P into Pm
	if(is_Zero_point(P))
		Point_Assign(Pm, P);
	else{
		memcpy( Pm->X, P.X, sizeof(uint64_t) * WORD_N * 2);
		memcpy( Pm->Z, P.Z, sizeof(uint64_t) * WORD_N * 2);
		NEG_Fp2(Pm->Y, P.Y);
	}
}

// ---------------------------------------------------------------
// Pt_MDBL() function ::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
//        m : an integer number less than e, where [2^e]P = infinity point,
//        a : the coefficient a of a supersingular elliptic curve E
// output:
//      P2m : [2^m]P
void Pt_MDBL(Point * P2m, Point P,  uint64_t m, uint64_t a[2][WORD_N]){
	//Computes [2^m] P and store into P2m
	uint64_t i;
	Point Q; 
	Point_Assign(&Q, P);
	for(i=0 ; i < m; i++){
		Pt_DBL(&Q, Q, a);
	}
	Point_Assign(P2m, Q);
}

// ---------------------------------------------------------------
// DBLADD() function :::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
//        m : an integer number less than 2^e, where [2^e]P = infinity point,
//        a : the coefficient a of a supersingular elliptic curve E
// output:
//       mP : [m]P
void DBLADD(Point *mP, Point P, uint64_t m[WORD_N] , uint64_t a[2][WORD_N])
{
	Point Q, N;
	// Q <- (1,1,0) i.e. \infty
	memcpy( ((&Q)->X)[0],  unoM, sizeof(uint64_t) * WORD_N);
	memcpy( ((&Q)->X)[1], zeroM, sizeof(uint64_t) * WORD_N);
	memcpy( ((&Q)->Y)[0],  unoM, sizeof(uint64_t) * WORD_N);
	memcpy( ((&Q)->Y)[1], zeroM, sizeof(uint64_t) * WORD_N);
	memcpy( ((&Q)->Z)[0], zeroM, sizeof(uint64_t) * WORD_N);
	memcpy( ((&Q)->Z)[1], zeroM, sizeof(uint64_t) * WORD_N);
	// N <- P
	Point_Assign(&N, P);

	uint64_t v[WORD_N];
	memcpy(v, m, sizeof(uint64_t) * WORD_N);
	while( compara(v, zeroM, WORD_N) != 0){
		if (v[0]&0x0000000000000001){
			Pt_ADD(&Q, Q, N, a);
		}
		Pt_DBL(&N, N, a);
		SHR(v, v);
	};

	Point_Assign(mP, Q);
}

// ---------------------------------------------------------------
// Pt_TRPL() function ::::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
//        a : the coefficient a of a supersingular elliptic curve E
// output:
//       P2 : [3]P
void Pt_TRPL(Point *Q, Point P, uint64_t a[2][WORD_N])
{
	if (is_Zero_point(P)){
		Point_Assign(Q, P);
	}
	else{
		uint64_t XX[2][WORD_N],
					YY[2][WORD_N], YYYY[2][WORD_N],
					ZZ[2][WORD_N],
					M[2][WORD_N], MM[2][WORD_N],
					E[2][WORD_N], EE[2][WORD_N],
					T[2][WORD_N], U[2][WORD_N], 
					AUX1[2][WORD_N], AUX2[2][WORD_N];

		SQR_Fp2(XX, P.X);		//   XX:=X1^2;
		SQR_Fp2(YY, P.Y);		//   YY:=Y1^2;
		SQR_Fp2(ZZ, P.Z);		//   ZZ:=Z1^2;
		SQR_Fp2(YYYY, YY);		// YYYY:=YY^2;
		ADD_Fp2(AUX1, XX, XX);
		ADD_Fp2(AUX1, AUX1, XX);
		SQR_Fp2(AUX2, ZZ);
		MUL_Fp2(M, a, AUX2);
		ADD_Fp2(M, M, AUX1);		//    M:=3*XX+a*ZZ^2;
		SQR_Fp2(MM, M);		//   MM:=M^2;
		ADD_Fp2(AUX1, P.X, YY);
		SQR_Fp2(E, AUX1);
		SUB_Fp2(E, E, XX);
		SUB_Fp2(E, E, YYYY);
		ADD_Fp2(AUX1, E, E);
		ADD_Fp2(AUX1, AUX1, E);
		ADD_Fp2(E, AUX1, AUX1);
		SUB_Fp2(E, E, MM);		//    E:=6*((X1+YY)^2-XX-YYYY)-MM;
		SQR_Fp2(EE, E);		//   EE:=E^2;
		ADD_Fp2(T, YYYY, YYYY);
		ADD_Fp2(T, T, T);
		ADD_Fp2(T, T, T);
		ADD_Fp2(T, T, T);		//    T:=16*YYYY;
		ADD_Fp2(AUX1, M, E);
		SQR_Fp2(U, AUX1);
		SUB_Fp2(U, U, MM);
		SUB_Fp2(U, U, EE);
		SUB_Fp2(U, U, T);		//    U:=(M+E)^2-MM-EE-T;
		MUL_Fp2(Q->X, P.X, EE);
		MUL_Fp2(AUX1, YY, U);
		ADD_Fp2(AUX1, AUX1, AUX1);
		ADD_Fp2(AUX1, AUX1, AUX1);
		SUB_Fp2(Q->X, Q->X, AUX1);
		ADD_Fp2(Q->X, Q->X, Q->X);
		ADD_Fp2(Q->X, Q->X, Q->X);	//   X3:=4*(X1*EE-4*YY*U);
		SUB_Fp2(AUX1, T, U);
		MUL_Fp2(AUX2, U, AUX1);
		MUL_Fp2(AUX1, E, EE);
		SUB_Fp2(AUX2, AUX2, AUX1);
		MUL_Fp2(Q->Y, P.Y, AUX2);
		ADD_Fp2(Q->Y, Q->Y, Q->Y);
		ADD_Fp2(Q->Y, Q->Y, Q->Y);
		ADD_Fp2(Q->Y, Q->Y, Q->Y);	//   Y3:=8*Y1*(U*(T-U)-E*EE);
		ADD_Fp2(AUX1, P.Z, E);
		SQR_Fp2(Q->Z, AUX1);
		SUB_Fp2(Q->Z, Q->Z, ZZ);
		SUB_Fp2(Q->Z, Q->Z, EE);	//   Z3:=(Z1+E)^2-ZZ-EE;
	}
}	// 5M + 10S + 1D + 15add + 1times3 + 2times4 + 1times6 + 1times8 + 1times16:

// ---------------------------------------------------------------
// Pt_MTRPL() function :::::::::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
//        m : an integer number less than e, where [3^e]P = infinity point,
//        a : the coefficient a of a supersingular elliptic curve E
// output:
//      P2m : [3^m]P
void Pt_MTRPL(Point * P3m, Point P,  uint64_t m, uint64_t a[2][WORD_N]){
	//Computes [3^m] P and store into P3m
	uint64_t i;
	Point Q; 
	Point_Assign(&Q, P);
	for(i=0 ; i < m; i++){
		Pt_TRPL(&Q, Q, a);
	}
		Point_Assign(P3m, Q);
}

// ---------------------------------------------------------------
// DLP_order_two() function ::::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
//        Q : a point of a supersingular elliptic curve,
//     EXPN : the minimal integer such that [2^EXPN]P = infinity point,
//        a : the coefficient a of a supersingular elliptic curve E
// output:
//      DLP : an integer such that Q = [DLP]P
void   DLP_order_two(uint64_t DLP[WORD_N], Point local_Q, Point local_P, uint64_t EXPN, uint64_t a[2][WORD_N]){
	int64_t acc_two = 1;
	memcpy(DLP, zeroM, sizeof(uint64_t) * WORD_N);

	Point P_aux, Q_acc, Q_aux;
	int local_i;

	Pt_MDBL(&P_aux, local_P, EXPN - 1, a);

	for(local_i = 0; local_i < EXPN; local_i++){
		DBLADD(&Q_aux, local_P,    DLP, a);
		Pt_Neg(&Q_acc, Q_aux);
		Pt_ADD(&Q_aux, local_Q, Q_acc, a);
		Pt_MDBL(&Q_acc, Q_aux, EXPN - 1 - local_i, a);

		// We ask if the point Q_acc is not the infinity point
		if( is_equal_point(Q_acc, P_aux))
			DLP[0] += acc_two;

		acc_two *= 2;
	}
}

// ---------------------------------------------------------------
// DLP_order_three() function ::::::::::::::::::::::::::::::::::::
// inputs:
//        P : a point of a supersingular elliptic curve,
//        Q : a point of a supersingular elliptic curve,
//     EXPN : the minimal integer such that [3^EXPN]P = infinity point,
//        a : the coefficient a of a supersingular elliptic curve E
// output:
//      DLP : an integer such that Q = [DLP]P
void DLP_order_three(uint64_t DLP[WORD_N], Point local_Q, Point local_P, uint64_t EXPN, uint64_t a[2][WORD_N])
{
	int64_t acc_three = 1;
	memcpy(DLP, zeroM, sizeof(uint64_t) * WORD_N);

	Point P_aux, P2_aux, Q_acc, Q_aux;
	int local_i;

	Pt_MTRPL(&P_aux, local_P, EXPN - 1, a);
	Pt_DBL(&P2_aux, P_aux, a);

	for(local_i = 0; local_i < EXPN; local_i++){
		DBLADD(&Q_aux, local_P,    DLP, a);
		Pt_Neg(&Q_acc, Q_aux);
		Pt_ADD(&Q_aux, local_Q, Q_acc, a);
		Pt_MTRPL(&Q_acc, Q_aux, EXPN - 1 - local_i, a);

		// We ask if the point Q_acc is not the infinity point
		if( is_equal_point(Q_acc, P_aux))
			DLP[0] += acc_three;
		else if( is_equal_point(Q_acc, P2_aux))
			DLP[0] += (2 * acc_three);

		acc_three *= 3;
	}
}
