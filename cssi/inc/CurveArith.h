#ifndef ELLIPTICCURVE
#define ELLIPTICCURVE

#include "FieldArith.h"

typedef struct ECurve{
    uint64_t A[2][WORD_N];
    uint64_t B[2][WORD_N];
} Curve;

typedef struct Afin{
    uint64_t X[2][WORD_N];
    uint64_t Y[2][WORD_N];
    uint64_t Z[2][WORD_N];
} Point;

// Public curve parameters
#include "setup_EC.h"

// ------------------------------------------------
void jInvariant(uint64_t J[2][WORD_N], Curve E);
void Point_Assign(Point *Q, Point P);
int is_equal_point(Point P, Point Q);
int is_Zero_point(Point P);

// ------------------------------------------------
void  Pt_Neg(Point *Pm, Point P);
void  Pt_DBL(Point *P2, Point P0, uint64_t a[2][WORD_N]);
void  Pt_ADD(Point *R, Point P, Point Q, uint64_t a[2][WORD_N]);
void Pt_TRPL(Point *Q, Point P, uint64_t a[2][WORD_N]);

void  Pt_MDBL(Point *P2m, Point P,  uint64_t m, uint64_t a[2][WORD_N]);
void   DBLADD(Point *mP, Point P, uint64_t m[WORD_N] , uint64_t a[2][WORD_N]);
void Pt_MTRPL(Point * P3m, Point P,  uint64_t m, uint64_t a[2][WORD_N]);

void   DLP_order_two(uint64_t DLP[WORD_N], Point local_Q, Point local_P, uint64_t EXPN, uint64_t a[2][WORD_N]);
void DLP_order_three(uint64_t DLP[WORD_N], Point local_Q, Point local_P, uint64_t EXPN, uint64_t a[2][WORD_N]);

#endif
