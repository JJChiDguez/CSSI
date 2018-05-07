#ifndef ISOGENY
#define ISOGENY

#include "CurveArith.h"

typedef struct Isom{
    uint64_t Z[2][WORD_N];
} Isomorphism;

typedef struct Isog{
    uint64_t Gxza[2][WORD_N];
    uint64_t Kx[2][WORD_N];
    uint64_t Ky[2][WORD_N];
    uint64_t KyKy4[2][WORD_N];
    uint64_t Kz[2][WORD_N];
    uint64_t KzKz[2][WORD_N];
    uint64_t KzKzKz[2][WORD_N];
} Isogeny;

// ------------------------------------------------
void get_isomorphism(Isomorphism *I, Curve E0, Curve E1);
void eval_isomorphism(Point *Q, Isomorphism I, Point P);

void get_2_isog(Curve *EA, Isogeny *f, Point K, Curve E0);
void eval_2_isog(Point *PA, Isogeny f, Point P0);

void get_3_isog(Curve *EA, Isogeny *f, Point K, Curve E0);
void eval_3_isog(Point *Q, Isogeny f, Point P);

void spliting_scalars_Floor(Point SPLITS[Log2_E], uint64_t ORDS[Log2_E], uint64_t *lenght, void (*local_SM)(), uint64_t a[2][WORD_N]);
void  spliting_scalars_Ceil(Point SPLITS[Log2_E], uint64_t ORDS[Log2_E], uint64_t *lenght, void (*local_SM)(), uint64_t a[2][WORD_N]);

void get_isogenous_curve(Curve *E_k, Point R, void (*local_SM)(), void (*local_get)(), void (*local_eval)(), Curve E0, uint64_t level);
void eval_isogeny(Point *W0, Point *W1, Point R, void (*local_SM)(), void (*local_get)(), void (*local_eval)(), Curve E0, uint64_t level);

#endif
