#ifndef NAIVE
#define NAIVE

#include "IsogLib.h"

//	Hash table element: an structure which fields are
// the j-invariant and the scalar (kernel)
typedef struct naive{
	 uint64_t  j_inv[2][WORD_N];
	 uint64_t             point;
} T_naive;

// ---- VARAIABLES TO BE USED

uint64_t DEGREE;
uint64_t e_DIVIDED_BY_2;
uint64_t VAR_N_DIVIDED_BY_DEGREE_PLUS_1;
uint64_t VAR_N;

Curve curve_E0;
Point point_P, point_Q,		// Generators of the torsion group E0(F_p^2)[d'^(e'/2)], and
		point_Pt, point_Qt,	// Generators of the torsion group E0(F_p^2)[d^(e/2)], and
		d_times_Pt,				
		order_d_point_Qt;		// [d^(e/2 - 1)] point_Qt
		
Curve curve_E1;
Point point_phi_P, point_phi_Q,	//	Images of point_Pt and point_Qt
		point_St, point_Tt,			// Generators of the torsion group E1(F_p^2)[d^(e/2)]
		d_times_St,
		order_d_point_Tt;				// [d^(e/2 - 1)] point_Tt
		
void (*Pt_DSM)(), (*Pt_MDSM)();		//	Scalar multiplication by d and d^i.
void (*PohligHellman_alg)();			// Pohlig-Hellman algorithm
void (*get_isog)(), (*eval_isog)();	// Isogeny functions

uint64_t NUMBER_OF_CORES;
uint64_t *RUNNING_TIME;
uint64_t *CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0;
T_naive **TREE_LEAVES_ROOTED_AT_E0;							// Each thread has a sorted linked list 

// -----

// Sorting and searching functions
void swap(T_naive *a, T_naive *b);
int partition (T_naive arr[], int low, int high);
void quickSort(T_naive arr[], int low, int high);
int64_t binarySearch(T_naive local_seq[], uint64_t local_jinv[2][WORD_N], uint64_t low, uint64_t high);

uint64_t Leaves_of_the_tree_rooted_at_E0(uint64_t kth_proc);
uint64_t Leaves_of_the_tree_rooted_at_E1(uint64_t COLLISION[2][WORD_N], uint8_t *BOOL, uint64_t kth_proc);

#endif
