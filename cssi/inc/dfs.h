#ifndef DFS
#define DFS

#include "IsogLib.h"

// Tree structure
typedef struct dfs{
	Curve Ei;
	Point Pi;
	Point Qi;
} T_dfs;

// Leaves of the tree
typedef struct naive{
	 uint64_t  j_inv[2][WORD_N];
	 uint64_t             point;
} T_naive;

// ---- VARAIABLES TO BE USED

uint64_t DEGREE;
uint64_t e_DIVIDED_BY_2;
uint64_t VAR_N_DIVIDED_BY_DEGREE_PLUS_1;
uint64_t VAR_N;
uint64_t LEVEL;
uint64_t e_DIVIDED_BY_2_MINUS_LEVEL;

Curve curve_E0;
Point point_P, point_Q,		// Generators of the torsion group E0(F_p^2)[d'^(e'/2)], and
		point_Pt, point_Qt;	// Generators of the torsion group E0(F_p^2)[d^(e/2)]
		
Curve curve_E1;
Point point_phi_P, point_phi_Q,	//	Images of point_Pt and point_Qt
		point_St, point_Tt;			// Generators of the torsion group E1(F_p^2)[d^(e/2)]
		
void (*Pt_DSM)(), (*Pt_MDSM)();		//	Scalar multiplication by d and d^i.
void (*PohligHellman_alg)();			// Pohlig-Hellman algorithm
void (*get_isog)(), (*eval_isog)();	// Isogeny functions

uint64_t NUMBER_OF_CORES;
uint64_t *CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0;
T_naive **TREE_LEAVES_ROOTED_AT_E0;							// Each thread has a sorted linked list 

Point KERNELS_AT_E0[3], KERNEL_AT_E0[3];
T_dfs INITIAL_NODES_AT_E0[3];			// This will determine the three 2-isogenous curves to E0
Isogeny PHI_AT_E0[3], PHI_LOCAL;

Point KERNELS_AT_E1[3], KERNEL_AT_E1[3];
T_dfs INITIAL_NODES_AT_E1[3];			// This will determine the three 2-isogenous curves to E0
Isogeny PSI_AT_E1[3], PSI_LOCAL;
Curve ISOMORPHIC_CURVE_TO_E1;
uint64_t *CURRENT_DOUBLINGS;
uint64_t *CURRENT_ADDITIONS;
uint64_t *CURRENT_ISOG_EVAL;
uint64_t *CURRENT_ISOG_COMP;
		
// Sorting and searching functions
void swap(T_naive *a, T_naive *b);
int partition (T_naive arr[], int low, int high);
void quickSort(T_naive arr[], int low, int high);
int64_t binarySearch(T_naive local_seq[], uint64_t local_jinv[2][WORD_N], uint64_t low, uint64_t high);

// Evaluating 2^(e/2)-isogenies
Curve isogenies_decomposition_at_E0(uint64_t X, Point *output_0, Point *output_1, Point input_0, Point input_1);
Curve isogenies_decomposition_at_E1(uint64_t X, Point *output_0, Point *output_1, Point input_0, Point input_1);
// Initial node for each processor
void build_initial_nodes(T_dfs NODES[], uint64_t Xs[], T_dfs NODE, uint64_t X, uint64_t level, uint64_t kth_proc);
// Depth First Search applied for computing 2^(e/2)-isogenies
void Depth_First_Search_at_E0(T_dfs NODE, uint64_t X, uint64_t level, uint64_t kth_proc);
void Depth_First_Search_at_E1(uint8_t *FLAG, uint64_t COLLISION[2], T_dfs NODE, uint64_t X, uint64_t level, uint64_t kth_proc);

#endif
