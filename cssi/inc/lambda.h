#ifndef LAMBDA
#define LAMBDA

#include "IsogLib.h"

typedef struct node{
	struct node *next;
	uint64_t point[2];
} T_node;

// Function that create a new Node
void createNode(T_node *newNode, uint64_t new_point[2]);
// Function that gets the middle of a linked link
T_node *middleNode(T_node * startNode, T_node *endNode);
// Binary search in a linked list
T_node *binarySearch(T_node *head_ref, uint64_t current_point[2]);
// Function that inserts (using a binary search) in a sorted linked list
void sortedInsertion(T_node **head_ref, T_node *new_node);

typedef struct lambda{
	uint64_t Seed_point;	// X_0 : initial point
	uint64_t Last_point;	// X_l : distinguished  point
	uint64_t pathlength;	//   l : trail length
} T_lambda;

// ---- VARAIABLES TO BE USED

uint64_t DEGREE;
uint64_t e_DIVIDED_BY_2;
uint64_t VAR_M_MINUS_ONE;
uint64_t VAR_M;
uint64_t BITS_OF_e_DIVIDED_BY_2_PLUS_2;

Curve curve_E0;
Point point_P, point_Q,		// Generators of the torsion group E0(F_p^2)[d'^(e'/2)], and
		point_Pt, point_Qt,	// Generators of the torsion group E0(F_p^2)[d^(e/2)], and
		order_d_point_Qt;		// [d^(e/2 - 1)] point_P
		
Curve curve_E1;
Point point_phi_P, point_phi_Q,	//	Images of point_Pt and point_Qt
		point_St, point_Tt,			// Generators of the torsion group E1(F_p^2)[d^(e/2)]
		order_d_point_Tt;				// [d^(e/2 - 1)] point_T
		
void (*Pt_DSM)(), (*Pt_MDSM)();		//	Scalar multiplication by d and d^i.
void (*PohligHellman_alg)();			// Pohlig-Hellman algorithm
void (*get_isog)(), (*eval_isog)();	// Isogeny functions

uint64_t OMEGA_MINUS_ONE;
uint64_t BITS_OF_OMEGA;
uint64_t OMEGA;
uint64_t BETA;
double THETA;
double BITS_OF_R;

uint64_t NUMBER_OF_CORES;
uint64_t NONCE;
uint64_t MAXIMUM_NUMBER_OF_FUNCTIONS_PERMITED;
			
uint64_t BETA_TIMES_W;
uint64_t TRAIL_LENGTH;
uint64_t DISTINGUISHABILITY;

omp_lock_t writelock;	// Lock used to avoid overwriting and reading conflicts

// Hash table to be used in VW algorithm
T_lambda *HASH_TABLE_FOR_THE_DISTINGUISHED_POINTS_REACHED;	// This hash table will have OMEGA elements

// The following variables are required for counting the different collisions reached
T_node **DIFFERENT_COLLISIONS_REACHED;					// Each thread has a sorted linked list of the different collisions reached
T_node **ADDRESS_USED_FOR_THE_LIST_OF_COLLISIONS;
uint64_t *CURRENT_SIZE_OF_THE_LIST_OF_COLLISIONS;	// Current size of the list DIFFERENT_COLLISIONS_REACHED
uint64_t *CURRENT_NUMBER_OF_TOTAL_COLLISIONS;

// The following variables are required for the running time
uint64_t *FUNCTION_EVALUATIONS_FOR_GENERATING_DISTINGUISHED_POINTS;	// Each thread has a partial running time
uint64_t *FUNCTION_EVALUATIONS_FOR_LOCATING_COLLISIONS;					// Each thread has a partial running time
							  
// ----

// Reduction into a t-bits number
uint64_t MD5_tLSbits(uint8_t TYPE, uint64_t point, uint64_t mask, uint8_t t_bits);

// ---
// function h maps  an element of {0,1,2}x{0,1}^(e/2) into an order-2^(e/2) point
void FUNCTION_h(Point *v_kernel, 
					 Curve v_Curve, Point v_Pt, Point v_Qt, Point v_Qt_ord_D, 
					 uint64_t v_scalar[WORD_N], uint8_t seed);

void FUNCTION_gn(uint64_t output_point[WORD_N], uint64_t local_jinv[2][WORD_N]);
void FUNCTION_fn(uint64_t local_jinv[2][WORD_N], uint64_t new_scalar[WORD_N], uint64_t current_scalar[WORD_N]);
void get_collision(uint64_t COLLISION[2][WORD_N], T_lambda OLD, T_lambda NEW);
uint64_t VW_algorithm(uint64_t GOLDEN_COLLISION[2][WORD_N], uint8_t *BOOL, uint64_t kth_proc);
//---

#endif
