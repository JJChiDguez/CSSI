#include "lambda.h"

// ---------------------------------------------------------------
// ---------------------------------------------------------------
// The following functions are required for storing the collisions
// (without repetitions)

// ---------------------------------------------------------------
// createNode() function :::::::::::::::::::::::::::::::::::::::::
// inputs:
//         new_point : a list of 64-t bits integers, which has size
//                     equals two
// output:
//         newNode : newNode->point = new_point
void createNode(T_node *newNode, uint64_t new_point[2])
{
	newNode->point[0] = new_point[0];
	newNode->point[1] = new_point[1];
	newNode->next = NULL;
};

// ---------------------------------------------------------------
// middleNode() function :::::::::::::::::::::::::::::::::::::::::
// inputs:
//         startnode : an element of a sorted linked list of elements 
//                     of type T_node
//         endNode : an element of a list of elements of type T_node
// output:
//         the middle node between startNode and endNode
T_node *middleNode(T_node * startNode, T_node *endNode)
{ 
	if( startNode == NULL ) 
	{
		//  If the linked list is empty 
		return NULL; 
	}
	T_node *slowPtr = startNode; 
	T_node *fastPtr = startNode->next; 
	while(fastPtr != endNode)
	{ 
		fastPtr = fastPtr->next;
		if(fastPtr != endNode) 
		{ 
			slowPtr = slowPtr->next ;
			fastPtr = fastPtr->next ;
			// Note that for each loop iteration, slowPtr moves just 
			// one location  while fastPtr moves two nodes at a time.  
		}
	}
	return slowPtr;
	//	At the end, the slowPtr will be  pointing to the middle node 
} 

// ---------------------------------------------------------------
// binarySearch() function :::::::::::::::::::::::::::::::::::::::
// inputs:
//         head : a sorted linked list of elements of type T_node
//         current_point : a list of 64-t bits integers, which has 
//                         equals two size
// output:
//         an element NODE of type T_node If NODE->point = current_point 
//                                           for some element of the
//                                           sorted linked list head
//         NULL otherwise.                
T_node *binarySearch(T_node *head, uint64_t current_point[2])
{ 
	T_node *current;
	T_node *startNode = head;
	T_node *endNode = NULL; 
	uint64_t current_id;
	if(head == NULL)
		return NULL;
	else
	{
		/* Locate the node before the point of insertion */
		current = middleNode(startNode, endNode);
		while( (current != NULL) && (startNode != endNode) )
		{
			if(current == NULL)
				return NULL;
			if(compara(current->point, current_point, 2) == 0)
				return current;
			else if(compara(current->point, current_point, 2) == -1)
				startNode = current->next;
			else
				endNode = current;
			
			if(startNode != endNode)
				current = middleNode(startNode, endNode);
		}
		// data not present 
		return NULL; 
	}
};

// ---------------------------------------------------------------
// middleNode() function :::::::::::::::::::::::::::::::::::::::::
// inputs:
//         head_ref : a sorted linked list of elements of type T_node
//         new_node : an element of type T_node
// output:
//         head_ref which has the element new_node (it was inserted 
//         in a sorted way)
void sortedInsertion(T_node **head_ref, T_node *new_node)
{
	/* Special case for the head end */
	if(*head_ref == NULL)
	{
		new_node->next = *head_ref;
		*head_ref = new_node;
	}
	else if(compara((*head_ref)->point, new_node->point, 2) >= 1)
	{
		new_node->next = *head_ref;
		*head_ref = new_node;
	}
	else
	{
		/* Locate the node before the point of insertion */
		T_node *current;
		T_node *startNode = *head_ref;
		T_node *endNode = NULL;
		current = middleNode(startNode, endNode);

		while( (current->next != NULL) && (startNode != endNode) )
		{
			if(compara(current->next->point, new_node->point, 2) == -1)
				startNode = current->next;
			else
				endNode = current;
				
			if(startNode != endNode)
				current = middleNode(startNode, endNode);
		}
		
		new_node->next = startNode->next;
		startNode->next = new_node;
	}
};

// ---------------------------------------------------------------
// ---------------------------------------------------------------

// ---------------------------------------------------------------
// MD5_tLSbits() function ::::::::::::::::::::::::::::::::::::::::
// inputs:
//          TYPE : an element of {0,1}^8
//         point : an element of {0,1}x{0,1,2}x{0,1}^(e/2 - 1)
//             t : number of bits to be required
//          mask : 2^t - 1
// output:
//         t-LEAST SIGNIFICANT BITS OF MD5(TYPE||point) output
uint64_t MD5_tLSbits(uint8_t TYPE, uint64_t point, uint64_t mask, uint8_t t_bits)
{
	// We transform the point into an string
	unsigned char string_input[ 2 * 16 + 1], v_string[16 + 1];
	string_input[0] = '\0';
	
	v_string[0] = '\0';
	sprintf(v_string, "%x", TYPE);
	v_string[strlen(v_string)] = '\0';
	strcat(string_input, v_string);
	string_input[strlen(string_input)] = '\0';
	
	v_string[0] = '\0';
	sprintf(v_string, "%jx", point);
	v_string[strlen(v_string)] = '\0';
	strcat(string_input, v_string);
	string_input[strlen(string_input)] = '\0';
			
	unsigned char v_md5[MD5_DIGEST_LENGTH];
	
	// Now, we apply MD5
	MD5(string_input, strlen(string_input), v_md5);
	
	// We transform the 32 least significant bits of MD5(point) into an 64-bits integer
	int i = 0, j = 0;
	uint64_t point_tLSbits_MD5 = 0x0000000000000000;
	while(j < t_bits)
	{
		point_tLSbits_MD5 ^= ( (uint64_t)(0xFF & v_md5[MD5_DIGEST_LENGTH - 1 - i]) ) << j;
		j += 8; i += 1;
	}
	return point_tLSbits_MD5 & mask;
};

// ---------------------------------------------------------------
// FUNCTION_h() function :::::::::::::::::::::::::::::::::::::::::
// inputs:
//               v_Curve : a supersingular elliptic curve
//         v_Pt and v_Qt : a generators of v_Curve(F_p^2)[2^(e/2)]
//            v_Qt_ord_D : [2^(e/2 - 1)]v_Qt
//              v_scalar : an element of {0,1}^(e/2 - 1)
//                  seed : an element of {0,1,2}
// output:
//         v_kernel : a order-2^(e/2) point in v_Curve(F_p^2)
void FUNCTION_h(Point *v_kernel, 
					 Curve v_Curve, Point v_Pt, Point v_Qt, Point v_Qt_ord_D, 
					 uint64_t v_scalar[WORD_N], uint8_t seed)
{
	if(seed == 0x00)
	{
		// v_Pt + [n]v_Qt
		DBLADD(v_kernel, v_Qt, v_scalar, v_Curve.A);
		Pt_ADD(v_kernel, *v_kernel, v_Pt, v_Curve.A);
	}
	else if(seed == 0x01)
	{
		// v_Pt + [d^(e/2 - 1) + n]v_Qt
		// Here, v_Qt2 must be equals [2^(e/2 - 1)]v_Qt;
		DBLADD(v_kernel, v_Qt, v_scalar, v_Curve.A);
		Pt_ADD(v_kernel, *v_kernel, v_Pt, v_Curve.A);
		Pt_ADD(v_kernel, *v_kernel, v_Qt_ord_D, v_Curve.A);
	}
	else if(seed == 0x02)
	{
		// [2*n]v_Pt + v_Qt
		DBLADD(v_kernel, v_Pt, v_scalar, v_Curve.A);
		Pt_DSM(v_kernel, *v_kernel, v_Curve.A);
		Pt_ADD(v_kernel, *v_kernel, v_Qt, v_Curve.A);
	}
};

// ---------------------------------------------------------------
// FUNCTION_gn() function ::::::::::::::::::::::::::::::::::::::::
// inputs:
//         local_jinv : an element of F_{p^2}
// output:
//         output_point : an element of {0,1}x{0,1,2}x{0,1}^(e/2 - 1)

void FUNCTION_gn(uint64_t output_point[WORD_N], uint64_t local_jinv[2][WORD_N])
{
	uint64_t local_i, local_j;
	
	uint64_t TYPE = 1;
	unsigned char string_input[(WORD_N * 2 + 3)* 16 + 1], v_string[16 + 1];
	string_input[0] = '\0';
	
	v_string[0] = '\0';
	sprintf(v_string, "%jX", TYPE);
	v_string[strlen(v_string)] = '\0';
	strcat(string_input, v_string);
	string_input[strlen(string_input)] = '\0';
		
	// We need to transform the j-invariant into a string
	
	// Real part of the j-invariant
	unsigned char REAL_PART[WORD_N * 2 * 16 + 1];
	print_hex(REAL_PART, local_jinv[0], WORD_N);
	strcat(string_input, REAL_PART);
	string_input[strlen(string_input)] = '\0';
	
	// Imaginary part of the j-invariant
	unsigned char IMAGINARY_PART[WORD_N * 2 * 16 + 1];
	print_hex(IMAGINARY_PART, local_jinv[1], WORD_N);
	strcat(string_input, IMAGINARY_PART);
	string_input[strlen(string_input)] = '\0';	
	
	// We append the NONCE to string_input
	// each NONCE determines a different function to be used
	v_string[0] = '\0';
	sprintf(v_string, "%jX", NONCE);
	v_string[strlen(v_string)] = '\0';
	strcat(string_input, v_string);
	string_input[strlen(string_input)] = '\0';
	
	uint64_t local_N = strlen(string_input);
	
	// Reduction algorithm: g
	uint64_t local_counter = 0;
	unsigned char v_md5[MD5_DIGEST_LENGTH];
	
	// In this part we append the counter to (1||j-invariant || NONCE)
	string_input[local_N] = '\0';
	v_string[0] = '\0';
	sprintf(v_string, "%jX", local_counter);
	v_string[strlen(v_string)] = '\0';
	strcat(string_input, v_string);
	string_input[strlen(string_input)] = '\0';
	
	// Now, we apply MD5
	MD5(string_input, strlen(string_input), v_md5);
	local_counter += 1;
		
	while( (uint8_t)(0x03 & v_md5[MD5_DIGEST_LENGTH - 1]) == 0x03)
	{
		// In this part we append the counter to (j-invariant || NONCE)
		string_input[local_N] = '\0';
		v_string[0] = '\0';
		sprintf(v_string, "%jX", local_counter);
		v_string[strlen(v_string)] = '\0';
		strcat(string_input, v_string);
		string_input[strlen(string_input)] = '\0';
		
		// Now, we compute MD5(j-invariant||NONCE||counter)
		MD5(string_input, strlen(string_input), v_md5);
		local_counter += 1;
	}
	
	// We save the t least significant bit of MD5(j-invariant||NONCE||counter)
	local_i = 0, local_j = 0;
	uint64_t point_tLSbits_MD5 = 0x0000000000000000;
	while(local_j < BITS_OF_e_DIVIDED_BY_2_PLUS_2)
	{
		point_tLSbits_MD5 ^= ( (uint64_t)(0xFF & v_md5[MD5_DIGEST_LENGTH - 1 - local_i]) ) << local_j;
		local_j += 8; local_i += 1;
	}
	
	output_point[0] = point_tLSbits_MD5 & VAR_M_MINUS_ONE;
};


// ---------------------------------------------------------------
// FUNCTION_fn() function ::::::::::::::::::::::::::::::::::::::::
// inputs:
//         current_scalar : an element of {0,1}x{0,1,2}x{0,1}^(e/2 - 1)
// output:
//         local_jinv : the j-invariant associated with current_scalar
//         new_scalar : an element of {0,1}x{0,1,2}x{0,1}^(e/2 - 1)

void FUNCTION_fn(uint64_t local_jinv[2][WORD_N], uint64_t new_scalar[WORD_N], uint64_t current_scalar[WORD_N])
{
	Curve local_curve;
	Point pt_kernel;
	uint8_t local_seed = (uint8_t)(current_scalar[0] & 0x0000000000000007);
	uint64_t local_scalar[WORD_N];
	memcpy(local_scalar, current_scalar, sizeof(uint64_t) * WORD_N);
	// Which side: E or E_isog
	if( (local_seed & 0x04) != 0x00 )
	{
		memcpy(local_curve.A, curve_E1.A, sizeof(uint64_t) * WORD_N * 2);
		memcpy(local_curve.B, curve_E1.B, sizeof(uint64_t) * WORD_N * 2);
		// ---
		local_scalar[0] = local_scalar[0] >> 3;
		FUNCTION_h(&pt_kernel, local_curve, point_St, point_Tt, order_d_point_Tt, local_scalar, local_seed & 0x03);
	}
	else{
		memcpy(local_curve.A, curve_E0.A, sizeof(uint64_t) * WORD_N * 2);
		memcpy(local_curve.B, curve_E0.B, sizeof(uint64_t) * WORD_N * 2);
		// ---
		local_scalar[0] = local_scalar[0] >> 3;
		FUNCTION_h(&pt_kernel, local_curve, point_Pt, point_Qt, order_d_point_Qt, local_scalar, local_seed & 0x03);
	}
		
	// In this part, we compute the isogeny and its j-invariant
	Curve new_curve;
	get_isogenous_curve(&new_curve, pt_kernel, Pt_MDSM, get_isog, eval_isog, local_curve, e_DIVIDED_BY_2);
	jInvariant(local_jinv, new_curve);
			
	// Now, we get the new state F(j-invariant)
	FUNCTION_gn(new_scalar, local_jinv);
};


// ---------------------------------------------------------------
// get_collision() function ::::::::::::::::::::::::::::::::::::::
// inputs:
//         OLD : a triple (OLD_{DISTINGUISHED POINT}, OLD_{INITIAL STATE}, OLD_{TRAIL_LENGTH})
//         NEW : a triple (NEW_{DISTINGUISHED POINT}, NEW_{INITIAL STATE}, NEW_{TRAIL_LENGTH})
// Here, the follosing four elements belong to {0,1}x{0,1,2}x{0,1}^(e/2 - 1)
//      OLD_{DISTINGUISHED POINT}, OLD_{INITIAL STATE},
//      NEW_{DISTINGUISHED POINT}, and NEW_{INITIAL STATE}.
// In addition, 
//      OLD_{DISTINGUISHED POINT} = FUNCTION_fn^(OLD_{TRAIL_LENGTH}) (OLD_{INITIAL STATE}),
//      NEW_{DISTINGUISHED POINT} = FUNCTION_fn^(NEW_{TRAIL_LENGTH}) (NEW_{INITIAL STATE}),
// output:
//         COLLISION : a collision, i.e., an element of 
//                     [{0,1}x{0,1,2}x{0,1}^(e/2 - 1)] x [{0,1}x{0,1,2}x{0,1}^(e/2 - 1)]

void get_collision(uint64_t COLLISION[2][WORD_N], T_lambda OLD, T_lambda NEW)
{
	uint64_t scalar_0[WORD_N], jinv_0[2][WORD_N];
	uint64_t scalar_1[WORD_N], jinv_1[2][WORD_N];
	int64_t local_i;
	
	memcpy(jinv_0[0], zeroM, sizeof(uint64_t) * WORD_N);
	memcpy(jinv_0[1], unoM, sizeof(uint64_t) * WORD_N);
	memcpy(jinv_1[0], unoM, sizeof(uint64_t) * WORD_N);
	memcpy(jinv_1[1], zeroM, sizeof(uint64_t) * WORD_N);
		
	int64_t local_len, local_dif;
	
	// We need to see which trajectory is the shortest
	memcpy(scalar_0, zeroM, sizeof(uint64_t) * WORD_N);
	memcpy(scalar_1, zeroM, sizeof(uint64_t) * WORD_N);
		
	if(OLD.pathlength <= NEW.pathlength)
	{
		local_len = OLD.pathlength;
		local_dif = NEW.pathlength - OLD.pathlength;
		
		scalar_0[0] = NEW.Seed_point;
		scalar_1[0] = OLD.Seed_point;
	}
	else
	{
		local_len = NEW.pathlength;
		local_dif = OLD.pathlength - NEW.pathlength;

		scalar_1[0] = NEW.Seed_point;
		scalar_0[0] = OLD.Seed_point;
	}

	// Now, we need to move into the largest trajectory until
	// obtain the same convergence distance in both trails.
	
	uint64_t scalar_00[WORD_N], scalar_11[WORD_N];
	memcpy(scalar_00, zeroM, sizeof(uint64_t) * WORD_N);
	memcpy(scalar_11, zeroM, sizeof(uint64_t) * WORD_N);
	uint8_t seed_00, seed_11;
		
	// We move through the largest trajectory
	for(local_i = 0; local_i < local_dif; local_i++)
	{
		scalar_00[0] = scalar_0[0];
		FUNCTION_fn(jinv_0, scalar_0, scalar_00);
	}
	
	// If the smallest trail is contained in the longest trail then, we return 0.
	if( compara(scalar_0, scalar_1, WORD_N) == 0 )
	{
		COLLISION[0][0] = 0x0000000000000000;
		COLLISION[1][0] = 0x0000000000000000;
	}
	else
	{
		scalar_00[0] = scalar_0[0];
		FUNCTION_fn(jinv_0, scalar_0, scalar_00);
		
		scalar_11[0] = scalar_1[0];
		FUNCTION_fn(jinv_1, scalar_1, scalar_11);
			
		// We search the first collision between the trails	
		local_i = 0;
		while(local_i < local_len)
		{
			if(compare_in_Fp2(jinv_0, jinv_1) == 0)
				break;
				
			scalar_00[0] = scalar_0[0];
			scalar_11[0] = scalar_1[0];
				
			FUNCTION_fn(jinv_0, scalar_0, scalar_00);
			FUNCTION_fn(jinv_1, scalar_1, scalar_11);
			
			local_i += 1;
		};

		// We saved the collision
		COLLISION[0][0] = scalar_00[0];
		COLLISION[1][0] = scalar_11[0];
	}
};


// ---------------------------------------------------------------
// VW_algorithm() function ::::::::::::::::::::::::::::::::::::::
// inputs:
//         BOOL : this variable determines if the golden collision
//                has been found 
//         kth_proc : the core ID
// output:
//         GOLDEN_COLLISION : a collision, i.e., an element of 
//                     [{0,1}x{0,1,2}x{0,1}^(e/2 - 1)] x [{0,1}x{0,1,2}x{0,1}^(e/2 - 1)]
//         The RUNNING TIME (i.e., function evaluations)

uint64_t VW_algorithm(uint64_t GOLDEN_COLLISION[2][WORD_N], uint8_t *BOOL, uint64_t kth_proc)
{
	clock_t local_start, local_end;
	
	uint64_t current_n[WORD_N], saved_n[WORD_N];
	memcpy(current_n, zeroM, sizeof(uint64_t) * WORD_N);
	memcpy(  saved_n, zeroM, sizeof(uint64_t) * WORD_N);
	
	uint64_t current_jinv[2][WORD_N];
	
	T_lambda current_node, saved_node;
	uint64_t local_COLLISION[2][WORD_N];
	
	// The parameters to be used in the initial (random) state.
	uint64_t local_j, acc_j, position, ctr_j;
	double local_rnd;
  
	local_start = get_cycles();
	
	// The following three variable are used to make a copy of the golden collision
	uint64_t local_point = 0,
				local_colls = 0,
				diffe_colls = 0,
				FunEvals_col = 0,
				FunEvals_pts = 0;
	
	
	uint64_t new_co[2], local_tmp, local_idx;
	T_node *node_co;

	while( (*BOOL == 0x01) && (local_point < BETA_TIMES_W) )
	{
		// Initial random element in Fp2, this is used to generate the first state in {0,1}x{0,1,2}x{0,1}^(e/2 - 1)
		for(local_j = 0; local_j < WORD_N; local_j++)
			current_jinv[0][local_j] = ((uint64_t)rand() & 0x00000000FFFFFFFF) ^ ( ((uint64_t)rand() & 0x00000000FFFFFFFF) << 32);
			
		for(local_j = 0; local_j < WORD_N; local_j++)
			current_jinv[1][local_j] = ((uint64_t)rand() & 0x00000000FFFFFFFF) ^ ( ((uint64_t)rand() & 0x00000000FFFFFFFF) << 32);
		
		// The initial state is determined by the scalar saved_n and the seed saved_seed
		FUNCTION_gn(saved_n, current_jinv);
		
		// We copy the initial state. Such copy will be modified until to find a distinguished point
		current_n[0] = saved_n[0];
		
		// Now, we evaluate F (at most 10 * TRAIL_LENGHT times)
		ctr_j = 0;
		while( (int)ctr_j < (int)(10 * TRAIL_LENGTH) )
		{
			FUNCTION_fn(current_jinv, current_n, current_n);
			ctr_j += 1;
			position = ( (0x07 & current_n[0]) << (BITS_OF_e_DIVIDED_BY_2_PLUS_2 - 3) ) ^ (current_n[0] >> 3);
			// Is a distinguished point?
			if(MD5_tLSbits(2, position, 0x00000000FFFFFFFF, 32) <= DISTINGUISHABILITY)
				break;
		}
		FunEvals_pts += (uint64_t)ctr_j;	// #(Function evaluations for finding a distinguished point)
		
		// Is a distinguished point?
		if(MD5_tLSbits(2, position, 0x00000000FFFFFFFF, 32) <= DISTINGUISHABILITY)
		{
			local_point += 1;	// #(Distinguished points)
			
			// The new triple
			current_node.Seed_point = saved_n[0];
			current_node.Last_point = current_n[0];
			current_node.pathlength = ctr_j;
			
			// The old triple	
			local_idx = MD5_tLSbits(3, position, OMEGA_MINUS_ONE, BITS_OF_OMEGA);
			
			//omp_set_lock(&writelock);	// To avoid reading an writing conflicts between threads
			saved_node = HASH_TABLE_FOR_THE_DISTINGUISHED_POINTS_REACHED[local_idx];
			
			// Do we have a collision?			
			if( (saved_node.pathlength > 0) && (current_node.Last_point == saved_node.Last_point) )
			{
				local_colls += 1;	// #(Collisions)
								
				FunEvals_col += saved_node.pathlength;		// #(Function evaluations for finding a collision)
				FunEvals_col += current_node.pathlength;	// #(Function evaluations for finding a collision)
				get_collision(local_COLLISION, saved_node, current_node);
				
				// We check if this collisions a new different from the previous saved				
				new_co[0] = local_COLLISION[0][0];
				new_co[1] = local_COLLISION[1][0];
				// --- We are saving the collision in lexicographic order
				if(new_co[0] > new_co[1])
				{
					local_tmp = new_co[0];
					new_co[0] = new_co[1];
					new_co[1] = local_tmp;
				}
				
				// We search for a possible repetition
				node_co = NULL;
				for(local_j = 0; local_j < NUMBER_OF_CORES; local_j++)
				{
					node_co = binarySearch(DIFFERENT_COLLISIONS_REACHED[local_j], new_co);
					if(node_co != NULL)
						break;
				}
				
				if(node_co == NULL)
				{
					diffe_colls += 1;	// #(Different collisions)
					node_co = &(ADDRESS_USED_FOR_THE_LIST_OF_COLLISIONS[kth_proc][CURRENT_SIZE_OF_THE_LIST_OF_COLLISIONS[kth_proc]]);
					createNode(node_co, new_co);
					sortedInsertion(&DIFFERENT_COLLISIONS_REACHED[kth_proc], node_co);
					CURRENT_SIZE_OF_THE_LIST_OF_COLLISIONS[kth_proc] += 1;
				}

				// Is it the golden collision?
				if( (local_COLLISION[0][0] & 0x04) != (local_COLLISION[1][0] & 0x04) )
				{
					// If so, we break (We've finished)
					GOLDEN_COLLISION[0][0] = local_COLLISION[0][0];
					GOLDEN_COLLISION[1][0] = local_COLLISION[1][0];
					*BOOL = 0x00;
					break;
				}
			}
			
			// We overwrite the distinguished point
			HASH_TABLE_FOR_THE_DISTINGUISHED_POINTS_REACHED[local_idx] = current_node;
			//omp_unset_lock(&writelock);	// To avoid reading an writing conflicts between threads
		}
		
	}
	local_end = get_cycles();
	
	// Logs: data about each version of f
	// --- Collisions	
	CURRENT_NUMBER_OF_TOTAL_COLLISIONS[kth_proc] = local_colls;
	CURRENT_SIZE_OF_THE_LIST_OF_COLLISIONS[kth_proc] = diffe_colls;
	
	// Running-time for distinguished points and collisions	
	FUNCTION_EVALUATIONS_FOR_GENERATING_DISTINGUISHED_POINTS[kth_proc] += FunEvals_pts;
	FUNCTION_EVALUATIONS_FOR_LOCATING_COLLISIONS[kth_proc] += FunEvals_col;
	
	return local_end - local_start;
};
