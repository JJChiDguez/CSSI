#include "lambda.h"

int main (int argc, char *argv[])
{
	// INPUTS:
	// argv[1] <- number of cores to be used
	// argv[2] <- A (Alice) or B (bob)
	// argv[3] <- Log(2, w)
	// argv[4] <- \beta
	// argv[5] <- number of functions to be used
	
	if( (argc == 5) || (argc == 6) ){
		
		// Variables for time execution
		time_t t;
		srand((unsigned) time(&t));
		clock_t start, end;
		// ---
		printf("clear;\n\n");
		printf("\"SETTING TO BE USED IN THE ATTACK:\";\n\n");
		print_NUM(p, 2, 0, "p");
		print_NUM(base_mont, 2, 0, "R");

		printf("_, Rinv, _ := XGCD(R,p);\n");

		printf("\n");
		printf("F_p     := GF(p);\n");
		printf("P_p<t>  := PolynomialRing(F_p);\n");
		printf("F_q<i> := ext<F_p | t^2 + 1>;\n");
		printf("P_q<X> := PolynomialRing(F_q);\n\n");
		printf("load \"./magma/balanced_chain_of_d_isogenies.m\";\n\n");

		// ---------------------------------------------------------------------
		// ---------------------------------------------------------------------

		NUMBER_OF_CORES = strtol(argv[1], NULL, 10);
		char which_case = (char)argv[2][0];

		MAXIMUM_NUMBER_OF_FUNCTIONS_PERMITED = 0;
		if(argc == 6)
		{
			MAXIMUM_NUMBER_OF_FUNCTIONS_PERMITED = strtol(argv[5], NULL, 10);
			printf("\"We are going to run with %jd random functions!\";\n", MAXIMUM_NUMBER_OF_FUNCTIONS_PERMITED);
		}
				
		// Variable to be used for the construction of the hash table
		uint64_t ctr_i, var_i;
				
		// Points
		Point Pt, Qt, P, Q, phi_P, phi_Q, St, Tt;
		Point local_Pt, local_Qt,
				local_St, local_Tt;
				
		// Choosing the side to be attacked
		memcpy(curve_E0.A, E_0.A, sizeof(uint64_t) * WORD_N * 2);
		memcpy(curve_E0.B, E_0.B, sizeof(uint64_t) * WORD_N * 2);
		
		printf("\n");
		switch(which_case) {
			case 'A' :
				printf("\"WE'LL ATTACK FROM THE ALICE SIDE.\";\n\"THUS, WE'LL FIND THE 2^%d-ISOGENY phi : E -> E_isog\";\n\n", expn_a);
				Pt_DSM  = &Pt_DBL;
				Pt_MDSM = &Pt_MDBL;
				get_isog = &get_2_isog;
				eval_isog = &eval_2_isog;
				PohligHellman_alg = &DLP_order_two;

				Point_Assign(&Pt, P_a); Point_Assign(&Qt, Q_a);
				Point_Assign(&P, P_b); Point_Assign(&Q, Q_b);
				Point_Assign(&phi_P, phi_P_b); Point_Assign(&phi_Q, phi_Q_b);
				
				memcpy(curve_E1.A, E_A.A, sizeof(uint64_t) * WORD_N * 2);
				memcpy(curve_E1.B, E_A.B, sizeof(uint64_t) * WORD_N * 2);
				
				Point_Assign(&St, S_a);  Point_Assign(&Tt, T_a);
				
				e_DIVIDED_BY_2 = expn_a / 2; DEGREE = 2;
				break;
			case 'B' :
				printf("[error]: This case is not implemented.\n");
				printf("\t You should to extend the function h_i.\n");
				return 2;
				
				printf("\"WE'LL ATTACK FROM THE BOB SIDE.\";\n\"THUS, WE'LL FIND THE 3^%d-ISOGENY phi : E -> E_isog\";\n\n", expn_b);
				Pt_DSM  = &Pt_TRPL;
				Pt_MDSM = &Pt_MTRPL;
				get_isog = &get_3_isog;
				eval_isog = &eval_3_isog;
				PohligHellman_alg = &DLP_order_three;

				Point_Assign(&Pt, P_b); Point_Assign(&Qt, Q_b);
				Point_Assign(&P, P_a); Point_Assign(&Q, Q_a);
				Point_Assign(&phi_P, psi_P_a); Point_Assign(&phi_Q, psi_Q_a);
				
				memcpy(curve_E1.A, E_B.A, sizeof(uint64_t) * WORD_N * 2);
				memcpy(curve_E1.B, E_B.B, sizeof(uint64_t) * WORD_N * 2);
				
				Point_Assign(&St, S_b);  Point_Assign(&Tt, T_b);
				
				e_DIVIDED_BY_2 = expn_b / 2; DEGREE = 3;
				break;
			default :
				printf("Invalid case: It should be A (Alice) or B (Bob)!\n" );
				return 2;
		}

		/* ---------------------------------------------------------- */
		/* ---------------------------------------------------------- */

		printf("\"PUBLIC CURVES TO BE ATTACKED:\";\n\n");
		print_ELE( curve_E0.A, WORD_N, "E_coef_A");
		print_ELE( curve_E0.B, WORD_N, "E_coef_B");
		printf("E := EllipticCurve([F_q| E_coef_A, E_coef_B]);\n");
		
		printf("\n");
		print_ELE( curve_E1.A, WORD_N, "E_isog_coef_A");
		print_ELE( curve_E1.B, WORD_N, "E_isog_coef_B");
		printf("E_isog := EllipticCurve([F_q| E_isog_coef_A, E_isog_coef_B]);\n");
		
		
		printf("\n\"PUBLIC POINTS TO BE USED:\";\n\n");
		print_ELE( Pt.X, WORD_N, "Pt_X");
		print_ELE( Pt.Y, WORD_N, "Pt_Y");
		print_ELE( Pt.Z, WORD_N, "Pt_Z");
		printf("Pt := E![Pt_X / (Pt_Z^2), Pt_Y / (Pt_Z^3)];\n");
		
		print_ELE( Qt.X, WORD_N, "Qt_X");
		print_ELE( Qt.Y, WORD_N, "Qt_Y");
		print_ELE( Qt.Z, WORD_N, "Qt_Z");
		printf("Qt := E![Qt_X / (Qt_Z^2), Qt_Y / (Qt_Z^3)];\n");
		printf("\n");
		
		print_ELE( P.X, WORD_N, "P_X");
		print_ELE( P.Y, WORD_N, "P_Y");
		print_ELE( P.Z, WORD_N, "P_Z");
		printf("P := E![P_X / (P_Z^2), P_Y / (P_Z^3)];\n");
		
		print_ELE( Q.X, WORD_N, "Q_X");
		print_ELE( Q.Y, WORD_N, "Q_Y");
		print_ELE( Q.Z, WORD_N, "Q_Z");
		printf("Q := E![Q_X / (Q_Z^2), Q_Y / (Q_Z^3)];\n");
		
		printf("\n");

		print_ELE( phi_P.X, WORD_N, "phiP_X");
		print_ELE( phi_P.Y, WORD_N, "phiP_Y");
		print_ELE( phi_P.Z, WORD_N, "phiP_Z");
		printf("phiP := E_isog![phiP_X / (phiP_Z^2), phiP_Y / (phiP_Z^3)];\n");
		
		print_ELE( phi_Q.X, WORD_N, "phiQ_X");
		print_ELE( phi_Q.Y, WORD_N, "phiQ_Y");
		print_ELE( phi_Q.Z, WORD_N, "phiQ_Z");
		printf("phiQ := E_isog![phiQ_X / (phiQ_Z^2), phiQ_Y / (phiQ_Z^3)];\n");
		
		printf("\n\"A BASIS FOR E'[%jd^%jd] = <St, Tt>:\";\n\n", DEGREE, e_DIVIDED_BY_2 * 2);
		
		print_ELE( St.X, WORD_N, "St_X");
		print_ELE( St.Y, WORD_N, "St_Y");
		print_ELE( St.Z, WORD_N, "St_Z");
		printf("St := E_isog![St_X / (St_Z^2), St_Y / (St_Z^3)];\n");

		print_ELE( Tt.X, WORD_N, "Tt_X");
		print_ELE( Tt.Y, WORD_N, "Tt_Y");
		print_ELE( Tt.Z, WORD_N, "Tt_Z");
		printf("Tt := E_isog![Tt_X / (Tt_Z^2), Tt_Y / (Tt_Z^3)];\n");

		printf("\n\"SETUP OF THE FINDING CLAW PROBLEM:\";\n");
		printf("\"This program found a CLAW (X,Y) where\";\n");
		printf("\"\\tX is a [%jd^%jd]-order point on E,\";\n", DEGREE, e_DIVIDED_BY_2);
		printf("\"\\tY is a [%jd^%jd]-order point on E', and\";\n", DEGREE, e_DIVIDED_BY_2);
		printf("\"\\tE / <X> is isomorphic to E' / <Y>.\";\n");
		printf("\n");
		
		/* ---------------------------------------------------------- */
		// Reducing into Montgomery representation	
		// Side of E0
		
		MOD_MAPSTO_MONT(Pt.X); MOD_MAPSTO_MONT(Pt.Y); MOD_MAPSTO_MONT(Pt.Z);
		MOD_MAPSTO_MONT(Qt.X); MOD_MAPSTO_MONT(Qt.Y); MOD_MAPSTO_MONT(Qt.Z);
		MOD_MAPSTO_MONT(curve_E0.A); MOD_MAPSTO_MONT(curve_E0.B);
		
		// Side of E1
		MOD_MAPSTO_MONT(phi_P.X); MOD_MAPSTO_MONT(phi_P.Y); MOD_MAPSTO_MONT(phi_P.Z);
		MOD_MAPSTO_MONT(phi_Q.X); MOD_MAPSTO_MONT(phi_Q.Y); MOD_MAPSTO_MONT(phi_Q.Z);
		MOD_MAPSTO_MONT(St.X); MOD_MAPSTO_MONT(St.Y); MOD_MAPSTO_MONT(St.Z);
		MOD_MAPSTO_MONT(Tt.X); MOD_MAPSTO_MONT(Tt.Y); MOD_MAPSTO_MONT(Tt.Z);
		MOD_MAPSTO_MONT(curve_E1.A); MOD_MAPSTO_MONT(curve_E1.B);
		
		// Computing d^(e/2)Pt and d^(e/2)Qt
		Pt_MDSM(&point_Pt, Pt,  e_DIVIDED_BY_2, curve_E0.A);
		Pt_MDSM(&point_Qt, Qt,  e_DIVIDED_BY_2, curve_E0.A);
		Pt_MDSM(&order_d_point_Qt, point_Qt,  e_DIVIDED_BY_2 - 1, curve_E0.A);
						
		// Computing d^(e/2)St and d^(e/2)Tt
		Pt_MDSM(&point_St, St,  e_DIVIDED_BY_2, curve_E1.A);
		Pt_MDSM(&point_Tt, Tt,  e_DIVIDED_BY_2, curve_E1.A);
		Pt_MDSM(&order_d_point_Tt, point_Tt,  e_DIVIDED_BY_2 - 1, curve_E1.A);
		
		/* ---------------------------------------------------------- */
		
		// Number of bits for each "point" in the domain and image of FUNCTION_fn()
		BITS_OF_e_DIVIDED_BY_2_PLUS_2 = ceil(log(DEGREE + 1) / log(2)) + ceil(e_DIVIDED_BY_2 * log(DEGREE) / log(2));		
		VAR_M_MINUS_ONE = (uint64_t)(pow(2, BITS_OF_e_DIVIDED_BY_2_PLUS_2 ) - 1);
		
		BITS_OF_OMEGA  = strtol(argv[3], NULL, 10);
		OMEGA  = (uint64_t)pow(2.0, strtol(argv[3], NULL, 10));
		OMEGA_MINUS_ONE = (OMEGA - 1);
		
		BETA = strtol(argv[4], NULL, 10);
		BETA_TIMES_W = (uint64_t)(BETA * OMEGA);
		
		
		VAR_M = (uint64_t)((DEGREE + 1)*pow(DEGREE, e_DIVIDED_BY_2));
		THETA = 2.25 * sqrt((double)OMEGA / (double)VAR_M);	
		TRAIL_LENGTH = (uint64_t)ceil(1.0 / THETA);
		BITS_OF_R  = log(1.0 / THETA)/log(2);		
		DISTINGUISHABILITY = (uint64_t)ceil(pow(2, 32 - BITS_OF_R ));
		BETA_TIMES_W = (uint64_t)ceil((double)BETA_TIMES_W/(double)NUMBER_OF_CORES);
		
		printf("\"Parameters:\"; \n");
		printf("\" r : %lf, 2^32/2^r : %jd, w : 2^%jd : %jd, \\beta : %jd.\";\n", BITS_OF_R, DISTINGUISHABILITY, BITS_OF_OMEGA, OMEGA, BETA);
		printf("\"Each function will be used until it reaches (%jd x %jd) distinguished points,\";\n", NUMBER_OF_CORES, BETA_TIMES_W);
		printf("\"and the maximum trail lenght will be (10 x %jd).\";\n", TRAIL_LENGTH);
		printf("\"mask for w : 0x%016jX, Log(2,w) : %jd, Log(2, #DP) : %lf \";\n\n", OMEGA_MINUS_ONE, BITS_OF_OMEGA, log(VAR_M) / log(2) -  BITS_OF_R );
		
		
		T_lambda global_lists[OMEGA];	// Each thread will use an address of this ARRAY
		
		// MEMORY ALLOCATION FOR MEASURING THE RUNNING TIMES AND NUMBER OF COLLISIONS
		HASH_TABLE_FOR_THE_DISTINGUISHED_POINTS_REACHED = (T_lambda *)malloc( sizeof(T_lambda) * OMEGA);
		
		uint64_t LOCAL_LIST_SIZE = (uint64_t)ceil((double)(2*OMEGA) / (double)NUMBER_OF_CORES);
		ADDRESS_USED_FOR_THE_LIST_OF_COLLISIONS = (T_node **)malloc( NUMBER_OF_CORES * sizeof(T_node *) );
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
			ADDRESS_USED_FOR_THE_LIST_OF_COLLISIONS[ctr_i] = (T_node *)malloc( LOCAL_LIST_SIZE * sizeof(T_node) );
		
		DIFFERENT_COLLISIONS_REACHED = (T_node **)malloc( NUMBER_OF_CORES * sizeof(T_node *) );
		CURRENT_SIZE_OF_THE_LIST_OF_COLLISIONS = (uint64_t *)malloc( sizeof(uint64_t) * NUMBER_OF_CORES );
		CURRENT_NUMBER_OF_TOTAL_COLLISIONS = (uint64_t *)malloc( sizeof(uint64_t) * NUMBER_OF_CORES );
		
		FUNCTION_EVALUATIONS_FOR_GENERATING_DISTINGUISHED_POINTS = (uint64_t *)malloc( sizeof(uint64_t) * NUMBER_OF_CORES );
		FUNCTION_EVALUATIONS_FOR_LOCATING_COLLISIONS = (uint64_t *)malloc( sizeof(uint64_t) * NUMBER_OF_CORES );					
		
		// WE SET THE RUNNING TIMES TO ZERO
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
		{			
			FUNCTION_EVALUATIONS_FOR_GENERATING_DISTINGUISHED_POINTS[ctr_i] = 0;
			FUNCTION_EVALUATIONS_FOR_LOCATING_COLLISIONS[ctr_i] = 0;
		}
		
		
		// We used OpenMP for parallel construction and searh of the hash table.
		uint8_t STOP = 0x01;
		omp_set_num_threads(NUMBER_OF_CORES);
		
		uint64_t acc_CC = 0, total_CC = 0;
		
		printf("NONCES := [];\n");
		printf("seq_diff_collisions := [];\n");
		printf("seq_total_collision := [];\n\n");
		
		int count = 1, local_ctr;
		uint64_t SOLUTION[2][WORD_N];
		memcpy(SOLUTION[0], zeroM, sizeof(uint64_t) * WORD_N);
		memcpy(SOLUTION[1], zeroM, sizeof(uint64_t) * WORD_N);
			
		uint64_t TOTAL_DIFFERENT_COLLISIONS = 0,
					TOTAL_NUMBER_OF_COLLISIONS = 0;
		
		uint64_t TOTAL_FUN_EVALS_FOR_POINTS = 0,
					TOTAL_FUN_EVALS_FOR_COLLIS = 0;
		
		uint64_t isog_required[NUMBER_OF_CORES];
		
		//omp_init_lock(&writelock);
		while(STOP == 0x01)
		{
			// WE RANDOMLY SELECTED A NONCE FROM {0,1}^64
			NONCE = ((uint64_t)rand() & 0x00000000FFFFFFFF) ^ ( ((uint64_t)rand() & 0x00000000FFFFFFFF) << 32);
			
			// EACH ENTRY OF THE HASH TABLE IS INITIALIZATES TO ZERO
			for(ctr_i = 0; ctr_i < OMEGA; ctr_i++)
			{
				HASH_TABLE_FOR_THE_DISTINGUISHED_POINTS_REACHED[ctr_i].Seed_point = 0x0000000000000000;
				HASH_TABLE_FOR_THE_DISTINGUISHED_POINTS_REACHED[ctr_i].Last_point = 0x0000000000000000;
				HASH_TABLE_FOR_THE_DISTINGUISHED_POINTS_REACHED[ctr_i].pathlength = 0;
			}
			
			// WE INITIALIZES THE NUMBER OF COLLISIONS TO ZERO
			for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
			{
				// Collisions
				DIFFERENT_COLLISIONS_REACHED[ctr_i] = NULL;
				CURRENT_SIZE_OF_THE_LIST_OF_COLLISIONS[ctr_i] = 0;
				CURRENT_NUMBER_OF_TOTAL_COLLISIONS[ctr_i] = 0;
			}
			
			// van Oorschot-Wiener algorithm
			#pragma omp parallel shared(SOLUTION,STOP) private(ctr_i)
			{
				ctr_i = omp_get_thread_num();
				isog_required[ctr_i] = VW_algorithm(SOLUTION, &STOP, ctr_i);
			}
			
			TOTAL_DIFFERENT_COLLISIONS = 0;
			TOTAL_NUMBER_OF_COLLISIONS = 0;
			
			acc_CC = 0;
			//#pragma omp parallel for reduction(+:acc_CC)
			for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
			{
				// Clock Cycles and Running-time for the current version of the function fn
				acc_CC += isog_required[ctr_i];
				// Collisions
				TOTAL_DIFFERENT_COLLISIONS += CURRENT_SIZE_OF_THE_LIST_OF_COLLISIONS[ctr_i];
				TOTAL_NUMBER_OF_COLLISIONS += CURRENT_NUMBER_OF_TOTAL_COLLISIONS[ctr_i];
			}
			total_CC += acc_CC;
			
			// Distinguished points and Collisions
			printf("NONCES[%05d] := 0x%016jX; seq_diff_collisions[%05d] := %jd; seq_total_collision[%05d] := %jd;\n", count, NONCE, count, TOTAL_DIFFERENT_COLLISIONS, count, TOTAL_NUMBER_OF_COLLISIONS);
		
			count += 1;
			if( (count > MAXIMUM_NUMBER_OF_FUNCTIONS_PERMITED) && (MAXIMUM_NUMBER_OF_FUNCTIONS_PERMITED > 0) )
				STOP = 0x00;
		}
		//omp_destroy_lock(&writelock);
		
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
		{
			TOTAL_FUN_EVALS_FOR_POINTS += FUNCTION_EVALUATIONS_FOR_GENERATING_DISTINGUISHED_POINTS[ctr_i];
			TOTAL_FUN_EVALS_FOR_COLLIS += FUNCTION_EVALUATIONS_FOR_LOCATING_COLLISIONS[ctr_i];
		}
		
		
		printf("DIFFERENT_FUNCTIONS_USED := %d;\n", count - 1);
		printf("TOTAL_FUN_EVALS_FOR_POINTS := %jd;\n", TOTAL_FUN_EVALS_FOR_POINTS);
		printf("TOTAL_FUN_EVALS_FOR_COLLIS := %jd;\n", TOTAL_FUN_EVALS_FOR_COLLIS);
		printf("TOTAL_FUNCTION_EVALUATIONS := %jd;\n", TOTAL_FUN_EVALS_FOR_POINTS + TOTAL_FUN_EVALS_FOR_COLLIS);
		printf("CLOCK_CYCLES_SEARCH := %jd;\n\n", total_CC);
		
		printf("dist_pts_run := Sprintf(\"running_time_dist_pnts.append(%%o)\", TOTAL_FUN_EVALS_FOR_POINTS);\n");
		printf("cols_running := Sprintf(\"running_time_collision.append(%%o)\", TOTAL_FUN_EVALS_FOR_COLLIS);\n");
		
		printf("diffunctions := Sprintf(\"diffunctions_lambda.append(%%o)\", DIFFERENT_FUNCTIONS_USED);\n");
		printf("running_time := Sprintf(\"running_time_lambda.append(%%o)\", TOTAL_FUNCTION_EVALUATIONS);\n");
		printf("clock_cycles := Sprintf(\"clock_cycles_lambda.append(%%o)\", CLOCK_CYCLES_SEARCH);\n");
		
		printf("PrintFile(\"./lambda_results.py\", dist_pts_run);\n");	// Running-time for generating all distinguished points (for all functions)
		printf("PrintFile(\"./lambda_results.py\", cols_running);\n");	// Running-time for locating all collisions (for all functions)
		
		printf("PrintFile(\"./lambda_results.py\", diffunctions);\n");	// # of version of the function fn
		printf("PrintFile(\"./lambda_results.py\", running_time);\n");	// Total running-time
		printf("PrintFile(\"./lambda_results.py\", clock_cycles);\n");	// Clock Cycles
		
		// Collisions
		printf("diff_collisions := Sprintf(\"c.append(%%o)\", seq_diff_collisions);\n");
		printf("total_collision := Sprintf(\"t.append(%%o)\", seq_total_collision);\n");
		
		printf("PrintFile(\"./lambda_results.py\", diff_collisions);\n");	// # Different collisions
		printf("PrintFile(\"./lambda_results.py\", total_collision);\n");	// # collisions
						
		// If the goal was to find the golden collision and wasn't to do experiments
		if(MAXIMUM_NUMBER_OF_FUNCTIONS_PERMITED == 0)
		{
			// ----------------------------------------
			uint64_t COEF_one[WORD_N], COEF_two[WORD_N], Priv_COEF[WORD_N];
			memcpy(COEF_one, zeroM, sizeof(uint64_t) * WORD_N);
			memcpy(COEF_two, zeroM, sizeof(uint64_t) * WORD_N);
			
			Point Pt_0, Qt_0, St_1, Tt_1, Pt_copy, Qt_copy;
			
			if( (SOLUTION[0][0] & 0x04) == 0x00 )
			{
				// ------------------------------------------- //
				// ------------ Side of local_E -------------- //
				COEF_one[0] = SOLUTION[0][0] >> 3;
				printf("\n\"kernel from E:\";\n");
				//if( (SOLUTION.sideANDform & 0x03) == 0x00 )
				if( (SOLUTION[0][0] & 0x03) == 0x00 )
				{
					printf("\tK_0 := (%jd^%jd) * Pt + 0x%016jX * (%jd^%jd) * Qt;\n", DEGREE, e_DIVIDED_BY_2, COEF_one[0], DEGREE, e_DIVIDED_BY_2);
					// The points to evaluate for solving the DLP (used to get the kernel of the d^e-isogeny)
					Point_Assign(&Pt_copy, Pt);
					Point_Assign(&Qt_copy, Qt);
					// Generators for the kernel
					Pt_MDSM(&Pt_0, Pt, e_DIVIDED_BY_2, curve_E0.A);
					Pt_MDSM(&Qt_0, Qt, e_DIVIDED_BY_2, curve_E0.A);
				}
				else if( (SOLUTION[0][0] & 0x03) == 0x01 )
				{
					printf("\tK_0 := (%jd^%jd) * Pt + (%jd^%jd + 0x%016jX) * (%jd^%jd) * Qt;\n", DEGREE, e_DIVIDED_BY_2, DEGREE, e_DIVIDED_BY_2 - 1, COEF_one[0], DEGREE, e_DIVIDED_BY_2);
					COEF_one[0] += (uint64_t)pow(DEGREE, e_DIVIDED_BY_2 - 1);
					// The points to evaluate for solving the DLP (used to get the kernel of the d^e-isogeny)
					Point_Assign(&Pt_copy, Pt);
					Point_Assign(&Qt_copy, Qt);
					// Generators for the kernel
					Pt_MDSM(&Pt_0, Pt, e_DIVIDED_BY_2, curve_E0.A);
					Pt_MDSM(&Qt_0, Qt, e_DIVIDED_BY_2, curve_E0.A);
				}			
				else if( ( (SOLUTION[0][0] & 0x03) == 0x02 ) && (DEGREE == 2) )
				{
					printf("\tK_0 := (%jd * 0x%016jX) * (%jd^%jd) * Pt + (%jd^%jd) * Qt;\n", DEGREE, COEF_one[0], DEGREE, e_DIVIDED_BY_2, DEGREE, e_DIVIDED_BY_2);
					COEF_one[0] *= DEGREE;
					// The points to evaluate for solving the DLP (used to get the kernel of the d^e-isogeny)
					Point_Assign(&Pt_copy, Qt);
					Point_Assign(&Qt_copy, Pt);
					// Generators for the kernel
					Pt_MDSM(&Pt_0, Qt, e_DIVIDED_BY_2, curve_E0.A);
					Pt_MDSM(&Qt_0, Pt, e_DIVIDED_BY_2, curve_E0.A);
				}
				// ------------------------------------------- //
				// ---------- Side of local_E_isog ----------- //
				COEF_two[0] = SOLUTION[1][0] >> 3;
				printf("\n\"kernel from E_isog:\";\n");
				if( (SOLUTION[1][0] & 0x03) == 0x00 )
				{
					printf("\tK_1 := (%jd^%jd) * St + 0x%016jX * (%jd^%jd) * Tt;\n", DEGREE, e_DIVIDED_BY_2, COEF_two[0], DEGREE, e_DIVIDED_BY_2);
					Pt_MDSM(&St_1, St, e_DIVIDED_BY_2, curve_E1.A);
					Pt_MDSM(&Tt_1, Tt, e_DIVIDED_BY_2, curve_E1.A);
				}
				else if( (SOLUTION[1][0] & 0x03) == 0x01 )
				{
					printf("\tK_1 := (%jd^%jd) * St + (%jd^%jd + 0x%016jX) * (%jd^%jd) * Tt;\n", DEGREE, e_DIVIDED_BY_2, DEGREE, e_DIVIDED_BY_2 - 1, COEF_two[0], DEGREE, e_DIVIDED_BY_2);
					COEF_two[0] += (uint64_t)pow(DEGREE, e_DIVIDED_BY_2 - 1);
					Pt_MDSM(&St_1, St, e_DIVIDED_BY_2, curve_E1.A);
					Pt_MDSM(&Tt_1, Tt, e_DIVIDED_BY_2, curve_E1.A);
				}
				else if( ( (SOLUTION[1][0] & 0x03) == 0x02 ) && (DEGREE == 2) )
				{
					printf("\tK_1 := (%jd * 0x%016jX) * (%jd^%jd) * St + (%jd^%jd) * Tt;\n", DEGREE, COEF_two[0], DEGREE, e_DIVIDED_BY_2, DEGREE, e_DIVIDED_BY_2);
					COEF_two[0] *= DEGREE;
					Pt_MDSM(&St_1, Tt, e_DIVIDED_BY_2, curve_E1.A);
					Pt_MDSM(&Tt_1, St, e_DIVIDED_BY_2, curve_E1.A);
				}
			}
			else
			{
				// ------------------------------------------- //
				// ------------ Side of local_E -------------- //
				COEF_one[0] = SOLUTION[1][0] >> 3;
				printf("\n\"kernel from E:\";\n");
				if( (SOLUTION[1][0] & 0x03) == 0x00 )
				{
					printf("\tK_0 := (%jd^%jd) * Pt + 0x%016jX * (%jd^%jd) * Qt;\n", DEGREE, e_DIVIDED_BY_2, COEF_one[0], DEGREE, e_DIVIDED_BY_2);
					// The points to evaluate for solving the DLP (used to get the kernel of the d^e-isogeny)
					Point_Assign(&Pt_copy, Pt);
					Point_Assign(&Qt_copy, Qt);
					// Generators for the kernel
					Pt_MDSM(&Pt_0, Pt, e_DIVIDED_BY_2, curve_E0.A);
					Pt_MDSM(&Qt_0, Qt, e_DIVIDED_BY_2, curve_E0.A);
				}
				else if( (SOLUTION[1][0] & 0x03) == 0x01 )
				{
					printf("\tK_0 := (%jd^%jd) * Pt + (%jd^%jd + 0x%016jX) * (%jd^%jd) * Qt;\n", DEGREE, e_DIVIDED_BY_2, DEGREE, e_DIVIDED_BY_2 - 1, COEF_one[0], DEGREE, e_DIVIDED_BY_2);
					COEF_one[0] += (uint64_t)pow(DEGREE, e_DIVIDED_BY_2 - 1);
					// The points to evaluate for solving the DLP (used to get the kernel of the d^e-isogeny)
					Point_Assign(&Pt_copy, Pt);
					Point_Assign(&Qt_copy, Qt);
					// Generators for the kernel
					Pt_MDSM(&Pt_0, Pt, e_DIVIDED_BY_2, curve_E0.A);
					Pt_MDSM(&Qt_0, Qt, e_DIVIDED_BY_2, curve_E0.A);
				}			
				else if( ( (SOLUTION[1][0] & 0x03) == 0x02 ) && (DEGREE == 2) )
				{
					printf("\tK_0 := (%jd * 0x%016jX) * (%jd^%jd) * Pt + (%jd^%jd) * Qt;\n", DEGREE, COEF_one[0], DEGREE, e_DIVIDED_BY_2, DEGREE, e_DIVIDED_BY_2);
					COEF_one[0] *= DEGREE;
					// The points to evaluate for solving the DLP (used to get the kernel of the d^e-isogeny)
					Point_Assign(&Pt_copy, Qt);
					Point_Assign(&Qt_copy, Pt);
					// Generators for the kernel
					Pt_MDSM(&Pt_0, Qt, e_DIVIDED_BY_2, curve_E0.A);
					Pt_MDSM(&Qt_0, Pt, e_DIVIDED_BY_2, curve_E0.A);
				}
				// ------------------------------------------- //
				// ---------- Side of local_E_isog ----------- //
				COEF_two[0] = SOLUTION[0][0] >> 3;
				printf("\n\"kernel from E_isog:\";\n");
				if( (SOLUTION[0][0] & 0x03) == 0x00 )
				{
					printf("\tK_1 := (%jd^%jd) * St + 0x%016jX * (%jd^%jd) * Tt;\n", DEGREE, e_DIVIDED_BY_2, COEF_two[0], DEGREE, e_DIVIDED_BY_2);
					Pt_MDSM(&St_1, St, e_DIVIDED_BY_2, curve_E1.A);
					Pt_MDSM(&Tt_1, Tt, e_DIVIDED_BY_2, curve_E1.A);
				}
				else if( (SOLUTION[0][0] & 0x03) == 0x01 )
				{
					printf("\tK_1 := (%jd^%jd) * St + (%jd^%jd + 0x%016jX) * (%jd^%jd) * Tt;\n", DEGREE, e_DIVIDED_BY_2, DEGREE, e_DIVIDED_BY_2 - 1, COEF_two[0], DEGREE, e_DIVIDED_BY_2);
					COEF_two[0] += (uint64_t)pow(DEGREE, e_DIVIDED_BY_2 - 1);
					Pt_MDSM(&St_1, St, e_DIVIDED_BY_2, curve_E1.A);
					Pt_MDSM(&Tt_1, Tt, e_DIVIDED_BY_2, curve_E1.A);
				}
				else if( ( (SOLUTION[0][0] & 0x03) == 0x02 ) && (DEGREE == 2) )
				{
					printf("\tK_1 := (%jd * 0x%016jX) * (%jd^%jd) * St + (%jd^%jd) * Tt;\n", DEGREE, COEF_two[0], DEGREE, e_DIVIDED_BY_2, DEGREE, e_DIVIDED_BY_2);
					COEF_two[0] *= DEGREE;
					Pt_MDSM(&St_1, Tt, e_DIVIDED_BY_2, curve_E1.A);
					Pt_MDSM(&Tt_1, St, e_DIVIDED_BY_2, curve_E1.A);
				}			
			}

			Curve C_one, C_two, C_three;
			Point K_one, K_two;
			Isomorphism iota;

			// ------
			DBLADD(&K_one, Qt_0, COEF_one, curve_E0.A);	// K_one <- [COEF_one]Qt_0
			Pt_ADD(&K_one, Pt_0,    K_one, curve_E0.A);	// K_one <- Pt_0 + [COEF_one]Qt_0
			get_isogenous_curve(&C_one, K_one, Pt_MDSM, get_isog, eval_isog, curve_E0, e_DIVIDED_BY_2);
			eval_isogeny(&Pt_copy,  &Qt_copy, K_one, Pt_MDSM, get_isog, eval_isog, curve_E0, e_DIVIDED_BY_2);

			// ------
			DBLADD(&K_two, Tt_1, COEF_two, curve_E1.A);	// K_two <- [COEF_two]Tt_1
			Pt_ADD(&K_two, St_1,    K_two, curve_E1.A);	// K_two <- St_1 + [COEF_two]Tt_1
			get_isogenous_curve(&C_two, K_two, Pt_MDSM, get_isog, eval_isog, curve_E1, e_DIVIDED_BY_2);
			eval_isogeny(&St_1,  &Tt_1, K_two, Pt_MDSM, get_isog, eval_isog, curve_E1, e_DIVIDED_BY_2);
			
			// ------
			get_isomorphism(&iota, C_two, C_one);
			eval_isomorphism(&Tt_1, iota, Tt_1);
			
			get_isogenous_curve(&C_three, Tt_1, Pt_MDSM, get_isog, eval_isog, C_one, e_DIVIDED_BY_2);
			eval_isogeny(&Pt_copy, &Qt_copy, Tt_1, Pt_MDSM, get_isog, eval_isog, C_one, e_DIVIDED_BY_2);
			
			// ------
			PohligHellman_alg(Priv_COEF, Pt_copy, Qt_copy, e_DIVIDED_BY_2 * 2, C_three.A);
			print_NUM(Priv_COEF, 2, 0, "Priv_COEF");

			printf("C_0 := d_pow_e_isogenous_curve(K_0, E, %jd, %jd);\n", DEGREE, e_DIVIDED_BY_2);
			printf("C_1 := d_pow_e_isogenous_curve(K_1, E_isog, %jd, %jd);\n", DEGREE, e_DIVIDED_BY_2);
			printf("jInvariant(C_0) eq jInvariant(C_1);\n");
		}
		
		// ---
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
			free(ADDRESS_USED_FOR_THE_LIST_OF_COLLISIONS[ctr_i]);
		free(ADDRESS_USED_FOR_THE_LIST_OF_COLLISIONS);
		
		free(DIFFERENT_COLLISIONS_REACHED);	
		free(CURRENT_SIZE_OF_THE_LIST_OF_COLLISIONS);
		free(CURRENT_NUMBER_OF_TOTAL_COLLISIONS);
		
		free(FUNCTION_EVALUATIONS_FOR_GENERATING_DISTINGUISHED_POINTS);
		free(FUNCTION_EVALUATIONS_FOR_LOCATING_COLLISIONS);
		
		free(HASH_TABLE_FOR_THE_DISTINGUISHED_POINTS_REACHED);
		printf("exit;\n");
		return 0;
	}
	else{
		printf("[error]: not enough number of arguments. This code requires:\n");
		printf("\t 1) number of cores to be used,\n");
		printf("\t 2) A (Alice) or B (Bob),\n");
		printf("\t 3) Log(2, w) : Log(2, number of elements in memory that we can use), \n");
		printf("\t 4) beta where beta * w is the maximum number of distinguished points to be reached, and\n");
		printf("\t 5) number of functions to be used.\n");
		return 1;
	}

};
