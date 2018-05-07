#include "naive.h"

int main (int argc, char *argv[]){

	if(argc == 3){
	
		// Variables for time execution
		time_t t;
		srand((unsigned) time(&t));
		clock_t start, end;
		
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
		Pt_DSM(&d_times_Pt, point_Pt, curve_E0.A);
		Pt_MDSM(&order_d_point_Qt, point_Qt,  e_DIVIDED_BY_2 - 1, curve_E0.A);
						
		// Computing d^(e/2)St and d^(e/2)Tt
		Pt_MDSM(&point_St, St,  e_DIVIDED_BY_2, curve_E1.A);
		Pt_MDSM(&point_Tt, Tt,  e_DIVIDED_BY_2, curve_E1.A);
		Pt_DSM(&d_times_St, point_St, curve_E1.A);
		Pt_MDSM(&order_d_point_Tt, point_Tt,  e_DIVIDED_BY_2 - 1, curve_E1.A);
		
		/* ---------------------------------------------------------- */
		
		VAR_N_DIVIDED_BY_DEGREE_PLUS_1 = (uint64_t)pow(DEGREE, e_DIVIDED_BY_2 - 1);
		VAR_N = VAR_N_DIVIDED_BY_DEGREE_PLUS_1 * (DEGREE + 1);
		
		// ----
		uint64_t LOCAL_LIST_SIZE = (uint64_t)ceil((double)VAR_N_DIVIDED_BY_DEGREE_PLUS_1 / (double)NUMBER_OF_CORES);
		TREE_LEAVES_ROOTED_AT_E0 = (T_naive **)malloc( NUMBER_OF_CORES * sizeof(T_naive *) );
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
			TREE_LEAVES_ROOTED_AT_E0[ctr_i] = (T_naive *)malloc( (DEGREE + 1) * LOCAL_LIST_SIZE * sizeof(T_naive) );
				
		RUNNING_TIME = (uint64_t *)malloc( sizeof(uint64_t) * NUMBER_OF_CORES );
		CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0 = (uint64_t *)malloc( sizeof(uint64_t) * NUMBER_OF_CORES );
		
		uint64_t TOTAL_CLOCK_CYCLES_BUILD = 0,
					TOTAL_CLOCK_CYCLES_SORT = 0,
					TOTAL_CLOCK_CYCLES_SEARCH = 0;
		// ----
		
		printf("\"We'll use %jd cores, and %jd  nodes for each core.\";\n\"We'll start with the hash table construction.\";\n", NUMBER_OF_CORES, LOCAL_LIST_SIZE);
		// We use OpenMP for parallel construction of the hash table.
		uint8_t STOP = 0x01;
		omp_set_num_threads(NUMBER_OF_CORES);
		#pragma omp parallel private(ctr_i) reduction(+:TOTAL_CLOCK_CYCLES_BUILD)
		{
			ctr_i = omp_get_thread_num();
			TOTAL_CLOCK_CYCLES_BUILD += Leaves_of_the_tree_rooted_at_E0(ctr_i);
		};
		
		printf("CLOCK_CYCLES_BUILD := %jd;\n", TOTAL_CLOCK_CYCLES_BUILD);
		printf("\"We'll start with the sorting of the hash table.\";\n");
		clock_t local_start, local_end;
		
		// We use OpenMP for parallel sorting of the hash table.
		// Until here, RUNNING_TIME is the size of each sub-list
		omp_set_num_threads(NUMBER_OF_CORES);
		#pragma omp parallel shared(TREE_LEAVES_ROOTED_AT_E0,RUNNING_TIME) private(ctr_i,local_start,local_end) reduction(+:TOTAL_CLOCK_CYCLES_SORT)
		{
			ctr_i = omp_get_thread_num();
			local_start = get_cycles();
			quickSort(TREE_LEAVES_ROOTED_AT_E0[ctr_i], 0, CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[ctr_i]);
			local_end = get_cycles();
			TOTAL_CLOCK_CYCLES_SORT += (local_end - local_start);
		};
		printf("CLOCK_CYCLES_SORT := %jd;\n", TOTAL_CLOCK_CYCLES_SORT);
		// We use OpenMP for parallel finding of the collision
		printf("\"We'll start with the collision search.\";\n");
		uint64_t SOLUTION[2][WORD_N];
		memcpy(SOLUTION[0], zeroM, sizeof(uint64_t) * WORD_N);
		memcpy(SOLUTION[1], zeroM, sizeof(uint64_t) * WORD_N);
		
		uint64_t acc_isog = 0;
		omp_set_num_threads(NUMBER_OF_CORES);
		#pragma omp parallel shared(SOLUTION,STOP) private(ctr_i) reduction(+:TOTAL_CLOCK_CYCLES_SEARCH)
		{
			ctr_i = omp_get_thread_num();
			TOTAL_CLOCK_CYCLES_SEARCH += Leaves_of_the_tree_rooted_at_E1(SOLUTION, &STOP, ctr_i);
		}
		
		printf("CLOCK_CYCLES_SEARCH := %jd;\n", TOTAL_CLOCK_CYCLES_SEARCH);
		uint64_t TOTAL_RUNNING_TIME = 0;
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
			TOTAL_RUNNING_TIME += RUNNING_TIME[ctr_i];
			
		printf("\nTOTAL_ISOGENY_COMPUTATIONS := %jd;\n", TOTAL_RUNNING_TIME);
		
		printf("clock_cycles := Sprintf(\"clock_cycles_naive.append(%%o)\", CLOCK_CYCLES_SEARCH + CLOCK_CYCLES_SORT + CLOCK_CYCLES_BUILD);\n");
		printf("running_time := Sprintf(\"running_time_naive.append(%%o)\", TOTAL_ISOGENY_COMPUTATIONS);\n");
		printf("PrintFile(\"./naive_results.py\", clock_cycles);\n");
		printf("PrintFile(\"./naive_results.py\", running_time);\n");
		// ----------------------------------------
		uint64_t COEF_one[WORD_N], COEF_two[WORD_N], Priv_COEF[WORD_N];
		memcpy(COEF_one, zeroM, sizeof(uint64_t) * WORD_N);
		COEF_one[0] = SOLUTION[0][0];
		memcpy(COEF_two, zeroM, sizeof(uint64_t) * WORD_N);
		COEF_two[0] = SOLUTION[1][0];
		
		// Recovery of the kernels from the left and right sides.
		Point Pt_0, Qt_0, St_1, Tt_1, Pt_copy, Qt_copy;

		// ------------------------------------------- //
		// ------------ Side of local_E -------------- //
		COEF_one[0] = SOLUTION[0][0] >> 2;
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
		else if( ((SOLUTION[0][0] & 0x03) == 0x02 ) && (DEGREE == 2) )
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
		COEF_two[0] = SOLUTION[1][0] >> 2;
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
		else if( ((SOLUTION[1][0] & 0x03) == 0x02 ) && (DEGREE == 2) )
		{
			printf("\tK_1 := (%jd * 0x%016jX) * (%jd^%jd) * St + (%jd^%jd) * Tt;\n", DEGREE, COEF_two[0], DEGREE, e_DIVIDED_BY_2, DEGREE, e_DIVIDED_BY_2);
			COEF_two[0] *= DEGREE;
			Pt_MDSM(&St_1, Tt, e_DIVIDED_BY_2, curve_E1.A);
			Pt_MDSM(&Tt_1, St, e_DIVIDED_BY_2, curve_E1.A);
		}
				
		
		// ----------------------------------------
		// Mapping the generators (isogeny evaluation)
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
		printf("exit;\n");

		// ---
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
			free(TREE_LEAVES_ROOTED_AT_E0[ctr_i]);
		
		free(TREE_LEAVES_ROOTED_AT_E0);
		free(RUNNING_TIME);
		free(CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0);
		
		return 0;
	}
	else{
		printf("Invalid case: It should take 2 args.\n\t\t 1) number of cores to be used, and\n\t\t 2) A (Alice) or B (Bob).\n" );
		return 1; 
	}
	
};
