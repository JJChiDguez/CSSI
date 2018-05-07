#include "dfs.h"

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

		LEVEL = strtol(argv[1], NULL, 10);
		NUMBER_OF_CORES = 3 * (uint64_t)pow(2, LEVEL);
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
						
		// Computing d^(e/2)St and d^(e/2)Tt
		Pt_MDSM(&point_St, St,  e_DIVIDED_BY_2, curve_E1.A);
		Pt_MDSM(&point_Tt, Tt,  e_DIVIDED_BY_2, curve_E1.A);
		
		/* ---------------------------------------------------------- */
		if( LEVEL > e_DIVIDED_BY_2 )
		{
			printf("[error]: The number of cores : 3*2^(%jd) is greater than the number of 2^(%jd)-isogenies : 3*2^(%jd)\n", LEVEL, e_DIVIDED_BY_2, e_DIVIDED_BY_2);
			return 3;
		}
			
		e_DIVIDED_BY_2_MINUS_LEVEL = e_DIVIDED_BY_2 - LEVEL;
				
		VAR_N_DIVIDED_BY_DEGREE_PLUS_1 = (uint64_t)pow(DEGREE, e_DIVIDED_BY_2 - 1);
		VAR_N = VAR_N_DIVIDED_BY_DEGREE_PLUS_1 * (DEGREE + 1);
		
		uint64_t LOCAL_LIST_SIZE = (uint64_t)pow(2, e_DIVIDED_BY_2_MINUS_LEVEL - 1);
		TREE_LEAVES_ROOTED_AT_E0 = (T_naive **)malloc( NUMBER_OF_CORES * sizeof(T_naive *) );
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
			TREE_LEAVES_ROOTED_AT_E0[ctr_i] = (T_naive *)malloc( (DEGREE + 1) * LOCAL_LIST_SIZE * sizeof(T_naive) );
		
		CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0 = (uint64_t *)malloc( sizeof(uint64_t) * NUMBER_OF_CORES );
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
			CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[ctr_i] = 0;
			
		CURRENT_DOUBLINGS = (uint64_t *)malloc( sizeof(uint64_t) * NUMBER_OF_CORES );
		CURRENT_ADDITIONS = (uint64_t *)malloc( sizeof(uint64_t) * NUMBER_OF_CORES );
		CURRENT_ISOG_EVAL = (uint64_t *)malloc( sizeof(uint64_t) * NUMBER_OF_CORES );
		CURRENT_ISOG_COMP = (uint64_t *)malloc( sizeof(uint64_t) * NUMBER_OF_CORES );
		
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
		{
			CURRENT_DOUBLINGS[ctr_i] = 0;
			CURRENT_ADDITIONS[ctr_i] = 0;
			CURRENT_ISOG_COMP[ctr_i] = 0;
			CURRENT_ISOG_EVAL[ctr_i] = 0;
		}
		
		uint64_t TOTAL_DOUBLINGS = 0,
					TOTAL_ADDITIONS = 0,
					TOTAL_ISOG_COMP = 0,
					TOTAL_ISOG_EVAL = 0;
		
		uint64_t TOTAL_CLOCK_CYCLES_BUILD = 0,
					TOTAL_CLOCK_CYCLES_SORT = 0,
					TOTAL_CLOCK_CYCLES_SEARCH = 0;
							
		// -------
		// -------
		uint64_t Xs_AT_E0[NUMBER_OF_CORES];
		T_dfs NODES_AT_E0[NUMBER_OF_CORES];
      
      		// [2]point_Pt
      Pt_DSM(&KERNELS_AT_E0[0], point_Pt, curve_E0.A);
      // [2]point_Qt
      Pt_DSM(&KERNELS_AT_E0[1], point_Qt, curve_E0.A);
      // point_Pt + point_Qt
      Pt_ADD(&KERNELS_AT_E0[2], point_Pt, point_Qt, curve_E0.A);
      
      // ----- CASE RELATED WITH THE POINT point_Pt
      Pt_MDSM(&KERNEL_AT_E0[0], KERNELS_AT_E0[0], e_DIVIDED_BY_2 - 2, curve_E0.A);
      // We compute the isogeny and the isogenous curve
      get_isog(&(INITIAL_NODES_AT_E0[0].Ei), &PHI_AT_E0[0], KERNEL_AT_E0[0], curve_E0);
      // phi(point_Pt)
      eval_isog(&(INITIAL_NODES_AT_E0[0].Pi), PHI_AT_E0[0], point_Pt);
      // phi([2]point_Qt)
      eval_isog(&(INITIAL_NODES_AT_E0[0].Qi), PHI_AT_E0[0], KERNELS_AT_E0[1]);
		
      // ----- CASE RELATED WITH THE POINT point_Qt
      Pt_MDSM(&KERNEL_AT_E0[1], KERNELS_AT_E0[1], e_DIVIDED_BY_2 - 2, curve_E0.A);
      // We compute the isogeny and the isogenous curve
      get_isog(&(INITIAL_NODES_AT_E0[1].Ei), &PHI_AT_E0[1], KERNEL_AT_E0[1], curve_E0);
      // phi([2]point_Pt)
      eval_isog(&(INITIAL_NODES_AT_E0[1].Qi), PHI_AT_E0[1], KERNELS_AT_E0[0]);
      // phi(point_Qt)
      eval_isog(&(INITIAL_NODES_AT_E0[1].Pi), PHI_AT_E0[1], point_Qt);
		
      // ----- CASE RELATED WITH THE POINT point_Pt + point_Qt
      // The Y coordinate of Px + Qx is equals zero
		memcpy(KERNEL_AT_E0[2].Y[0], zeroM, sizeof(uint64_t) * WORD_N);
		memcpy(KERNEL_AT_E0[2].Y[1], zeroM, sizeof(uint64_t) * WORD_N);
		// The Z coordinate of Px + Qx is equals Pi.Z*Qi.Z
		MUL_Fp2(KERNEL_AT_E0[2].Z, KERNEL_AT_E0[0].Z, KERNEL_AT_E0[1].Z);
		// The X coordinate of Px + Qx is equals Pi.X*(Qi.Z)^2 + (Pi.Z^2)*Qi.X
		uint64_t XpZqZq[2][WORD_N], XqZpZp[2][WORD_N], ZpZp[2][WORD_N], ZqZq[2][WORD_N];
		SQR_Fp2(ZqZq, KERNEL_AT_E0[1].Z);
		MUL_Fp2(XpZqZq, KERNEL_AT_E0[0].X, ZqZq);
		SQR_Fp2(ZpZp, KERNEL_AT_E0[0].Z);
		MUL_Fp2(XqZpZp, KERNEL_AT_E0[1].X, ZpZp);
		ADD_Fp2(KERNEL_AT_E0[2].X, XpZqZq, XqZpZp);
		NEG_Fp2(KERNEL_AT_E0[2].X, KERNEL_AT_E0[2].X);
				
      // We compute the isogeny and the isogenous curve
      get_isog(&(INITIAL_NODES_AT_E0[2].Ei), &PHI_AT_E0[2], KERNEL_AT_E0[2], curve_E0);
      // phi([2]point_Pt)
      eval_isog(&(INITIAL_NODES_AT_E0[2].Qi), PHI_AT_E0[2], KERNELS_AT_E0[0]);
      // phi(point_Pt + point_Qt)
      eval_isog(&(INITIAL_NODES_AT_E0[2].Pi), PHI_AT_E0[2], KERNELS_AT_E0[2]);		
      
      // ---
      TOTAL_ADDITIONS += 1;
      TOTAL_ISOG_COMP += 3;
      TOTAL_ISOG_EVAL += 6;
      TOTAL_DOUBLINGS += (e_DIVIDED_BY_2 - 1) * 2;
      
		// -------
		// -------
		
		omp_set_num_threads(3);
		#pragma omp parallel private(ctr_i) shared(NODES_AT_E0,Xs_AT_E0)
		{
			ctr_i = omp_get_thread_num();
			build_initial_nodes(NODES_AT_E0, Xs_AT_E0, INITIAL_NODES_AT_E0[ctr_i], ctr_i, 1, ctr_i);
		};
				
		printf("\"We'll use %jd cores, and %jd  nodes for each core.\";\n\"We'll start with the hash table construction.\";\n", NUMBER_OF_CORES, LOCAL_LIST_SIZE);
		// We use OpenMP for parallel construction of the hash table.
		clock_t local_start, local_end;
				
		omp_set_num_threads(NUMBER_OF_CORES);
		#pragma omp parallel shared(NODES_AT_E0, Xs_AT_E0) private(ctr_i,local_start,local_end) reduction(+:TOTAL_CLOCK_CYCLES_BUILD)
		{
			ctr_i = omp_get_thread_num();
			local_start = get_cycles();
			Depth_First_Search_at_E0(NODES_AT_E0[ctr_i], Xs_AT_E0[ctr_i], LEVEL + 1, ctr_i);
			local_end = get_cycles();
			TOTAL_CLOCK_CYCLES_BUILD += (local_end - local_start);
		};
		
		printf("CLOCK_CYCLES_BUILD := %jd;\n", TOTAL_CLOCK_CYCLES_BUILD);
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
		{
			TOTAL_DOUBLINGS += CURRENT_DOUBLINGS[ctr_i],
			TOTAL_ADDITIONS += CURRENT_ADDITIONS[ctr_i],
			TOTAL_ISOG_COMP += CURRENT_ISOG_COMP[ctr_i];
			TOTAL_ISOG_EVAL += CURRENT_ISOG_EVAL[ctr_i];
		}
		printf("DBLS_AT_E0 := %jd;\tADDS_AT_E0 := %jd;\tCOMP_AT_E0 := %jd;\tEVAL_AT_E0 := %jd;\n", TOTAL_DOUBLINGS, TOTAL_ADDITIONS, TOTAL_ISOG_COMP,TOTAL_ISOG_EVAL);
		printf("\"We'll start with the sorting of the hash table.\";\n");
		
		// We use OpenMP for parallel sorting of the hash table.
		// Until here, RUNNING_TIME is the size of each sub-list
		omp_set_num_threads(NUMBER_OF_CORES);
		#pragma omp parallel shared(TREE_LEAVES_ROOTED_AT_E0) private(ctr_i,local_start,local_end) reduction(+:TOTAL_CLOCK_CYCLES_SORT)
		{
			ctr_i = omp_get_thread_num();
			local_start = get_cycles();
			quickSort(TREE_LEAVES_ROOTED_AT_E0[ctr_i], 0, CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[ctr_i]);
			local_end = get_cycles();
			TOTAL_CLOCK_CYCLES_SORT += (local_end - local_start);
		};
		printf("CLOCK_CYCLES_SORT := %jd;\n", TOTAL_CLOCK_CYCLES_SORT);
		
		
		// -------
		// -------
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
		{
			CURRENT_DOUBLINGS[ctr_i] = 0;
			CURRENT_ADDITIONS[ctr_i] = 0;
			CURRENT_ISOG_COMP[ctr_i] = 0;
			CURRENT_ISOG_EVAL[ctr_i] = 0;
		}
		
		TOTAL_DOUBLINGS = 0,
		TOTAL_ADDITIONS = 0,
		TOTAL_ISOG_COMP = 0,
		TOTAL_ISOG_EVAL = 0;
					
		uint64_t Xs_AT_E1[NUMBER_OF_CORES];
		T_dfs NODES_AT_E1[NUMBER_OF_CORES];
       
 		// [2]point_St
      Pt_DSM(&KERNELS_AT_E1[0], point_St, curve_E1.A);
      // [2]point_Tt
      Pt_DSM(&KERNELS_AT_E1[1], point_Tt, curve_E1.A);
      // point_St + point_Tt
      Pt_ADD(&KERNELS_AT_E1[2], point_St, point_Tt, curve_E1.A);
      
      // ----- CASE RELATED WITH THE POINT point_St
      Pt_MDSM(&KERNEL_AT_E1[0], KERNELS_AT_E1[0], e_DIVIDED_BY_2 - 2, curve_E1.A);
      // We compute the isogeny and the isogenous curve
      get_isog(&(INITIAL_NODES_AT_E1[0].Ei), &PSI_AT_E1[0], KERNEL_AT_E1[0], curve_E1);
      // phi(point_St)
      eval_isog(&(INITIAL_NODES_AT_E1[0].Pi), PSI_AT_E1[0], point_St);
      // phi([2]point_Tt)
      eval_isog(&(INITIAL_NODES_AT_E1[0].Qi), PSI_AT_E1[0], KERNELS_AT_E1[1]);
		
      // ----- CASE RELATED WITH THE POINT point_Qt
      Pt_MDSM(&KERNEL_AT_E1[1], KERNELS_AT_E1[1], e_DIVIDED_BY_2 - 2, curve_E1.A);
      // We compute the isogeny and the isogenous curve
      get_isog(&(INITIAL_NODES_AT_E1[1].Ei), &PSI_AT_E1[1], KERNEL_AT_E1[1], curve_E1);
      // phi([2]point_St)
      eval_isog(&(INITIAL_NODES_AT_E1[1].Qi), PSI_AT_E1[1], KERNELS_AT_E1[0]);
      // phi(point_Tt)
      eval_isog(&(INITIAL_NODES_AT_E1[1].Pi), PSI_AT_E1[1], point_Tt);
		
      // ----- CASE RELATED WITH THE POINT point_PS + point_Tt
      // The Y coordinate of Px + Qx is equals zero
		memcpy(KERNEL_AT_E1[2].Y[0], zeroM, sizeof(uint64_t) * WORD_N);
		memcpy(KERNEL_AT_E1[2].Y[1], zeroM, sizeof(uint64_t) * WORD_N);
		// The Z coordinate of Px + Qx is equals Pi.Z*Qi.Z
		MUL_Fp2(KERNEL_AT_E1[2].Z, KERNEL_AT_E1[0].Z, KERNEL_AT_E1[1].Z);
		// The X coordinate of Px + Qx is equals Pi.X*(Qi.Z)^2 + (Pi.Z^2)*Qi.X
		SQR_Fp2(ZqZq, KERNEL_AT_E1[1].Z);
		MUL_Fp2(XpZqZq, KERNEL_AT_E1[0].X, ZqZq);
		SQR_Fp2(ZpZp, KERNEL_AT_E1[0].Z);
		MUL_Fp2(XqZpZp, KERNEL_AT_E1[1].X, ZpZp);
		ADD_Fp2(KERNEL_AT_E1[2].X, XpZqZq, XqZpZp);
		NEG_Fp2(KERNEL_AT_E1[2].X, KERNEL_AT_E1[2].X);
				
      // We compute the isogeny and the isogenous curve
      get_isog(&(INITIAL_NODES_AT_E1[2].Ei), &PSI_AT_E1[2], KERNEL_AT_E1[2], curve_E1);
      // phi([2]point_St)
      eval_isog(&(INITIAL_NODES_AT_E1[2].Qi), PSI_AT_E1[2], KERNELS_AT_E1[0]);
      // phi(point_St + point_Tt)
      eval_isog(&(INITIAL_NODES_AT_E1[2].Pi), PSI_AT_E1[2], KERNELS_AT_E1[2]);		
           
      // ---
      TOTAL_ADDITIONS += 1;
      TOTAL_ISOG_COMP += 3;
      TOTAL_ISOG_EVAL += 6;
      TOTAL_DOUBLINGS += (e_DIVIDED_BY_2 - 1) * 2;      
      
		// -------
		// -------
		
		omp_set_num_threads(3);
		#pragma omp parallel private(ctr_i) shared(NODES_AT_E1,Xs_AT_E1)
		{
			ctr_i = omp_get_thread_num();
			build_initial_nodes(NODES_AT_E1, Xs_AT_E1, INITIAL_NODES_AT_E1[ctr_i], ctr_i, 1, ctr_i);
		};

		// We use OpenMP for parallel finding of the collision
		printf("\"We'll start with the collision search.\";\n");
		uint8_t STOP = 0x00;
		uint64_t SOLUTION[2];		
		uint64_t acc_isog = 0;
		
		omp_set_num_threads(NUMBER_OF_CORES);
		#pragma omp parallel shared(SOLUTION,STOP) private(ctr_i) reduction(+:TOTAL_CLOCK_CYCLES_SEARCH)
		{
			ctr_i = omp_get_thread_num();
			local_start = get_cycles();
			Depth_First_Search_at_E1(&STOP, SOLUTION, NODES_AT_E1[ctr_i], Xs_AT_E1[ctr_i], LEVEL + 1, ctr_i);
			local_end = get_cycles();
			TOTAL_CLOCK_CYCLES_SEARCH += (local_end - local_start);
		}
		
		printf("CLOCK_CYCLES_SEARCH := %jd;\n", TOTAL_CLOCK_CYCLES_SEARCH);
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
		{
			TOTAL_DOUBLINGS += CURRENT_DOUBLINGS[ctr_i],
			TOTAL_ADDITIONS += CURRENT_ADDITIONS[ctr_i],
			TOTAL_ISOG_COMP += CURRENT_ISOG_COMP[ctr_i];
			TOTAL_ISOG_EVAL += CURRENT_ISOG_EVAL[ctr_i];
		}
		printf("DBLS_AT_E1 := %jd;\tADDS_AT_E1 := %jd;\tCOMP_AT_E1 := %jd;\tEVAL_AT_E1 := %jd;\n", TOTAL_DOUBLINGS, TOTAL_ADDITIONS, TOTAL_ISOG_COMP,TOTAL_ISOG_EVAL);
		
		
		printf("clock_cycles := Sprintf(\"clock_cycles_dfs.append(%%o)\", CLOCK_CYCLES_SEARCH + CLOCK_CYCLES_SORT + CLOCK_CYCLES_BUILD);\n");
		printf("running_time_dbls := Sprintf(\"running_time_dfs_dbls.append(%%o)\", DBLS_AT_E0 + DBLS_AT_E1);\n");
		printf("running_time_adds := Sprintf(\"running_time_dfs_adds.append(%%o)\", ADDS_AT_E0 + ADDS_AT_E1);\n");
		printf("running_time_comp := Sprintf(\"running_time_dfs_comp.append(%%o)\", COMP_AT_E0 + COMP_AT_E1);\n");
		printf("running_time_eval := Sprintf(\"running_time_dfs_eval.append(%%o)\", EVAL_AT_E0 + EVAL_AT_E1);\n");
		
		printf("PrintFile(\"./dfs_results.py\", clock_cycles);\n");
		printf("PrintFile(\"./dfs_results.py\", running_time_dbls);\n");					
		printf("PrintFile(\"./dfs_results.py\", running_time_adds);\n");					
		printf("PrintFile(\"./dfs_results.py\", running_time_comp);\n");					
		printf("PrintFile(\"./dfs_results.py\", running_time_eval);\n");					
		// Kernel recovery
		Point Pt_0, Qt_0, St_1, Tt_1, ZERO_POINT;
		
		Curve curve_E0g = isogenies_decomposition_at_E0(SOLUTION[0], &Pt_0, &Qt_0, Pt, Qt);
		Curve curve_E1g = isogenies_decomposition_at_E1(SOLUTION[1], &St_1, &Tt_1, point_St, point_Tt);
						
		Pt_MDSM(&ZERO_POINT, Tt_1,  e_DIVIDED_BY_2 - 1, curve_E1g.A);
				
		Isomorphism iota;
		
		if( is_Zero_point(ZERO_POINT) == 1)
			Point_Assign(&Tt_1, St_1);
				
		get_isomorphism(&iota, curve_E1g, curve_E0g);
		eval_isomorphism(&Tt_1, iota, Tt_1);
		
		Curve LAST_CURVE;
		get_isogenous_curve(&LAST_CURVE, Tt_1, Pt_MDSM, get_isog, eval_isog, curve_E0g, e_DIVIDED_BY_2);
		eval_isogeny(&Pt_0, &Qt_0, Tt_1, Pt_MDSM, get_isog, eval_isog, curve_E0g, e_DIVIDED_BY_2);
		
		Pt_MDSM(&ZERO_POINT, Qt_0,  (e_DIVIDED_BY_2*2) - 1, LAST_CURVE.A);
		
		if( is_Zero_point(ZERO_POINT) == 1)
		{
			Point_Assign(&ZERO_POINT, Qt_0);
			Point_Assign(&Qt_0, Pt_0);
			Point_Assign(&Pt_0, ZERO_POINT);
			
			print_ELE(Qt.X ,WORD_N, "K_P_X");
			print_ELE(Qt.Y ,WORD_N, "K_P_Y");
			print_ELE(Qt.Z ,WORD_N, "K_P_Z");
			
			print_ELE(Pt.X ,WORD_N, "K_Q_X");
			print_ELE(Pt.Y ,WORD_N, "K_Q_Y");
			print_ELE(Pt.Z ,WORD_N, "K_Q_Z");
			
		}
		else
		{
			print_ELE(Pt.X ,WORD_N, "K_P_X");
			print_ELE(Pt.Y ,WORD_N, "K_P_Y");
			print_ELE(Pt.Z ,WORD_N, "K_P_Z");
			
			print_ELE(Qt.X ,WORD_N, "K_Q_X");
			print_ELE(Qt.Y ,WORD_N, "K_Q_Y");
			print_ELE(Qt.Z ,WORD_N, "K_Q_Z");
		}
		
		
		printf("K_P := E![(K_P_X  * Rinv)/ ( (K_P_Z  * Rinv)^2), (K_P_Y  * Rinv) / ((K_P_Z * Rinv)^3)];\n");
		printf("K_Q := E![(K_Q_X  * Rinv)/ ( (K_Q_Z  * Rinv)^2), (K_Q_Y  * Rinv) / ((K_Q_Z * Rinv)^3)];\n");

		// ------
		uint64_t Priv_COEF[WORD_N];
		PohligHellman_alg(Priv_COEF, Pt_0, Qt_0, e_DIVIDED_BY_2 * 2, LAST_CURVE.A);
		print_NUM(Priv_COEF, 2, 0, "Priv_COEF");
		printf("K := K_P - Priv_COEF * K_Q;\n");
		
		printf("C_isog := d_pow_e_isogenous_curve(K, E, %jd, %jd);\n", DEGREE, e_DIVIDED_BY_2 * 2);
		printf("jInvariant(C_isog) eq jInvariant(E_isog);\n");
		
		printf("exit;\n");
		// ---
		for(ctr_i = 0; ctr_i < NUMBER_OF_CORES; ctr_i++)
			free(TREE_LEAVES_ROOTED_AT_E0[ctr_i]);
		
		free(TREE_LEAVES_ROOTED_AT_E0);
		free(CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0);
		free(CURRENT_DOUBLINGS);
		free(CURRENT_ADDITIONS);
		free(CURRENT_ISOG_COMP);
		free(CURRENT_ISOG_EVAL);
		return 0;
	}
	else{
		printf("Invalid case: It should take 2 args.\n\t\t 1) i : 3*2^i cores to be used, and\n\t\t 2) A (Alice) or B (Bob).\n" );
		return 1; 
	}
	
};
