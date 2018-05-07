#include "dfs.h"

// ---------------------------------------------------------------
// swap() function :::::::::::::::::::::::::::::::::::::::::
// A utility function to swap two elements
void swap(T_naive *a, T_naive *b)
{
   T_naive t;
   memcpy(t.j_inv, a->j_inv, sizeof(uint64_t) * WORD_N * 2);
   t.point = a->point;

   memcpy(a->j_inv, b->j_inv, sizeof(uint64_t) * WORD_N * 2);
   a->point = b->point;
	
   memcpy(b->j_inv, t.j_inv, sizeof(uint64_t) * WORD_N * 2);
   b->point = t.point;
};


// ---------------------------------------------------------------
// partition() function :::::::::::::::::::::::::::::::::::::::::
// This function takes last element as pivot, places the pivot 
// element at its correct position in sorted array, and places 
// all smaller (smaller than pivot) to left of pivot and all 
// greater elements to right of pivot 
int partition (T_naive arr[], int low, int high)
{
   uint64_t pivot[2][WORD_N];
   memcpy(pivot, arr[high].j_inv, sizeof(uint64_t) * WORD_N * 2);    // pivot
   int j, i = (low - 1);  // Index of smaller element

   for(j = low; j <= high- 1; j++)
   {
      // If current element is smaller than or
      // equal to pivot
      if(compare_in_Fp2(arr[j].j_inv, pivot) <= 0)
      {
	 i++;    // increment index of smaller element
	 swap(&arr[i], &arr[j]);
      }
   }
   swap(&arr[i + 1], &arr[high]);
   return (i + 1);
};


// ---------------------------------------------------------------
// quickSort() function :::::::::::::::::::::::::::::::::::::::::
// inputs:
//         arr[] : Array to be sorted,
//           low : Starting index,
//          high : Ending index 
// output:
//        the sorted list arr[]
void quickSort(T_naive arr[], int low, int high)
{
   if (low < high)
   {
      /* pi is partitioning index, arr[p] is now
       * at right place */
      int pi = partition(arr, low, high);

      // Separately sort elements before
      // partition and after partition
      quickSort(arr, low, pi - 1);
      quickSort(arr, pi + 1, high);
   }
};

// ---------------------------------------------------------------
// binarySearch() function :::::::::::::::::::::::::::::::::::::::
// inputs:
//         local_seq [] : a sorted list of elements of type T_naive
//           local_jinv : the j-invariant to be located in the list
// output:
//         the position in the list which corresponds to the j-invariant, or
//         NULL if the j-invariant is not in the list 
int64_t binarySearch(T_naive local_seq[], uint64_t local_jinv[2][WORD_N], uint64_t low, uint64_t high)
{
   int mid;
   while((int)high >= (int)low)
   {
      mid = (low + high)/2;
      if( compare_in_Fp2(local_seq[mid].j_inv, local_jinv) == 0)
	 return mid;
      else if( compare_in_Fp2(local_seq[mid].j_inv, local_jinv) == -1)
	 low = mid + 1;
      else
	 high = mid - 1;
   }

   return -1;
};

// ---------------------------------------------------------------
// isogenies_decomposition_at_E0() function :::::::::::::::::::::::::::
//     This function construct the Leaves of the tree rooted at E0
// inputs:
//            X : the trajectory of a given leave of the tree rooted at E0,
//      input_0 : a point to be evaluated by the degree-2^e/2 isogeny which
//                corresponds to the trajectory X
//      input_1 : a point to be evaluated by the degree-2^e/2 isogeny which
//                corresponds to the trajectory X
// output:
//      output_0 : the evaluation of input_0 under the degree-2^e/2 isogeny
//                 which correspond to the trajectory X
//      output_1 : the evaluation of input_1 under the degree-2^e/2 isogeny
//                 which correspond to the trajectory X
Curve isogenies_decomposition_at_E0(uint64_t X, Point *output_0, Point *output_1, Point input_0, Point input_1)
{
   T_dfs INITIAL_NODE;
   Point kernel, Pi_plus_Qi, TMP_Pi, TMP_Qi;
   uint64_t itr_k, mask;

   mask = 0x0000000000000003 << (e_DIVIDED_BY_2 - 1);
   mask = (mask & X) >> (e_DIVIDED_BY_2 - 1);
   
   // Copying the first order-2^(e/2 - 1) points in the 2-isogenous curve
   Point_Assign(&(INITIAL_NODE.Pi), INITIAL_NODES_AT_E0[mask].Pi);
   Point_Assign(&(INITIAL_NODE.Qi), INITIAL_NODES_AT_E0[mask].Qi);
   // Copying the first 2-isogenous elliptic curve
   memcpy(INITIAL_NODE.Ei.A, INITIAL_NODES_AT_E0[mask].Ei.A, sizeof(uint64_t) * WORD_N * 2 );
   memcpy(INITIAL_NODE.Ei.B, INITIAL_NODES_AT_E0[mask].Ei.B, sizeof(uint64_t) * WORD_N * 2 );
   
   eval_isog(output_0, PHI_AT_E0[mask], input_0);
   eval_isog(output_1, PHI_AT_E0[mask], input_1);
   
   mask = 0x0000000000000001 << (e_DIVIDED_BY_2 - 1);
   for(itr_k = 1; itr_k < (e_DIVIDED_BY_2 - 1); itr_k++)
   {
      // Pi + Qi
      Pt_ADD(&Pi_plus_Qi, INITIAL_NODE.Pi, INITIAL_NODE.Qi, INITIAL_NODE.Ei.A);
      // [2](Pi + Qi)
      Pt_DSM(&TMP_Qi, Pi_plus_Qi, INITIAL_NODE.Ei.A);
      // [2]Pi
      Pt_DSM(&TMP_Pi, INITIAL_NODE.Pi, INITIAL_NODE.Ei.A);
      
      mask = mask >> 1;
      if( (mask & X) == 0x0000000000000000)
      {
	 // ----- CASE RELATED WITH THE POINT Pi
	 // We compute the kernel of phi : <[2^(e/2 - 1 - level)]Pi>
	 Pt_MDSM(&kernel, TMP_Pi, e_DIVIDED_BY_2 - 2 - itr_k, INITIAL_NODE.Ei.A);
	 // We compute the isogeny and the isogenous curve
	 get_isog(&(INITIAL_NODE.Ei), &PHI_LOCAL, kernel, INITIAL_NODE.Ei);
	 // phi(Pi)
	 eval_isog(&(INITIAL_NODE.Pi), PHI_LOCAL, INITIAL_NODE.Pi);
	 // phi([2](Pi + Qi))
	 eval_isog(&(INITIAL_NODE.Qi), PHI_LOCAL, TMP_Qi);
      }
      else
      {
	 // ----- CASE RELATED WITH THE POINT Pi + Qi
	 Pt_MDSM(&kernel, TMP_Qi, e_DIVIDED_BY_2 - 2 - itr_k, INITIAL_NODE.Ei.A);
	 // We compute the isogeny and the isogenous curve
	 get_isog(&(INITIAL_NODE.Ei), &PHI_LOCAL, kernel, INITIAL_NODE.Ei);
	 // phi([2]Pi)
	 eval_isog(&(INITIAL_NODE.Qi), PHI_LOCAL, TMP_Pi);
	 // phi(Pi + Qi)
	 eval_isog(&(INITIAL_NODE.Pi), PHI_LOCAL, Pi_plus_Qi);
      }
      
      eval_isog(output_0, PHI_LOCAL, *output_0);
      eval_isog(output_1, PHI_LOCAL, *output_1);
   }
   
   // The Y coordinate of Px + Qx is equals zero
   memcpy(Pi_plus_Qi.Y[0], zeroM, sizeof(uint64_t) * WORD_N);
   memcpy(Pi_plus_Qi.Y[1], zeroM, sizeof(uint64_t) * WORD_N);
   // The Z coordinate of Px + Qx is equals Pi.Z*Qi.Z
   MUL_Fp2(Pi_plus_Qi.Z, INITIAL_NODE.Pi.Z, INITIAL_NODE.Qi.Z);
   // The X coordinate of Px + Qx is equals Pi.X*(Qi.Z)^2 + (Pi.Z^2)*Qi.X
   uint64_t XpZqZq[2][WORD_N], XqZpZp[2][WORD_N], ZpZp[2][WORD_N], ZqZq[2][WORD_N];
   SQR_Fp2(ZqZq, INITIAL_NODE.Qi.Z);
   MUL_Fp2(XpZqZq, INITIAL_NODE.Pi.X, ZqZq);
   SQR_Fp2(ZpZp, INITIAL_NODE.Pi.Z);
   MUL_Fp2(XqZpZp, INITIAL_NODE.Qi.X, ZpZp);
   ADD_Fp2(Pi_plus_Qi.X, XpZqZq, XqZpZp);
   NEG_Fp2(Pi_plus_Qi.X, Pi_plus_Qi.X);
   
   mask = mask >> 1;
   if( (mask & X) == 0x0000000000000000)
   {
      // ----- CASE RELATED WITH THE POINT Pi
      // We compute the isogeny and the isogenous curve
      get_isog(&(INITIAL_NODE.Ei), &PHI_LOCAL, INITIAL_NODE.Pi, INITIAL_NODE.Ei);
   }
   else
   {
      // ----- CASE RELATED WITH THE POINT Pi + Qi
      // We compute the isogeny and the isogenous curve
      get_isog(&(INITIAL_NODE.Ei), &PHI_LOCAL, Pi_plus_Qi, INITIAL_NODE.Ei);
   }
   
   eval_isog(output_0, PHI_LOCAL, *output_0);
   eval_isog(output_1, PHI_LOCAL, *output_1);
   
   return INITIAL_NODE.Ei;
};

// ---------------------------------------------------------------
// isogenies_decomposition_at_E1() function :::::::::::::::::::::::::::
//     This function construct the Leaves of the tree rooted at E1
// inputs:
//            X : the trajectory of a given leave of the tree rooted at E1,
//      input_0 : a point to be evaluated by the degree-2^e/2 isogeny which
//                corresponds to the trajectory X
//      input_1 : a point to be evaluated by the degree-2^e/2 isogeny which
//                corresponds to the trajectory X
// output:
//      output_0 : the evaluation of input_0 under the degree-2^e/2 isogeny
//                 which correspond to the trajectory X
//      output_1 : the evaluation of input_1 under the degree-2^e/2 isogeny
//                 which correspond to the trajectory X
Curve isogenies_decomposition_at_E1(uint64_t X, Point *output_0, Point *output_1, Point input_0, Point input_1)
{
   T_dfs INITIAL_NODE;
   Point kernel, Pi_plus_Qi, TMP_Pi, TMP_Qi;
   uint64_t itr_k, mask;

   mask = 0x0000000000000003 << (e_DIVIDED_BY_2 - 1);
   mask = (mask & X) >> (e_DIVIDED_BY_2 - 1);

   // Copying the first order-2^(e/2 - 1) points in the 2-isogenous curve
   Point_Assign(&(INITIAL_NODE.Pi), INITIAL_NODES_AT_E1[mask].Pi);
   Point_Assign(&(INITIAL_NODE.Qi), INITIAL_NODES_AT_E1[mask].Qi);
   // Copying the first 2-isogenous elliptic curve
   memcpy(INITIAL_NODE.Ei.A, INITIAL_NODES_AT_E1[mask].Ei.A, sizeof(uint64_t) * WORD_N * 2 );
   memcpy(INITIAL_NODE.Ei.B, INITIAL_NODES_AT_E1[mask].Ei.B, sizeof(uint64_t) * WORD_N * 2 );

   eval_isog(output_0, PSI_AT_E1[mask], input_0);
   eval_isog(output_1, PSI_AT_E1[mask], input_1);
   
   mask = 0x0000000000000001 << (e_DIVIDED_BY_2 - 1);

   for(itr_k = 1; itr_k < (e_DIVIDED_BY_2 - 1); itr_k++)
   {
      // Pi + Qi
      Pt_ADD(&Pi_plus_Qi, INITIAL_NODE.Pi, INITIAL_NODE.Qi, INITIAL_NODE.Ei.A);
      // [2](Pi + Qi)
      Pt_DSM(&TMP_Qi, Pi_plus_Qi, INITIAL_NODE.Ei.A);
      // [2]Pi
      Pt_DSM(&TMP_Pi, INITIAL_NODE.Pi, INITIAL_NODE.Ei.A);
      
      mask = mask >> 1;
      if( (mask & X) == 0x0000000000000000)
      {
	 // ----- CASE RELATED WITH THE POINT Pi
	 // We compute the kernel of phi : <[2^(e/2 - 1 - level)]Pi>
	 Pt_MDSM(&kernel, TMP_Pi, e_DIVIDED_BY_2 - 2 - itr_k, INITIAL_NODE.Ei.A);
	 // We compute the isogeny and the isogenous curve
	 get_isog(&(INITIAL_NODE.Ei), &PSI_LOCAL, kernel, INITIAL_NODE.Ei);
	 // phi(Pi)
	 eval_isog(&(INITIAL_NODE.Pi), PSI_LOCAL, INITIAL_NODE.Pi);
	 // phi([2](Pi + Qi))
	 eval_isog(&(INITIAL_NODE.Qi), PSI_LOCAL, TMP_Qi);
      }
      else
      {
	 // ----- CASE RELATED WITH THE POINT Pi + Qi
	 Pt_MDSM(&kernel, TMP_Qi, e_DIVIDED_BY_2 - 2 - itr_k, INITIAL_NODE.Ei.A);
	 // We compute the isogeny and the isogenous curve
	 get_isog(&(INITIAL_NODE.Ei), &PSI_LOCAL, kernel, INITIAL_NODE.Ei);
	 // phi([2]Pi)
	 eval_isog(&(INITIAL_NODE.Qi), PSI_LOCAL, TMP_Pi);
	 // phi(Pi + Qi)
	 eval_isog(&(INITIAL_NODE.Pi), PSI_LOCAL, Pi_plus_Qi);
      }
      
      eval_isog(output_0, PSI_LOCAL, *output_0);
      eval_isog(output_1, PSI_LOCAL, *output_1);
   }
   
   // The Y coordinate of Px + Qx is equals zero
   memcpy(Pi_plus_Qi.Y[0], zeroM, sizeof(uint64_t) * WORD_N);
   memcpy(Pi_plus_Qi.Y[1], zeroM, sizeof(uint64_t) * WORD_N);
   // The Z coordinate of Px + Qx is equals Pi.Z*Qi.Z
   MUL_Fp2(Pi_plus_Qi.Z, INITIAL_NODE.Pi.Z, INITIAL_NODE.Qi.Z);
   // The X coordinate of Px + Qx is equals Pi.X*(Qi.Z)^2 + (Pi.Z^2)*Qi.X
   uint64_t XpZqZq[2][WORD_N], XqZpZp[2][WORD_N], ZpZp[2][WORD_N], ZqZq[2][WORD_N];
   SQR_Fp2(ZqZq, INITIAL_NODE.Qi.Z);
   MUL_Fp2(XpZqZq, INITIAL_NODE.Pi.X, ZqZq);
   SQR_Fp2(ZpZp, INITIAL_NODE.Pi.Z);
   MUL_Fp2(XqZpZp, INITIAL_NODE.Qi.X, ZpZp);
   ADD_Fp2(Pi_plus_Qi.X, XpZqZq, XqZpZp);
   NEG_Fp2(Pi_plus_Qi.X, Pi_plus_Qi.X);
   
   mask = mask >> 1;
   if( (mask & X) == 0x0000000000000000)
   {
      // ----- CASE RELATED WITH THE POINT Pi
      // We compute the isogeny and the isogenous curve
      get_isog(&(INITIAL_NODE.Ei), &PSI_LOCAL, INITIAL_NODE.Pi, INITIAL_NODE.Ei);
   }
   else
   {
      // ----- CASE RELATED WITH THE POINT Pi + Qi
      // We compute the isogeny and the isogenous curve
      get_isog(&(INITIAL_NODE.Ei), &PSI_LOCAL, Pi_plus_Qi, INITIAL_NODE.Ei);
   }
   
   eval_isog(output_0, PSI_LOCAL, *output_0);
   eval_isog(output_1, PSI_LOCAL, *output_1);
   
   return INITIAL_NODE.Ei;
};

// ---------------------------------------------------------------
// build_initial_nodes() function :::::::::::::::::::::::::::
//     This function construct the Leaves of the tree rooted at E
// inputs:
//         NODE : the initial node in the tree for which we will 
//                compute all its leaves,
//            X : The intial seed correspoding with NODE, i.e.,
//                the trajectory which compute NODE,
//        level : the current level of the tree,
//     kth_proc : the process id, i.e., the number of the core to 
//                be used.
// output:
//      NODES[] : All the nodes at the level "level",
//         Xs[] : The seed for each node at the level "level", i.e.,
//                the corresponding trajectories
void build_initial_nodes(T_dfs NODES[], uint64_t Xs[], T_dfs NODE, uint64_t X, uint64_t level, uint64_t kth_proc)
{
   if(level == (LEVEL + 1))
   {
      Point_Assign(&(NODES[X].Pi), NODE.Pi);
      Point_Assign(&(NODES[X].Qi), NODE.Qi);
      memcpy(NODES[X].Ei.A, NODE.Ei.A, sizeof(uint64_t) * WORD_N * 2);
      memcpy(NODES[X].Ei.B, NODE.Ei.B, sizeof(uint64_t) * WORD_N * 2);
      Xs[X] = X;
   }
   else if(level < (LEVEL + 1))
   {
      T_dfs child;
      Isogeny phi;
      Point kernel, tmp, Pi_plus_Qi, TMP_Pi, TMP_Qi;
      // Pi + Qi
      Pt_ADD(&Pi_plus_Qi, NODE.Pi, NODE.Qi, NODE.Ei.A);
      // [2](Pi + Qi)
      Pt_DSM(&TMP_Qi, Pi_plus_Qi, NODE.Ei.A);
      // [2]Pi
      Pt_DSM(&TMP_Pi, NODE.Pi, NODE.Ei.A);
      
      // ----- CASE RELATED WITH THE POINT Pi
      // We compute the kernel of phi : <[2^(e/2 - 1 - level)]Pi>
      Pt_MDSM(&kernel, TMP_Pi, e_DIVIDED_BY_2 - 2 - level, NODE.Ei.A);      
      // We compute the isogeny and the isogenous curve
      get_isog(&(child.Ei), &phi, kernel, NODE.Ei);
      // phi(Pi)
      eval_isog(&(child.Pi), phi, NODE.Pi);
      // phi([2](Pi + Qi))
      eval_isog(&(child.Qi), phi, TMP_Qi);
      // We apply recursion for computing the leaves
      build_initial_nodes(NODES, Xs, child, X << 1, level + 1, kth_proc);		// SIDE OF Pi
      
      // ----- CASE RELATED WITH THE POINT Pi + Qi
      Pt_MDSM(&kernel, TMP_Qi, e_DIVIDED_BY_2 - 2 - level, NODE.Ei.A);      
      // We compute the isogeny and the isogenous curve
      get_isog(&(child.Ei), &phi, kernel, NODE.Ei);
      // phi([2]Pi)
      eval_isog(&(child.Qi), phi, TMP_Pi);
      // phi(Pi + Qi)
      eval_isog(&(child.Pi), phi, Pi_plus_Qi);
      // We apply recursion for computing the leaves
      build_initial_nodes(NODES, Xs, child, (X << 1) ^ 0x0000000000000001, level + 1, kth_proc);		// SIDE OF Qi
      
      // ---
      CURRENT_ADDITIONS[kth_proc] += 1;
      CURRENT_ISOG_COMP[kth_proc] += 2;
      CURRENT_ISOG_EVAL[kth_proc] += 4;
      CURRENT_DOUBLINGS[kth_proc] += (e_DIVIDED_BY_2 - 2 - level) * 2;
   }
   return ;
}

// ---------------------------------------------------------------
// Depth_First_Search_at_E0() function :::::::::::::::::::::::::::
//     This function constructs the Leaves of the tree rooted at E0
// inputs:
//         NODE : the initial node in the tree for which we will 
//                compute all its leaves,
//            X : The intial seed correspoding with NODE, i.e.,
//                the trajectory which compute NODE,
//        level : the current level of the tree,
//     kth_proc : the process id, i.e., the number of the core to 
//                be used.
// output:
//         Clock Cycles required for the construction of all the 
//         Leaves of the tree rooted at E0
void Depth_First_Search_at_E0(T_dfs NODE, uint64_t X, uint64_t level, uint64_t kth_proc)
{
   if(level == (e_DIVIDED_BY_2 - 1))
   {
      Curve Ei;
      Isogeny phi;
      Point Pi_plus_Qi;
      
      // ----- CASE RELATED WITH THE POINT Pi
      CURRENT_ISOG_COMP[kth_proc] += 1;
      
      get_isog(&Ei, &phi, NODE.Pi, NODE.Ei);
      jInvariant((&TREE_LEAVES_ROOTED_AT_E0[kth_proc][CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[kth_proc]])->j_inv, Ei);
      (&TREE_LEAVES_ROOTED_AT_E0[kth_proc][CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[kth_proc]])->point = X << 1;
      CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[kth_proc] += 1;
      
      // ----- CASE RELATED WITH THE POINT Pi + Qi
      CURRENT_ISOG_COMP[kth_proc] += 1;
      
      // The Y coordinate of Px + Qx is equals zero
      memcpy(Pi_plus_Qi.Y[0], zeroM, sizeof(uint64_t) * WORD_N);
      memcpy(Pi_plus_Qi.Y[1], zeroM, sizeof(uint64_t) * WORD_N);
      // The Z coordinate of Px + Qx is equals Pi.Z*Qi.Z
      MUL_Fp2(Pi_plus_Qi.Z, NODE.Pi.Z, NODE.Qi.Z);
      // The X coordinate of Px + Qx is equals Pi.X*(Qi.Z)^2 + (Pi.Z^2)*Qi.X
      uint64_t XpZqZq[2][WORD_N], XqZpZp[2][WORD_N], ZpZp[2][WORD_N], ZqZq[2][WORD_N];
      SQR_Fp2(ZqZq, NODE.Qi.Z);
      MUL_Fp2(XpZqZq, NODE.Pi.X, ZqZq);
      SQR_Fp2(ZpZp, NODE.Pi.Z);
      MUL_Fp2(XqZpZp, NODE.Qi.X, ZpZp);
      ADD_Fp2(Pi_plus_Qi.X, XpZqZq, XqZpZp);
      NEG_Fp2(Pi_plus_Qi.X, Pi_plus_Qi.X);
      
      get_isog(&Ei, &phi, Pi_plus_Qi, NODE.Ei);      
      jInvariant((&TREE_LEAVES_ROOTED_AT_E0[kth_proc][CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[kth_proc]])->j_inv, Ei);
      (&TREE_LEAVES_ROOTED_AT_E0[kth_proc][CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[kth_proc]])->point = (X << 1) ^ 0x0000000000000001;
      CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[kth_proc] += 1;
   }
   else if(level < (e_DIVIDED_BY_2 - 1))
   {
      T_dfs child;
      Isogeny phi;
      Point kernel, tmp, Pi_plus_Qi, TMP_Pi, TMP_Qi;
      // Pi + Qi
      Pt_ADD(&Pi_plus_Qi, NODE.Pi, NODE.Qi, NODE.Ei.A);
      // [2](Pi + Qi)
      Pt_DSM(&TMP_Qi, Pi_plus_Qi, NODE.Ei.A);
      // [2]Pi
      Pt_DSM(&TMP_Pi, NODE.Pi, NODE.Ei.A);
      
      // ----- CASE RELATED WITH THE POINT Pi
      // We compute the kernel of phi : <[2^(e/2 - 1 - level)]Pi>
      Pt_MDSM(&kernel, TMP_Pi, e_DIVIDED_BY_2 - 2 - level, NODE.Ei.A);
      // We compute the isogeny and the isogenous curve
      get_isog(&(child.Ei), &phi, kernel, NODE.Ei);
      // phi(Pi)
      eval_isog(&(child.Pi), phi, NODE.Pi);
      // phi([2](Pi + Qi))
      eval_isog(&(child.Qi), phi, TMP_Qi);
      // We apply recursion for computing the leaves
      Depth_First_Search_at_E0(child, X << 1, level + 1, kth_proc);				// SIDE OF Pi
      
      // ----- CASE RELATED WITH THE POINT Pi + Qi
      Pt_MDSM(&kernel, TMP_Qi, e_DIVIDED_BY_2 - 2 - level, NODE.Ei.A);
      // We compute the isogeny and the isogenous curve
      get_isog(&(child.Ei), &phi, kernel, NODE.Ei);
      // phi([2]Pi)
      eval_isog(&(child.Qi), phi, TMP_Pi);
      // phi(Pi + Qi)
      eval_isog(&(child.Pi), phi, Pi_plus_Qi);
      // We apply recursion for computing the leaves
      Depth_First_Search_at_E0(child,(X << 1) ^ 0x0000000000000001, level + 1, kth_proc);	// SIDE OF Qi
      
      // ---
      CURRENT_ADDITIONS[kth_proc] += 1;
      CURRENT_ISOG_COMP[kth_proc] += 2;
      CURRENT_ISOG_EVAL[kth_proc] += 4;
      CURRENT_DOUBLINGS[kth_proc] += (e_DIVIDED_BY_2 - 1 - level) * 2;
   }
   return ;
}

// ---------------------------------------------------------------
// Depth_First_Search_at_E1() function :::::::::::::::::::::::::::
//     This function computes the Leaves of the tree rooted at E1
// inputs:
//         BOOL : this variable determines if the collision has
//                been found,
//         NODE : the initial node in the tree for which we will 
//                compute all its leaves,
//            X : The intial seed correspoding with NODE, i.e.,
//                the trajectory which compute NODE,
//        level : the current level of the tree,
//     kth_proc : the process id, i.e., the number of the core to 
//                be used.
// output:
//         Clock Cycles required for the computation of all the 
//         Leaves of the tree rooted at E1,
//         COLLISION : a collision, i.e., an element of 
//                     [{0,1}x{0,1,2}x{0,1}^(e/2 - 1)] x [{0,1}x{0,1,2}x{0,1}^(e/2 - 1)]
void Depth_First_Search_at_E1(uint8_t *FLAG, uint64_t COLLISION[2], T_dfs NODE, uint64_t X, uint64_t level, uint64_t kth_proc)
{
   if(*FLAG != 0x01)
   {
      if(level == (e_DIVIDED_BY_2 - 1))
      {
	 uint64_t JINV[2][WORD_N], local_i;
	 int64_t local_collision;
	 Curve Ei;
	 Isogeny phi;
	 Point Pi_plus_Qi;
	 
	 // ----- CASE RELATED WITH THE POINT Pi
	 CURRENT_ISOG_COMP[kth_proc] += 1;
	 
	 get_isog(&Ei, &phi, NODE.Pi, NODE.Ei);
	 jInvariant(JINV, Ei);
	 for(local_i = 0; local_i < NUMBER_OF_CORES; local_i++)
	 {
	    local_collision = binarySearch(TREE_LEAVES_ROOTED_AT_E0[local_i], JINV, 0, CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[local_i]);
	    if( local_collision != -1 )
	    {
	       COLLISION[0] = TREE_LEAVES_ROOTED_AT_E0[local_i][local_collision].point;
	       COLLISION[1] = X << 1;			 
	       *FLAG = 0x01;
	       return;
	    }
	 }
	 
	 // ----- CASE RELATED WITH THE POINT Pi + Qi
	 CURRENT_ISOG_COMP[kth_proc] += 1;
	 // The Y coordinate of Px + Qx is equals zero
	 memcpy(Pi_plus_Qi.Y[0], zeroM, sizeof(uint64_t) * WORD_N);
	 memcpy(Pi_plus_Qi.Y[1], zeroM, sizeof(uint64_t) * WORD_N);
	 // The Z coordinate of Px + Qx is equals Pi.Z*Qi.Z
	 MUL_Fp2(Pi_plus_Qi.Z, NODE.Pi.Z, NODE.Qi.Z);
	 // The X coordinate of Px + Qx is equals Pi.X*(Qi.Z)^2 + (Pi.Z^2)*Qi.X
	 uint64_t XpZqZq[2][WORD_N], XqZpZp[2][WORD_N], ZpZp[2][WORD_N], ZqZq[2][WORD_N];
	 SQR_Fp2(ZqZq, NODE.Qi.Z);
	 MUL_Fp2(XpZqZq, NODE.Pi.X, ZqZq);
	 SQR_Fp2(ZpZp, NODE.Pi.Z);
	 MUL_Fp2(XqZpZp, NODE.Qi.X, ZpZp);
	 ADD_Fp2(Pi_plus_Qi.X, XpZqZq, XqZpZp);
	 NEG_Fp2(Pi_plus_Qi.X, Pi_plus_Qi.X);
	 
	 get_isog(&Ei, &phi, Pi_plus_Qi, NODE.Ei);      
	 jInvariant(JINV, Ei);
	 for(local_i = 0; local_i < NUMBER_OF_CORES; local_i++)
	 {
	    local_collision = binarySearch(TREE_LEAVES_ROOTED_AT_E0[local_i], JINV, 0, CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[local_i]);
	    if( local_collision != -1 )
	    {
	       COLLISION[0] = TREE_LEAVES_ROOTED_AT_E0[local_i][local_collision].point;
	       COLLISION[1] = (X << 1) ^ 0x0000000000000001;			 
	       *FLAG = 0x01;
	       return;
	    }
	 }
	 
      }
      else if(level < (e_DIVIDED_BY_2 - 1))
      {
	 T_dfs child;
	 Isogeny phi;
	 Point kernel, tmp, Pi_plus_Qi, TMP_Pi, TMP_Qi;
	 // Pi + Qi
	 Pt_ADD(&Pi_plus_Qi, NODE.Pi, NODE.Qi, NODE.Ei.A);
	 // [2](Pi + Qi)
	 Pt_DSM(&TMP_Qi, Pi_plus_Qi, NODE.Ei.A);
	 // [2]Pi
	 Pt_DSM(&TMP_Pi, NODE.Pi, NODE.Ei.A);
	 
	 // ----- CASE RELATED WITH THE POINT Pi
	 // We compute the kernel of phi : <[2^(e/2 - 1 - level)]Pi>
	 Pt_MDSM(&kernel, TMP_Pi, e_DIVIDED_BY_2 - 2 - level, NODE.Ei.A);
	 // We compute the isogeny and the isogenous curve
	 get_isog(&(child.Ei), &phi, kernel, NODE.Ei);
	 // phi(Pi)
	 eval_isog(&(child.Pi), phi, NODE.Pi);
	 // phi([2](Pi + Qi))
	 eval_isog(&(child.Qi), phi, TMP_Qi);
	 // We apply recursion for computing the leaves
	 Depth_First_Search_at_E1(FLAG, COLLISION, child, X << 1, level + 1, kth_proc);				// SIDE OF Pi
	 
	 // ----- CASE RELATED WITH THE POINT Pi + Qi
	 Pt_MDSM(&kernel, TMP_Qi, e_DIVIDED_BY_2 - 2 - level, NODE.Ei.A);
	 // We compute the isogeny and the isogenous curve
	 get_isog(&(child.Ei), &phi, kernel, NODE.Ei);
	 // phi([2]Pi)
	 eval_isog(&(child.Qi), phi, TMP_Pi);
	 // phi(Qi + Pi)
	 eval_isog(&(child.Pi), phi, Pi_plus_Qi);
	 // We apply recursion for computing the leaves
	 Depth_First_Search_at_E1(FLAG, COLLISION, child, (X << 1) ^ 0x0000000000000001, level + 1, kth_proc);	// SIDE OF Qi
	 
	 // ---
	 CURRENT_ADDITIONS[kth_proc] += 1;
	 CURRENT_ISOG_COMP[kth_proc] += 2;
	 CURRENT_ISOG_EVAL[kth_proc] += 4;
	 CURRENT_DOUBLINGS[kth_proc] += (e_DIVIDED_BY_2 - 1 - level) * 2;
      }
   }
   return ;
}
