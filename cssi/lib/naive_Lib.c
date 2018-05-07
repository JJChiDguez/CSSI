#include "naive.h"

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
// Leaves_of_the_tree_rooted_at_E0() function ::::::::::::::::::::
//     This function constructs the Leaves of the tree rooted at E0
// inputs:
//         kth_proc : the process id, i.e., the number of the core 
//                    to be used.
// output:
//         Clock Cycles required for the construction of all the 
//         Leaves of the tree rooted at E0
uint64_t Leaves_of_the_tree_rooted_at_E0(uint64_t kth_proc)
{
   clock_t local_start, local_end;
   local_start = get_cycles();
   
   uint64_t var_j, local_i;
   uint64_t MAXIMUM_NUMBER_OF_ISOGENY_TO_BE_COMPUTED = VAR_N_DIVIDED_BY_DEGREE_PLUS_1 / NUMBER_OF_CORES;
   // initial scalar for the kth_proc-th processor 
   uint64_t scalar[WORD_N];
   memcpy(scalar, zeroM, sizeof(uint64_t) * WORD_N);
   scalar[0] = MAXIMUM_NUMBER_OF_ISOGENY_TO_BE_COMPUTED * kth_proc;
   
   if(kth_proc == (NUMBER_OF_CORES - 1))
      MAXIMUM_NUMBER_OF_ISOGENY_TO_BE_COMPUTED += VAR_N_DIVIDED_BY_DEGREE_PLUS_1 % NUMBER_OF_CORES;
   
   // Three possible kernels
   Point point_0, point_1, point_2;
   // point_Pt + [scalar]point_Qt
   DBLADD(&point_0, point_Qt, scalar, curve_E0.A);
   Pt_ADD(&point_0, point_0, point_Pt, curve_E0.A);
   // point_Pt + [scalar + DEGREE^(e_DIVIDED_BY_2 - 1)]point_Qt
   Pt_ADD(&point_1, point_0, order_d_point_Qt, curve_E0.A);
   // [2 * scalar]point_Pt + point_Qt
   DBLADD(&point_2, point_Pt, scalar, curve_E0.A);
   Pt_DSM(&point_2, point_2, curve_E0.A);
   Pt_ADD(&point_2, point_2, point_Qt, curve_E0.A);
   // Three possible isogenous curves
   Curve curve_0, curve_1, curve_2;
   
   uint64_t ACC = 0;
   uint64_t COEF_K = scalar[0];
   
   for(var_j = 0; var_j < MAXIMUM_NUMBER_OF_ISOGENY_TO_BE_COMPUTED; var_j++)
   {
      get_isogenous_curve(&curve_0, point_0, Pt_MDSM, get_isog, eval_isog, curve_E0, e_DIVIDED_BY_2);
      jInvariant((&TREE_LEAVES_ROOTED_AT_E0[kth_proc][ACC])->j_inv, curve_0);
      (&TREE_LEAVES_ROOTED_AT_E0[kth_proc][ACC])->point = (COEF_K << 2) ^ 0x00;
      Pt_ADD(&point_0, point_0, point_Qt, curve_E0.A);
      ACC += 1;
      
      get_isogenous_curve(&curve_1, point_1, Pt_MDSM, get_isog, eval_isog, curve_E0, e_DIVIDED_BY_2);
      jInvariant((&TREE_LEAVES_ROOTED_AT_E0[kth_proc][ACC])->j_inv, curve_1);
      (&TREE_LEAVES_ROOTED_AT_E0[kth_proc][ACC])->point = (COEF_K << 2) ^ 0x01;
      Pt_ADD(&point_1, point_1, point_Qt, curve_E0.A);
      ACC += 1;
      
      get_isogenous_curve(&curve_2, point_2, Pt_MDSM, get_isog, eval_isog, curve_E0, e_DIVIDED_BY_2);
      jInvariant((&TREE_LEAVES_ROOTED_AT_E0[kth_proc][ACC])->j_inv, curve_2);
      (&TREE_LEAVES_ROOTED_AT_E0[kth_proc][ACC])->point = (COEF_K << 2) ^ 0x02;
      Pt_ADD(&point_2, point_2, d_times_Pt, curve_E0.A);
      ACC += 1;
      
      COEF_K += 1;
   }
   local_end = get_cycles();
   
   RUNNING_TIME[kth_proc] = ACC;
   CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[kth_proc] = ACC;
      
   return local_end - local_start;
};

// ---------------------------------------------------------------
// Leaves_of_the_tree_rooted_at_E1() function ::::::::::::::::::::
//      This function computes the Leaves of the tree rooted at E1
// inputs:
//             BOOL : this variable determines if the collision has
//                    been found,
//         kth_proc : the process id, i.e., the number of the core to be used
// output:
//         Clock Cycles required for the computation of all the 
//         Leaves of the tree rooted at E1,
//         COLLISION : a collision, i.e., an element of 
//                     [{0,1}x{0,1,2}x{0,1}^(e/2 - 1)] x [{0,1}x{0,1,2}x{0,1}^(e/2 - 1)]
uint64_t Leaves_of_the_tree_rooted_at_E1(uint64_t COLLISION[2][WORD_N], uint8_t *BOOL, uint64_t kth_proc)
{
   clock_t local_start, local_end;
   local_start = get_cycles();
   
   uint64_t var_j, local_i;
   uint64_t MAXIMUM_NUMBER_OF_ISOGENY_TO_BE_COMPUTED = VAR_N_DIVIDED_BY_DEGREE_PLUS_1 / NUMBER_OF_CORES;
   // initial scalar for the kth_proc-th processor 
   uint64_t scalar[WORD_N];
   memcpy(scalar, zeroM, sizeof(uint64_t) * WORD_N);
   scalar[0] = MAXIMUM_NUMBER_OF_ISOGENY_TO_BE_COMPUTED * kth_proc;
   
   if(kth_proc == (NUMBER_OF_CORES - 1))
      MAXIMUM_NUMBER_OF_ISOGENY_TO_BE_COMPUTED += VAR_N_DIVIDED_BY_DEGREE_PLUS_1 % NUMBER_OF_CORES;
   
   // Three possible kernels
   Point point_0, point_1, point_2;
   // point_Pt + [scalar]point_Qt
   DBLADD(&point_0, point_Tt, scalar, curve_E1.A);
   Pt_ADD(&point_0, point_0, point_St, curve_E1.A);
   // point_Pt + [scalar + DEGREE^(e_DIVIDED_BY_2 - 1)]point_Qt
   Pt_ADD(&point_1, point_0, order_d_point_Tt, curve_E1.A);
   // [2 * scalar]point_Pt + point_Qt
   DBLADD(&point_2, point_St, scalar, curve_E1.A);
   Pt_DSM(&point_2, point_2, curve_E1.A);
   Pt_ADD(&point_2, point_2, point_Tt, curve_E1.A);
   // Three possible isogenous curves
   Curve curve_0, curve_1, curve_2;
   // Three possible j-invariants
   uint64_t jinv_0[2][WORD_N], jinv_1[2][WORD_N], jinv_2[2][WORD_N];
   
   uint64_t ACC = 0;
   uint64_t COEF_K = scalar[0];
   
   int local_collision;
   var_j = 0;
   while( (var_j < MAXIMUM_NUMBER_OF_ISOGENY_TO_BE_COMPUTED) && (*BOOL != 0x00) )
   {
      // ------------------------------------------------------------------------
      get_isogenous_curve(&curve_0, point_0, Pt_MDSM, get_isog, eval_isog, curve_E1, e_DIVIDED_BY_2);
      jInvariant(jinv_0, curve_0);
      for(local_i = 0; local_i < NUMBER_OF_CORES; local_i++)
      {
	 local_collision = binarySearch(TREE_LEAVES_ROOTED_AT_E0[local_i], jinv_0, 0, CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[local_i]);
	 if( local_collision != -1 )
	    break;
      };
      if( local_collision != -1 )
      {
	 // A collision was found!
	 COLLISION[0][0] = TREE_LEAVES_ROOTED_AT_E0[local_i][local_collision].point;
	 COLLISION[1][0] = (COEF_K << 2) ^ 0x00;
	 *BOOL = 0x00;
	 break;
      }
      Pt_ADD(&point_0, point_0, point_Tt, curve_E1.A);
      ACC += 1;
      
      // ------------------------------------------------------------------------
      get_isogenous_curve(&curve_1, point_1, Pt_MDSM, get_isog, eval_isog, curve_E1, e_DIVIDED_BY_2);
      jInvariant(jinv_1, curve_1);
      for(local_i = 0; local_i < NUMBER_OF_CORES; local_i++)
      {
	 local_collision = binarySearch(TREE_LEAVES_ROOTED_AT_E0[local_i], jinv_1, 0, CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[local_i]);
	 if( local_collision != -1 )
	    break;
      };
      if( local_collision != -1 )
      {
	 // A collision was found!
	 COLLISION[0][0] = TREE_LEAVES_ROOTED_AT_E0[local_i][local_collision].point;
	 COLLISION[1][0] = (COEF_K << 2) ^ 0x01;
	 *BOOL = 0x00;
	 break;
      }
      Pt_ADD(&point_1, point_1, point_Tt, curve_E1.A);
      ACC += 1;
      
      // ------------------------------------------------------------------------
      get_isogenous_curve(&curve_2, point_2, Pt_MDSM, get_isog, eval_isog, curve_E1, e_DIVIDED_BY_2);
      jInvariant(jinv_2, curve_2);
      for(local_i = 0; local_i < NUMBER_OF_CORES; local_i++)
      {
	 local_collision = binarySearch(TREE_LEAVES_ROOTED_AT_E0[local_i], jinv_2, 0, CURRENT_SIZE_OF_THE_TREE_LEAVES_ROOTED_AT_E0[local_i]);
	 if( local_collision != -1 )
	    break;
      };
      if( local_collision != -1 )
      {
	 // A collision was found!
	 COLLISION[0][0] = TREE_LEAVES_ROOTED_AT_E0[local_i][local_collision].point;
	 COLLISION[1][0] = (COEF_K << 2) ^ 0x02;
	 *BOOL = 0x00;
	 break;
      }
      Pt_ADD(&point_2, point_2, d_times_St, curve_E1.A);
      ACC += 1;
      
      COEF_K += 1;
      var_j += 1;
   }
   local_end = get_cycles();
   RUNNING_TIME[kth_proc] += ACC;
   
   return local_end - local_start;
};
