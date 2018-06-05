#ifndef _CPU_POISSON_SOLVER_3D_OPEN_OCTREE_POISSON_H_
#define _CPU_POISSON_SOLVER_3D_OPEN_OCTREE_POISSON_H_


void CPU_SolveOpenOctreePoissonRedBlack3_MAC(float* mac_u, float* mac_v, float* mac_w, const bool* leaf0, const bool* leaf1, const bool* leaf2, const bool* leaf3, 
											 const int width, const int height, const int depth, const int maxIter,
											 const int level0_num_red, const int* level0_index_red, const int level0_info_len_red, const int* level0_neighborinfo_red,
											 const int level0_num_black, const int* level0_index_black, const int level0_info_len_black, const int* level0_neighborinfo_black,
											 const int level1_num_red, const int* level1_index_red, const int level1_info_len_red, const int* level1_neighborinfo_red, 
											 const int level1_num_black, const int* level1_index_black, const int level1_info_len_black, const int* level1_neighborinfo_black,
											 const int level2_num_red, const int* level2_index_red, const int level2_info_len_red, const int* level2_neighborinfo_red,
											 const int level2_num_black, const int* level2_index_black, const int level2_info_len_black, const int* level2_neighborinfo_black,
											 const int level3_num_red, const int* level3_index_red, const int level3_info_len_red, const int* level3_neighborinfo_red,
											 const int level3_num_black, const int* level3_index_black, const int level3_info_len_black, const int* level3_neighborinfo_black);


#endif