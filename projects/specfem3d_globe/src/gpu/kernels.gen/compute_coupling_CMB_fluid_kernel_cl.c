//note: please do not modify this file manually!
//      this file has been generated automatically by BOAST version 1.2.0
//      by: make boast_kernels

/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
*/

const char * compute_coupling_CMB_fluid_kernel_program = "\
inline void atomicAdd(volatile __global float *source, const float val) {\n\
  union {\n\
    unsigned int iVal;\n\
    float fVal;\n\
  } res, orig;\n\
  do {\n\
    orig.fVal = *source;\n\
    res.fVal = orig.fVal + val;\n\
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, orig.iVal, res.iVal) != orig.iVal);\n\
}\n\
#ifndef INDEX2\n\
#define INDEX2(isize,i,j) i + isize*j\n\
#endif\n\
#ifndef INDEX3\n\
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)\n\
#endif\n\
#ifndef INDEX4\n\
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))\n\
#endif\n\
#ifndef INDEX5\n\
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))\n\
#endif\n\
\n\
#ifndef NDIM\n\
#define NDIM 3\n\
#endif\n\
#ifndef NGLLX\n\
#define NGLLX 5\n\
#endif\n\
#ifndef NGLL2\n\
#define NGLL2 25\n\
#endif\n\
#ifndef NGLL3\n\
#define NGLL3 125\n\
#endif\n\
#ifndef NGLL3_PADDED\n\
#define NGLL3_PADDED 128\n\
#endif\n\
#ifndef N_SLS\n\
#define N_SLS 3\n\
#endif\n\
#ifndef IREGION_CRUST_MANTLE\n\
#define IREGION_CRUST_MANTLE 1\n\
#endif\n\
#ifndef IREGION_INNER_CORE\n\
#define IREGION_INNER_CORE 3\n\
#endif\n\
#ifndef IFLAG_IN_FICTITIOUS_CUBE\n\
#define IFLAG_IN_FICTITIOUS_CUBE 11\n\
#endif\n\
#ifndef R_EARTH_KM\n\
#define R_EARTH_KM 6371.0f\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_INNER_CORE\n\
#define COLORING_MIN_NSPEC_INNER_CORE 1000\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_OUTER_CORE\n\
#define COLORING_MIN_NSPEC_OUTER_CORE 1000\n\
#endif\n\
#ifndef BLOCKSIZE_TRANSFER\n\
#define BLOCKSIZE_TRANSFER 256\n\
#endif\n\
\n\
__kernel void compute_coupling_CMB_fluid_kernel(const __global float * displ_crust_mantle, __global float * accel_crust_mantle, const __global float * accel_outer_core, const __global int * ibool_crust_mantle, const __global int * ibelm_bottom_crust_mantle, const __global float * normal_top_outer_core, const __global float * jacobian2D_top_outer_core, const __global float * wgllwgll_xy, const __global int * ibool_outer_core, const __global int * ibelm_top_outer_core, const float RHO_TOP_OC, const float minus_g_cmb, const int GRAVITY, const int NSPEC2D_BOTTOM_CM){\n\
  int i;\n\
  int j;\n\
  int k;\n\
  int iface;\n\
  int k_corresp;\n\
  int iglob_oc;\n\
  int iglob_cm;\n\
  float pressure;\n\
  int ispec;\n\
  int ispec_selected;\n\
  float nx;\n\
  float ny;\n\
  float nz;\n\
  float weight;\n\
  i = get_local_id(0);\n\
  j = get_local_id(1);\n\
  iface = get_group_id(0) + (get_num_groups(0)) * (get_group_id(1));\n\
  if (iface < NSPEC2D_BOTTOM_CM) {\n\
    ispec = ibelm_bottom_crust_mantle[iface] - (1);\n\
    ispec_selected = ibelm_top_outer_core[iface] - (1);\n\
    k = 0;\n\
    k_corresp = NGLLX - (1);\n\
    iglob_oc = ibool_outer_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k_corresp, ispec_selected)] - (1);\n\
    nx = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 0, i, j, iface)];\n\
    ny = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 1, i, j, iface)];\n\
    nz = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 2, i, j, iface)];\n\
    weight = (jacobian2D_top_outer_core[INDEX3(NGLLX, NGLLX, i, j, iface)]) * (wgllwgll_xy[INDEX2(NGLLX, i, j)]);\n\
    iglob_cm = ibool_crust_mantle[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)] - (1);\n\
    if (GRAVITY) {\n\
      pressure = (RHO_TOP_OC) * ((minus_g_cmb) * ((displ_crust_mantle[(iglob_cm) * (3)]) * (nx) + (displ_crust_mantle[(iglob_cm) * (3) + 1]) * (ny) + (displ_crust_mantle[(iglob_cm) * (3) + 2]) * (nz)) - (accel_outer_core[iglob_oc]));\n\
    } else {\n\
      pressure = ( -(RHO_TOP_OC)) * (accel_outer_core[iglob_oc]);\n\
    }\n\
    atomicAdd(accel_crust_mantle + (iglob_cm) * (3) + 0, ((weight) * (nx)) * (pressure));\n\
    atomicAdd(accel_crust_mantle + (iglob_cm) * (3) + 1, ((weight) * (ny)) * (pressure));\n\
    atomicAdd(accel_crust_mantle + (iglob_cm) * (3) + 2, ((weight) * (nz)) * (pressure));\n\
  }\n\
}\n\
";
