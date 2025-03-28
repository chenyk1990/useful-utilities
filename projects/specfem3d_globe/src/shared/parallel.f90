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

! Dimitri Komatitsch, July 2014, CNRS Marseille, France:
! added the ability to run several calculations (several earthquakes)
! in an embarrassingly-parallel fashion from within the same run;
! this can be useful when using a very large supercomputer to compute
! many earthquakes in a catalog, in which case it can be better from
! a batch job submission point of view to start fewer and much larger jobs,
! each of them computing several earthquakes in parallel.
! To turn that option on, set parameter NUMBER_OF_SIMULTANEOUS_RUNS to a value greater than 1 in the Par_file.
! To implement that, we create NUMBER_OF_SIMULTANEOUS_RUNS MPI sub-communicators,
! each of them being labeled "my_local_mpi_comm_world", and we use them
! in all the routines in "src/shared/parallel.f90", except in MPI_ABORT() because in that case
! we need to kill the entire run.
! When that option is on, of course the number of processor cores used to start
! the code in the batch system must be a multiple of NUMBER_OF_SIMULTANEOUS_RUNS,
! all the individual runs must use the same number of processor cores,
! which as usual is NPROC in the input file DATA/Par_file,
! and thus the total number of processor cores to request from the batch system
! should be NUMBER_OF_SIMULTANEOUS_RUNS * NPROC.
! All the runs to perform must be placed in directories called run0001, run0002, run0003 and so on
! (with exactly four digits).

!-------------------------------------------------------------------------------------------------
!
! Parallel routines.  All MPI calls belong in this file!
!
!-------------------------------------------------------------------------------------------------

module my_mpi

! main parameter module for specfem simulations

  use mpi

  implicit none

  integer :: my_local_mpi_comm_world, my_local_mpi_comm_for_bcast

end module my_mpi

!-------------------------------------------------------------------------------------------------
!
! MPI wrapper functions
!
!-------------------------------------------------------------------------------------------------

  subroutine init_mpi()

  use my_mpi
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer :: myrank,ier

! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)
  if (ier /= 0 ) stop 'Error initializing MPI'

  ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS and BROADCAST_SAME_MESH_AND_MODEL are read before calling world_split()
  ! thus read the parameter file
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
  if (myrank == 0) then
    call open_parameter_file_from_master_only(ier)
    ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS and BROADCAST_SAME_MESH_AND_MODEL are read
    call read_value_integer(NUMBER_OF_SIMULTANEOUS_RUNS, 'NUMBER_OF_SIMULTANEOUS_RUNS', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter NUMBER_OF_SIMULTANEOUS_RUNS'
    call read_value_logical(BROADCAST_SAME_MESH_AND_MODEL, 'BROADCAST_SAME_MESH_AND_MODEL', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter BROADCAST_SAME_MESH_AND_MODEL'
    ! close parameter file
    call close_parameter_file()
  endif

  ! broadcast parameters read from master to all processes
  my_local_mpi_comm_world = MPI_COMM_WORLD

  call bcast_all_singlei(NUMBER_OF_SIMULTANEOUS_RUNS)
  call bcast_all_singlel(BROADCAST_SAME_MESH_AND_MODEL)

! create sub-communicators if needed, if running more than one earthquake from the same job
  call world_split()

  end subroutine init_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_mpi()

  use my_mpi

  implicit none

  integer :: ier

! close sub-communicators if needed, if running more than one earthquake from the same job
  call world_unsplit()

! stop all the MPI processes, and exit
  ! do NOT remove the barrier here, it is critical in order for the failsafe mechanism to work fine when it is activated
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  call MPI_FINALIZE(ier)
  if (ier /= 0) stop 'Error finalizing MPI'

  end subroutine finalize_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine abort_mpi()

  use my_mpi
  use constants, only: MAX_STRING_LEN,mygroup
  use shared_input_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,USE_FAILSAFE_MECHANISM

  implicit none

  integer :: my_local_rank,my_global_rank,ier
  logical :: run_file_exists
  character(len=MAX_STRING_LEN) :: filename

  ! get my local rank and my global rank (in the case of simultaneous jobs, for which we split
  ! the MPI communicator, they will be different; otherwise they are the same)
  call world_rank(my_local_rank)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_global_rank,ier)

  ! write a stamp file to disk to let the user know that the run failed
  if(NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
    ! notifies which run directory failed
    write(filename,"('run',i4.4,'_failed')") mygroup + 1
    inquire(file=trim(filename), exist=run_file_exists)
    if (run_file_exists) then
      open(unit=9765,file=trim(filename),status='old',position='append',action='write',iostat=ier)
    else
      open(unit=9765,file=trim(filename),status='new',action='write',iostat=ier)
    endif
    if (ier == 0) then
      write(9765,*) 'run ',mygroup+1,' with local rank ',my_local_rank,' and global rank ',my_global_rank,' failed'
      close(9765)
    endif

    ! notifies which rank failed
    write(filename,"('run_with_local_rank_',i8.8,'and_global_rank_',i8.8,'_failed')") my_local_rank,my_global_rank
    open(unit=9765,file=trim(filename),status='unknown',action='write')
    write(9765,*) 'run with local rank ',my_local_rank,' and global rank ',my_global_rank,' failed'
    close(9765)
  else
    ! note: we already output an OUTPUT_FILES/error_message***.txt file for each failed rank (single runs)
    ! debug
    !write(filename,"('run_with_local_rank_',i8.8,'_failed')") my_local_rank
    !open(unit=9765,file=filename,status='unknown',action='write')
    !write(9765,*) 'run with local rank ',my_local_rank,' failed'
    !close(9765)
  endif

  ! in case of a large number of simultaneous runs, if one fails we may want that one to just call MPI_FINALIZE() and wait
  ! until all the others are finished instead of calling MPI_ABORT(), which would instead kill all the runs,
  ! including all the successful ones
  if(USE_FAILSAFE_MECHANISM .and. NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
    ! do NOT remove the barrier here, it is critical in order to let other runs finish before calling MPI_FINALIZE
    call MPI_BARRIER(MPI_COMM_WORLD,ier)
    call MPI_FINALIZE(ier)
    if (ier /= 0) stop 'Error finalizing MPI'
  else
    ! note: MPI_ABORT does not return, it makes the program exit with an error code of 30
    call MPI_ABORT(MPI_COMM_WORLD,30,ier)
    stop 'error, program ended in exit_MPI'
  endif

  end subroutine abort_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine synchronize_all()

  use my_mpi

  implicit none

  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(my_local_mpi_comm_world, ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'

  end subroutine synchronize_all

!
!-------------------------------------------------------------------------------------------------
!

  subroutine synchronize_all_comm(comm)

  use my_mpi

  implicit none

  integer,intent(in) :: comm

  ! local parameters
  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(comm,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes for specified communicator'

  end subroutine synchronize_all_comm

!
!-------------------------------------------------------------------------------------------------
!

  double precision function wtime()

  use my_mpi

  implicit none

  wtime = MPI_WTIME()

  end function wtime

!
!-------------------------------------------------------------------------------------------------
!

  integer function null_process()

  use my_mpi

  implicit none

  null_process = MPI_PROC_NULL

  end function null_process

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_request(request,flag_result_test)

  use my_mpi

  implicit none

  integer :: request
  logical :: flag_result_test

  integer :: ier

  call MPI_TEST(request,flag_result_test,MPI_STATUS_IGNORE,ier)

  end subroutine test_request

!
!-------------------------------------------------------------------------------------------------
!

  subroutine wait_req(req)

  use my_mpi

  implicit none

  integer :: req

  integer :: ier

  call mpi_wait(req,MPI_STATUS_IGNORE,ier)

  end subroutine wait_req


!-------------------------------------------------------------------------------------------------
!
! MPI broadcasting helper
!
!-------------------------------------------------------------------------------------------------

!
!---- broadcast using the default communicator for the whole run
!

  subroutine bcast_iproc_i(buffer,iproc)

  use my_mpi

  implicit none

  integer :: iproc
  integer :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_INTEGER,iproc,my_local_mpi_comm_world,ier)

  end subroutine bcast_iproc_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_i(buffer, countval)

  use my_mpi

  implicit none

  integer :: countval
  integer, dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlei(buffer)

  use my_mpi

  implicit none

  integer :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlel(buffer)

  use my_mpi

  implicit none

  logical :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singlel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_cr(buffer, countval)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: countval
  real(kind=CUSTOM_REAL), dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlecr(buffer)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singlecr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_r(buffer, countval)

  use my_mpi

  implicit none

  integer :: countval
  real, dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_r

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singler(buffer)

  use my_mpi

  implicit none

  real :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_REAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singler

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_dp(buffer, countval)

  use my_mpi

  implicit none

  integer :: countval
  double precision, dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singledp(buffer)

  use my_mpi

  implicit none

  double precision :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singledp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_ch(buffer, countval)

  use my_mpi

  implicit none

  integer :: countval
  character(len=countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_ch

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_ch_array(buffer,ndim,countval)

  use my_mpi

  implicit none

  integer :: countval,ndim
  character(len=countval),dimension(ndim) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,ndim*countval,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_ch_array

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_ch_array2(buffer,ndim1,ndim2,countval)

  use my_mpi

  implicit none

  integer :: countval,ndim1,ndim2
  character(len=countval),dimension(ndim1,ndim2) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,ndim1*ndim2*countval,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_ch_array2

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_l(buffer, countval)

  use my_mpi

  implicit none

  integer :: countval
  logical,dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_l


!
!---- broadcast using the communicator to send the mesh and model to other simultaneous runs
!

  subroutine bcast_all_i_for_database(buffer, countval)

  use my_mpi
  use shared_parameters,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  integer :: buffer

  integer ier

  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_i_for_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_l_for_database(buffer, countval)

  use my_mpi
  use shared_parameters,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  logical :: buffer

  integer ier

  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_l_for_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_cr_for_database(buffer, countval)

  use my_mpi
  use constants,only: CUSTOM_REAL
  use shared_parameters,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  include "precision.h"

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  real(kind=CUSTOM_REAL) :: buffer

  integer ier

  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_cr_for_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_dp_for_database(buffer, countval)

  use my_mpi
  use shared_parameters,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  double precision :: buffer

  integer ier

  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_dp_for_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_r_for_database(buffer, countval)

  use my_mpi
  use shared_parameters,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  real :: buffer

  integer ier

  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_r_for_database

!
!---- broadcast using MPI_COMM_WORLD
!

!  subroutine bcast_all_singlei_world(buffer)
!  end subroutine bcast_all_singlei_world

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine bcast_all_singlel_world(buffer)
!  end subroutine bcast_all_singlel_world

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine bcast_all_singledp_world(buffer)
!  end subroutine bcast_all_singledp_world

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine bcast_all_string_world(buffer)
!  end subroutine bcast_all_string_world


!-------------------------------------------------------------------------------------------------
!
! MPI math helper
!
!-------------------------------------------------------------------------------------------------

  subroutine min_all_i(sendbuf, recvbuf)

  use my_mpi

  implicit none

  integer:: sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MIN,0,my_local_mpi_comm_world,ier)

  end subroutine min_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_i(sendbuf, recvbuf)

  use my_mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MAX,0,my_local_mpi_comm_world,ier)

  end subroutine max_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_allreduce_i(buffer,countval)

  use my_mpi

  implicit none

  integer :: countval
  integer,dimension(countval),intent(inout) :: buffer

  ! local parameters
  integer :: ier
  integer,dimension(countval) :: send

  ! seems not to be supported on all kind of MPI implementations...
  !! DK DK: yes, I confirm, using MPI_IN_PLACE is tricky
  !! DK DK (see the answer at http://stackoverflow.com/questions/17741574/in-place-mpi-reduce-crashes-with-openmpi
  !! DK DK      for how to use it right)
  !call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, countval, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)

  send(:) = buffer(:)

  call MPI_ALLREDUCE(send, buffer, countval, MPI_INTEGER, MPI_MAX, my_local_mpi_comm_world, ier)
  if (ier /= 0 ) stop 'Allreduce to get max values failed.'

  end subroutine max_allreduce_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MIN,0,my_local_mpi_comm_world,ier)

  end subroutine min_all_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine min_all_all_cr(sendbuf, recvbuf)
!  end subroutine min_all_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MAX,0,my_local_mpi_comm_world,ier)

  end subroutine max_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_allreduce_cr(sendbuf, recvbuf)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MAX,my_local_mpi_comm_world,ier)

  end subroutine max_allreduce_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine max_all_all_cr(sendbuf, recvbuf)
!  end subroutine max_all_all_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine min_all_dp(sendbuf, recvbuf)
!  end subroutine min_all_dp

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine max_all_dp(sendbuf, recvbuf)
!  end subroutine max_all_dp

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine max_all_all_dp(sendbuf, recvbuf)
!  end subroutine max_all_all_dp

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine maxloc_all_dp(sendbuf, recvbuf)
!  end subroutine maxloc_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine any_all_l(sendbuf, recvbuf)

  use my_mpi

  implicit none

  logical :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_LOGICAL,MPI_LOR,my_local_mpi_comm_world,ier)

  end subroutine any_all_l

!
!-------------------------------------------------------------------------------------------------
!

! MPI summations

  subroutine sum_all_i(sendbuf, recvbuf)

  use my_mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_i

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine sum_all_all_i(sendbuf, recvbuf)
!  end subroutine sum_all_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine sum_all_all_cr(sendbuf, recvbuf)
!  end subroutine sum_all_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_dp

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine sum_all_1Darray_dp(sendbuf, recvbuf, nx)
!  end subroutine sum_all_1Darray_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_3Darray_dp(sendbuf, recvbuf, nx,ny,nz)

  use my_mpi

  implicit none

  integer :: nx,ny,nz
  double precision, dimension(nx,ny,nz) :: sendbuf, recvbuf
  integer :: ier

  ! this works only if the arrays are contiguous in memory (which is always the case for static arrays, as used in the code)
  call MPI_REDUCE(sendbuf,recvbuf,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_3Darray_dp


!-------------------------------------------------------------------------------------------------
!
! Send/Receive MPI
!
!-------------------------------------------------------------------------------------------------

! asynchronuous send/receive

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf

  integer :: ier

  call MPI_ISEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,my_local_mpi_comm_world,req,ier)

  end subroutine isend_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine isend_i(sendbuf, sendcount, dest, sendtag, req)
!  end subroutine isend_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine isend_dp(sendbuf, sendcount, dest, sendtag, req)

  use my_mpi

  implicit none

  integer :: sendcount, dest, sendtag, req
  double precision, dimension(sendcount) :: sendbuf

  integer :: ier

  call MPI_ISEND(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,sendtag,my_local_mpi_comm_world,req,ier)

  end subroutine isend_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_IRECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag,my_local_mpi_comm_world,req,ier)

  end subroutine irecv_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine irecv_dp(recvbuf, recvcount, dest, recvtag, req)

  use my_mpi

  implicit none

  integer :: recvcount, dest, recvtag, req
  double precision, dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_IRECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,dest,recvtag,my_local_mpi_comm_world,req,ier)

  end subroutine irecv_dp

!
!-------------------------------------------------------------------------------------------------
!

! synchronuous/blocking send/receive

  subroutine recv_i(recvbuf, recvcount, dest, recvtag)

  use my_mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  integer,dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_cr(recvbuf, recvcount, dest, recvtag)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: dest,recvtag
  integer :: recvcount
  real(kind=CUSTOM_REAL),dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_dp(recvbuf, recvcount, dest, recvtag)

  use my_mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  double precision,dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_ch(recvbuf, recvcount, dest, recvtag)

  use my_mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  character(len=recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_CHARACTER,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_ch

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_singlei(recvbuf, dest, recvtag)

  use my_mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,1,MPI_INTEGER,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_singlel(recvbuf, dest, recvtag)

  use my_mpi

  implicit none

  integer :: dest,recvtag
  logical :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,1,MPI_LOGICAL,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_singlel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_ch(sendbuf, sendcount, dest, sendtag)

  use my_mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  character(len=sendcount) :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_CHARACTER,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_ch


!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_i(sendbuf, sendcount, dest, sendtag)

  use my_mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  integer,dimension(sendcount):: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_singlei(sendbuf, dest, sendtag)

  use my_mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,1,MPI_INTEGER,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_singlel(sendbuf, dest, sendtag)

  use my_mpi

  implicit none

  integer :: dest,sendtag
  logical :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,1,MPI_LOGICAL,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_singlel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_cr(sendbuf, sendcount, dest, sendtag)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: dest,sendtag
  integer :: sendcount
  real(kind=CUSTOM_REAL),dimension(sendcount):: sendbuf
  integer :: ier

  call MPI_SEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_dp(sendbuf, sendcount, dest, sendtag)

  use my_mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  double precision,dimension(sendcount):: sendbuf
  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sendrecv_cr(sendbuf, sendcount, dest, sendtag, &
                         recvbuf, recvcount, source, recvtag)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: sendcount, recvcount, dest, sendtag, source, recvtag
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_SENDRECV(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag, &
                    recvbuf,recvcount,CUSTOM_MPI_TYPE,source,recvtag, &
                    my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine sendrecv_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sendrecv_dp(sendbuf, sendcount, dest, sendtag, &
                         recvbuf, recvcount, source, recvtag)

  use my_mpi

  implicit none

  integer :: sendcount, recvcount, dest, sendtag, source, recvtag
  double precision, dimension(sendcount) :: sendbuf
  double precision, dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_SENDRECV(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,sendtag, &
                    recvbuf,recvcount,MPI_DOUBLE_PRECISION,source,recvtag, &
                    my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine sendrecv_dp

!-------------------------------------------------------------------------------------------------
!
! MPI gather helper
!
!-------------------------------------------------------------------------------------------------

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi

  implicit none

  integer :: sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_INTEGER, &
                  recvbuf,recvcount,MPI_INTEGER, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_singlei(sendbuf, recvbuf, NPROC)

  use my_mpi

  implicit none

  integer :: NPROC
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,1,MPI_INTEGER, &
                  recvbuf,1,MPI_INTEGER, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_cr(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: sendcnt, recvcount, NPROC
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                  recvbuf,recvcount,CUSTOM_MPI_TYPE, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi

  implicit none

  integer :: sendcnt, recvcount, NPROC
  double precision, dimension(sendcnt) :: sendbuf
  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
                  recvbuf,recvcount,MPI_DOUBLE_PRECISION, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gatherv_all_i(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  use my_mpi

  implicit none

  include "precision.h"

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcounttot) :: recvbuf

  integer :: ier

  call MPI_GATHERV(sendbuf,sendcnt,MPI_INTEGER, &
                  recvbuf,recvcount,recvoffset,MPI_INTEGER, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gatherv_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gatherv_all_cr(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcounttot) :: recvbuf

  integer :: ier

  call MPI_GATHERV(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                  recvbuf,recvcount,recvoffset,CUSTOM_MPI_TYPE, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gatherv_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gatherv_all_r(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  use my_mpi

  implicit none

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real, dimension(sendcnt) :: sendbuf
  real, dimension(recvcounttot) :: recvbuf

  integer :: ier

  call MPI_GATHERV(sendbuf,sendcnt,MPI_REAL, &
                  recvbuf,recvcount,recvoffset,MPI_REAL, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gatherv_all_r

!
!-------------------------------------------------------------------------------------------------
!

  subroutine all_gather_all_i(sendbuf, recvbuf, NPROC)

  use constants
  use my_mpi

  implicit none

  integer :: NPROC
  integer :: sendbuf
  integer, dimension(NPROC) :: recvbuf

  integer :: ier

  call MPI_Allgather(sendbuf,1,MPI_INTEGER, &
                  recvbuf,1,MPI_INTEGER, &
                  my_local_mpi_comm_world,ier)

  end subroutine all_gather_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine all_gather_all_r(sendbuf, sendcnt, recvbuf, recvcnt, recvoffset, dim1, NPROC)

  use constants
  use my_mpi

  implicit none

  integer :: sendcnt, dim1, NPROC

  real, dimension(NPROC) :: sendbuf
  real, dimension(dim1, NPROC) :: recvbuf

  integer, dimension(NPROC) :: recvoffset, recvcnt

  integer :: ier

  call MPI_Allgatherv(sendbuf,sendcnt,MPI_REAL, &
                  recvbuf,recvcnt,recvoffset,MPI_REAL, &
                  my_local_mpi_comm_world,ier)

  end subroutine all_gather_all_r

!
!-------------------------------------------------------------------------------------------------
!

  subroutine all_gather_all_ch(sendbuf, sendcnt, recvbuf, recvcnt, recvoffset, dim1, dim2, NPROC)

  use constants
  use my_mpi

  implicit none

  integer :: sendcnt, dim1, dim2, NPROC

  character(len=dim2), dimension(NPROC) :: sendbuf
  character(len=dim2), dimension(dim1, NPROC) :: recvbuf

  integer, dimension(NPROC) :: recvoffset, recvcnt

  integer :: ier

  call MPI_Allgatherv(sendbuf,sendcnt,MPI_CHARACTER, &
                  recvbuf,recvcnt,recvoffset,MPI_CHARACTER, &
                  my_local_mpi_comm_world,ier)

  end subroutine all_gather_all_ch

!
!-------------------------------------------------------------------------------------------------
!

  subroutine scatter_all_singlei(sendbuf, recvbuf, NPROC)

  use my_mpi

  implicit none

  integer :: NPROC
  integer, dimension(0:NPROC-1) :: sendbuf
  integer :: recvbuf

  integer :: ier

  call MPI_Scatter(sendbuf, 1, MPI_INTEGER, &
                   recvbuf, 1, MPI_INTEGER, &
                   0, my_local_mpi_comm_world, ier)

  end subroutine scatter_all_singlei

!-------------------------------------------------------------------------------------------------
!
! MPI world helper
!
!-------------------------------------------------------------------------------------------------

  subroutine world_size(sizeval)

  use my_mpi

  implicit none

  integer,intent(out) :: sizeval

  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(my_local_mpi_comm_world,sizeval,ier)
  if (ier /= 0 ) stop 'Error getting MPI world size'

  end subroutine world_size

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_rank(rank)

  use my_mpi

  implicit none

  integer,intent(out) :: rank

  ! local parameters
  integer :: ier

  call MPI_COMM_RANK(my_local_mpi_comm_world,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'

  end subroutine world_rank

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_rank_comm(rank,comm)

  use my_mpi

  implicit none

  integer,intent(out) :: rank
  integer,intent(in) :: comm

  ! local parameters
  integer :: ier

  call MPI_COMM_RANK(comm,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'

  end subroutine world_rank_comm

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_duplicate(comm)

  use my_mpi

  implicit none

  integer,intent(out) :: comm
  integer :: ier

  ! note: see http://www.mpich.org/static/docs/v3.1/www3/MPI_Comm_dup.html
  ! "..
  ! This routine is used to create a new communicator that has a new communication context but contains
  ! the same group of processes as the input communicator. Since all MPI communication is performed within
  ! a communicator (specifies as the group of processes plus the context), this routine provides an effective way
  ! to create a private communicator for use by a software module or library.
  !
  ! In particular, no library routine should use MPI_COMM_WORLD as the communicator;
  ! instead, a duplicate of a user-specified communicator should always be used."

  call MPI_COMM_DUP(my_local_mpi_comm_world,comm,ier)
  if (ier /= 0 ) stop 'Error duplicating my_local_mpi_comm_world communicator'

  end subroutine world_duplicate

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_get_comm(comm)

  use my_mpi

  implicit none

  integer,intent(out) :: comm

  comm = my_local_mpi_comm_world

  end subroutine world_get_comm

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_get_comm_self(comm)

  use my_mpi

  implicit none

  integer,intent(out) :: comm

  comm = MPI_COMM_SELF

  end subroutine world_get_comm_self


!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_get_info_null(info)

  use my_mpi

  implicit none

  integer,intent(out) :: info

  info = MPI_INFO_NULL

  end subroutine world_get_info_null

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_comm_free(comm)

  use my_mpi

  implicit none

  integer,intent(inout) :: comm

  ! local parameters
  integer :: ier

  call MPI_Comm_free(comm,ier)
  if (ier /= 0 ) stop 'Error freeing MPI communicator'

  end subroutine world_comm_free


!
!-------------------------------------------------------------------------------------------------
!

! create sub-communicators if needed, if running more than one earthquake from the same job.
! create a sub-communicator for each independent run;
! if there is a single run to do, then just copy the default communicator to the new one
  subroutine world_split()

  use my_mpi

  use constants,only: MAX_STRING_LEN,OUTPUT_FILES_BASE, &
    IMAIN,ISTANDARD_OUTPUT,mygroup,I_should_read_the_database

  use shared_parameters,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL,OUTPUT_FILES

  implicit none

  integer :: sizeval,myrank,ier,key,my_group_for_bcast,my_local_rank_for_bcast,NPROC

  character(len=MAX_STRING_LEN) :: path_to_add

  if (NUMBER_OF_SIMULTANEOUS_RUNS <= 0) stop 'NUMBER_OF_SIMULTANEOUS_RUNS <= 0 makes no sense'

  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeval,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mod(sizeval,NUMBER_OF_SIMULTANEOUS_RUNS) /= 0) then
    if (myrank == 0) print *,'Error: the number of MPI processes ',sizeval, &
                            ' is not a multiple of NUMBER_OF_SIMULTANEOUS_RUNS = ',NUMBER_OF_SIMULTANEOUS_RUNS
    stop 'the number of MPI processes is not a multiple of NUMBER_OF_SIMULTANEOUS_RUNS'
  endif

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. IMAIN == ISTANDARD_OUTPUT) &
    stop 'must not have IMAIN == ISTANDARD_OUTPUT when NUMBER_OF_SIMULTANEOUS_RUNS > 1 otherwise output to screen is mingled'

  OUTPUT_FILES = OUTPUT_FILES_BASE(1:len_trim(OUTPUT_FILES_BASE))

  if (NUMBER_OF_SIMULTANEOUS_RUNS == 1) then

    my_local_mpi_comm_world = MPI_COMM_WORLD

! no broadcast of the mesh and model databases to other runs in that case
    my_group_for_bcast = 0
    my_local_mpi_comm_for_bcast = MPI_COMM_NULL

  else

!--- create a subcommunicator for each independent run

    NPROC = sizeval / NUMBER_OF_SIMULTANEOUS_RUNS

!   create the different groups of processes, one for each independent run
    mygroup = myrank / NPROC
    key = myrank
    if (mygroup < 0 .or. mygroup > NUMBER_OF_SIMULTANEOUS_RUNS-1) stop 'invalid value of mygroup'

!   build the sub-communicators
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, mygroup, key, my_local_mpi_comm_world, ier)
    if (ier /= 0) stop 'error while trying to create the sub-communicators'

!   add the right directory for that run
!   (group numbers start at zero, but directory names start at run0001, thus we add one)
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    OUTPUT_FILES = path_to_add(1:len_trim(path_to_add))//OUTPUT_FILES(1:len_trim(OUTPUT_FILES))

!--- create a subcommunicator to broadcast the identical mesh and model databases if needed
    if (BROADCAST_SAME_MESH_AND_MODEL) then

      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
!     to broadcast the model, split along similar ranks per run instead
      my_group_for_bcast = mod(myrank,NPROC)
      key = myrank
      if (my_group_for_bcast < 0 .or. my_group_for_bcast > NPROC-1) stop 'invalid value of my_group_for_bcast'

!     build the sub-communicators
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, my_group_for_bcast, key, my_local_mpi_comm_for_bcast, ier)
      if (ier /= 0) stop 'error while trying to create the sub-communicators'

!     see if that process will need to read the mesh and model database and then broadcast it to others
      call MPI_COMM_RANK(my_local_mpi_comm_for_bcast,my_local_rank_for_bcast,ier)
      if (my_local_rank_for_bcast > 0) I_should_read_the_database = .false.

    else

! no broadcast of the mesh and model databases to other runs in that case
      my_group_for_bcast = 0
      my_local_mpi_comm_for_bcast = MPI_COMM_NULL

    endif

  endif

  end subroutine world_split

!
!-------------------------------------------------------------------------------------------------
!

! close sub-communicators if needed, if running more than one earthquake from the same job.
  subroutine world_unsplit()

  use my_mpi
  use shared_parameters,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer :: ier

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
    call MPI_COMM_FREE(my_local_mpi_comm_world,ier)
    if (BROADCAST_SAME_MESH_AND_MODEL) call MPI_COMM_FREE(my_local_mpi_comm_for_bcast,ier)
  endif

  end subroutine world_unsplit

