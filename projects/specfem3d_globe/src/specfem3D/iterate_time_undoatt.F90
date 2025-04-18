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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine iterate_time_undoatt()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_noise
  use specfem_par_movie
  use write_seismograms_mod, only: write_seismograms
  implicit none

  ! local parameters
  integer :: it_temp,seismo_current_temp
  integer :: i,j,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_displ_cm_store_buffer
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_displ_ic_store_buffer
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_displ_oc_store_buffer,b_accel_oc_store_buffer
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: b_noise_surface_movie_buffer
  double precision :: sizeval
  ! timing
  double precision, external :: wtime

  !----  create a Gnuplot script to display the energy curve in log scale
  if (OUTPUT_ENERGY .and. myrank == 0) then
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'plot_energy.gnu',status='unknown',action='write')
    write(IOUT_ENERGY,*) 'set term wxt'
    write(IOUT_ENERGY,*) '#set term postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) '#set output "energy.ps"'
    write(IOUT_ENERGY,*) 'set logscale y'
    write(IOUT_ENERGY,*) 'set xlabel "Time step number"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,'(a152)') '#plot "energy.dat" us 1:2 t ''Kinetic Energy'' w l lc 1, "energy.dat" us 1:3 &
                         &t ''Potential Energy'' w l lc 2, "energy.dat" us 1:4 t ''Total Energy'' w l lc 4'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) '#plot "energy.dat" us 1:2 t ''Kinetic Energy'' w l lc 1'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) '#plot "energy.dat" us 1:3 t ''Potential Energy'' w l lc 2'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) 'plot "energy.dat" us 1:4 t ''Total Energy'' w l lc 4'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    close(IOUT_ENERGY)
  endif

  ! open the file in which we will store the energy curve
  if (OUTPUT_ENERGY .and. myrank == 0) &
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'energy.dat',status='unknown',action='write')

!
!   s t a r t   t i m e   i t e r a t i o n s
!

  ! checks
  if (.not. UNDO_ATTENUATION) return

  ! number of time subsets for time loop
  NSUBSET_ITERATIONS = ceiling( dble(NSTEP)/dble(NT_DUMP_ATTENUATION) )

  ! checks
  if (NSUBSET_ITERATIONS <= 0) call exit_MPI(myrank,'Error invalid number of time subsets for undoing attenuation')

  ! user output
  if (SAVE_FORWARD .or. SIMULATION_TYPE == 3) then
    if (myrank == 0) then
      write(IMAIN,*) 'undoing attenuation:'
      write(IMAIN,*) '  total number of time subsets                     = ',NSUBSET_ITERATIONS
      write(IMAIN,*) '  wavefield snapshots at every NT_DUMP_ATTENUATION = ',NT_DUMP_ATTENUATION
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! allocates buffers
  if (SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
      ! crust/mantle
      ! buffer(NDIM,NGLOB_CRUST_MANTLE_ADJOINT,NT_DUMP_ATTENUATION) in MB
      sizeval = dble(NDIM) * dble(NGLOB_CRUST_MANTLE_ADJOINT) * dble(NT_DUMP_ATTENUATION) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
      write(IMAIN,*) '  size of crust/mantle wavefield buffer per slice = ', sngl(sizeval),'MB'
      ! outer core
      ! buffer(NGLOB_OUTER_CORE_ADJOINT,NT_DUMP_ATTENUATION) in MB
      sizeval = dble(2) * dble(NGLOB_OUTER_CORE_ADJOINT) * dble(NT_DUMP_ATTENUATION) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
      write(IMAIN,*) '  size of   outer core wavefield buffer per slice = ', sngl(sizeval),'MB'
      ! inner core
      ! buffer(NDIM,NGLOB_CRUST_MANTLE_ADJOINT,NT_DUMP_ATTENUATION) in MB
      sizeval = dble(NDIM) * dble(NGLOB_INNER_CORE_ADJOINT) * dble(NT_DUMP_ATTENUATION) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
      write(IMAIN,*) '  size of   inner core wavefield buffer per slice = ', sngl(sizeval),'MB'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    !! DK DK to Daniel, July 2013: in the case of GPU_MODE it would probably be better to leave these arrays on the host
    !! i.e. on the CPU, in order to be able to use all the (unused) memory of the host for them, since they are
    !! (purposely) huge and designed to use almost all the memory available (by carefully optimizing the
    !! value of NT_DUMP_ATTENUATION); when writing to these buffers, it will then be OK to use non-blocking writes
    !! from the device to the host, since we do not reuse the content of these buffers until much later, in a second part
    !! of the run; however when reading back from these buffers, the reads from host to device should then be blocking
    !! because we then use the values read immediately (one value at a time, but to reduce the total number of reads
    !! across the PCI-Express bus we could / should consider reading them 10 by 10 for instance (?) if that fits
    !! in the memory of the GPU
    allocate(b_displ_cm_store_buffer(NDIM,NGLOB_CRUST_MANTLE_ADJOINT,NT_DUMP_ATTENUATION),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_displ_cm_store_buffer')
    allocate(b_displ_oc_store_buffer(NGLOB_OUTER_CORE_ADJOINT,NT_DUMP_ATTENUATION),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_displ_oc_store_buffer')
    allocate(b_accel_oc_store_buffer(NGLOB_OUTER_CORE_ADJOINT,NT_DUMP_ATTENUATION),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_accel_oc_store_buffer')
    allocate(b_displ_ic_store_buffer(NDIM,NGLOB_INNER_CORE_ADJOINT,NT_DUMP_ATTENUATION),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_displ_ic_store_buffer')

    ! noise kernel for source strength (sigma_kernel) needs buffer for reconstructed noise_surface_movie array,
    ! otherwise we need file i/o which will considerably slow down performance
    if (NOISE_TOMOGRAPHY == 3) then
      allocate(b_noise_surface_movie_buffer(NDIM,NGLLX,NGLLY,NSPEC_TOP,NT_DUMP_ATTENUATION),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_noise_surface_movie_buffer')
    endif
  endif

  ! synchronize all processes to make sure everybody is ready to start time loop
  call synchronize_all()
  if (myrank == 0) write(IMAIN,*) 'All processes are synchronized before time loop'

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop (undoatt)...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! create an empty file to monitor the start of the simulation
  if (myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/starttimeloop.txt',status='unknown',action='write')
    write(IOUT,*) 'hello, starting time loop'
    close(IOUT)
  endif

  ! initializes variables for writing seismograms
  seismo_offset = it_begin-1
  seismo_current = 0

  ! initializes time increments
  it = 0

  ! get MPI starting time
  time_start = wtime()

  ! *********************************************************
  ! ************* MAIN LOOP OVER THE TIME STEPS *************
  ! *********************************************************

  if (EXACT_UNDOING_TO_DISK) call setup_exact_undoing_to_disk()

  ! loops over time subsets
  do iteration_on_subset = 1, NSUBSET_ITERATIONS

    ! wavefield storage
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      ! saves forward wavefields
      call save_forward_arrays_undoatt()

    else if (SIMULATION_TYPE == 3) then
      ! reads in last stored forward wavefield
      call read_forward_arrays_undoatt()

      ! note: after reading the restart files of displacement back from disk, recompute the strain from displacement;
      !       this is better than storing the strain to disk as well, which would drastically increase I/O volume
      ! computes strain based on current backward/reconstructed wavefield
      if (COMPUTE_AND_STORE_STRAIN) call it_compute_strain_att_backward()
    endif

    ! time loop within this iteration subset
    select case (SIMULATION_TYPE)
    case (1, 2)
      ! forward and adjoint simulations

      ! increment end of this subset
      if (iteration_on_subset < NSUBSET_ITERATIONS) then
        ! takes full length of subset
        it_subset_end = NT_DUMP_ATTENUATION
      else
        ! loops over remaining steps in last subset
        it_subset_end = NSTEP - (iteration_on_subset-1)*NT_DUMP_ATTENUATION
      endif
      ! checks end index
      if (it_subset_end > NT_DUMP_ATTENUATION) &
        call exit_MPI(myrank,'Error invalid buffer index for undoing attenuation')

      ! subset loop
      do it_of_this_subset = 1, it_subset_end

        it = it + 1

        ! simulation status output and stability check
        if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
          call check_stability()
        endif

        do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then

          if (USE_LDDRK) then
            ! update displacement using Runge-Kutta time scheme
            call update_displ_lddrk()
          else
            ! update displacement using Newmark time scheme
            call update_displ_Newmark()
          endif

          ! acoustic solver for outer core
          ! (needs to be done first, before elastic one)
          call compute_forces_acoustic()

          ! elastic solver for crust/mantle and inner core
          call compute_forces_viscoelastic()

        enddo ! istage

        ! write the seismograms with time shift
        call write_seismograms()

        ! outputs movie files
        call write_movie_output()

        ! first step of noise tomography, i.e., save a surface movie at every time step
        ! modified from the subroutine 'write_movie_surface'
        if (NOISE_TOMOGRAPHY == 1) then
          call noise_save_surface_movie()
        endif

        ! updates VTK window
        if (VTK_MODE) then
          call it_update_vtkwindow()
        endif

      enddo ! subset loop

    case (3)
      ! kernel simulations

      ! intermediate storage of it and seismo_current positions
      it_temp = it
      seismo_current_temp = seismo_current

      ! increment end of this subset
      if (iteration_on_subset == 1) then
        ! loops over remaining steps in last forward subset
        it_subset_end = NSTEP - (NSUBSET_ITERATIONS-1)*NT_DUMP_ATTENUATION
      else
        ! takes full length of subset
        it_subset_end = NT_DUMP_ATTENUATION
      endif
      ! checks end index
      if (it_subset_end > NT_DUMP_ATTENUATION) &
        call exit_MPI(myrank,'Error invalid buffer index for undoing attenuation')

      ! reconstructs forward wavefields based on last stored wavefield data

      ! note: we step forward in time here, starting from last snapshot.
      !       the newly computed, reconstructed forward wavefields (b_displ_..) get stored in buffers.

      ! subset loop
      do it_of_this_subset = 1, it_subset_end

        it = it + 1

        ! simulation status output and stability check
        if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
          call check_stability_backward()
        endif

        do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then

          if (USE_LDDRK) then
            ! update displacement using Runge-Kutta time scheme
            call update_displ_lddrk_backward()
          else
            ! update displacement using Newmark time scheme
            call update_displ_Newmark_backward()
          endif

          ! acoustic solver for outer core
          ! (needs to be done first, before elastic one)
          call compute_forces_acoustic_backward()

          ! elastic solver for crust/mantle and inner core
          call compute_forces_viscoelastic_backward()

        enddo ! istage

        ! transfers wavefields from GPU to CPU for buffering
        if (GPU_MODE) then
          ! daniel debug: check if these transfers could be made async to overlap
          call transfer_b_displ_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,b_displ_crust_mantle,Mesh_pointer)
          call transfer_b_displ_ic_from_device(NDIM*NGLOB_INNER_CORE,b_displ_inner_core,Mesh_pointer)
          call transfer_b_displ_oc_from_device(NGLOB_OUTER_CORE,b_displ_outer_core,Mesh_pointer)
          call transfer_b_accel_oc_from_device(NGLOB_OUTER_CORE,b_accel_outer_core,Mesh_pointer)
        endif

        ! stores wavefield in buffers
        b_displ_cm_store_buffer(:,:,it_of_this_subset) = b_displ_crust_mantle(:,:)
        b_displ_oc_store_buffer(:,it_of_this_subset) = b_displ_outer_core(:)
        b_accel_oc_store_buffer(:,it_of_this_subset) = b_accel_outer_core(:)
        b_displ_ic_store_buffer(:,:,it_of_this_subset) = b_displ_inner_core(:,:)

        ! for noise kernel
        if (NOISE_TOMOGRAPHY == 3) then
          b_noise_surface_movie_buffer(:,:,:,:,it_of_this_subset) = noise_surface_movie(:,:,:,:)
        endif
      enddo ! subset loop

      ! resets current it and seismo_current positions
      it = it_temp
      seismo_current = seismo_current_temp

      ! computes strain based on current adjoint wavefield
      if (COMPUTE_AND_STORE_STRAIN) call it_compute_strain_att()

      ! adjoint wavefield simulation
      do it_of_this_subset = 1, it_subset_end

        ! reads backward/reconstructed wavefield from buffers
        ! note: uses wavefield at corresponding time (NSTEP - it + 1 ), i.e. we have now time-reversed wavefields
        ! crust/mantle
        do j = 1,NGLOB_CRUST_MANTLE_ADJOINT
          do i = 1, NDIM
            b_displ_crust_mantle(i,j) = b_displ_cm_store_buffer(i,j,it_subset_end-it_of_this_subset+1)
          enddo
        enddo
        ! outer core
        do j = 1,NGLOB_OUTER_CORE_ADJOINT
            b_displ_outer_core(j) = b_displ_oc_store_buffer(j,it_subset_end-it_of_this_subset+1)
            b_accel_outer_core(j) = b_accel_oc_store_buffer(j,it_subset_end-it_of_this_subset+1)
        enddo
        ! inner core
        do j = 1,NGLOB_INNER_CORE_ADJOINT
          do i = 1, NDIM
            b_displ_inner_core(i,j) = b_displ_ic_store_buffer(i,j,it_subset_end-it_of_this_subset+1)
          enddo
        enddo

        ! for noise kernel
        if (NOISE_TOMOGRAPHY == 3) then
          noise_surface_movie(:,:,:,:) = b_noise_surface_movie_buffer(:,:,:,:,it_subset_end-it_of_this_subset+1)
        endif

        ! transfers wavefields from CPU to GPU
        if (GPU_MODE) then
          ! daniel debug: check if these transfers could be made async to overlap
          call transfer_b_displ_cm_to_device(NDIM*NGLOB_CRUST_MANTLE,b_displ_crust_mantle,Mesh_pointer)
          call transfer_b_displ_ic_to_device(NDIM*NGLOB_INNER_CORE,b_displ_inner_core,Mesh_pointer)
          call transfer_b_displ_oc_to_device(NGLOB_OUTER_CORE,b_displ_outer_core,Mesh_pointer)
          call transfer_b_accel_oc_to_device(NGLOB_OUTER_CORE,b_accel_outer_core,Mesh_pointer)
        endif

        it = it + 1

        ! simulation status output and stability check
        if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
          call check_stability()
        endif

        ! computes adjoint wavefield
        do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then

          if (USE_LDDRK) then
            ! update displacement using Runge-Kutta time scheme
            call update_displ_lddrk()
          else
            ! update displacement using Newmark time scheme
            call update_displ_Newmark()
          endif

          ! acoustic solver for outer core
          ! (needs to be done first, before elastic one)
          call compute_forces_acoustic()

          ! elastic solver for crust/mantle and inner core
          call compute_forces_viscoelastic()

        enddo ! istage

        ! write the seismograms with time shift
        call write_seismograms()

        ! kernel computation
        ! adjoint simulations: kernels
        call compute_kernels()

      enddo ! subset loop

    end select ! SIMULATION_TYPE

  !
  !---- end of time iteration loop
  !
  enddo   ! end of main time loop

  ! frees undo_attenuation buffers
  if (SIMULATION_TYPE == 3) then
    deallocate(b_displ_cm_store_buffer, &
               b_displ_oc_store_buffer, &
               b_accel_oc_store_buffer, &
               b_displ_ic_store_buffer)
    ! noise simulations
    if (NOISE_TOMOGRAPHY == 3) then
      deallocate(b_noise_surface_movie_buffer)
    endif
  endif

  ! close the huge file that contains a dump of all the time steps to disk
  if (EXACT_UNDOING_TO_DISK) call finish_exact_undoing_to_disk()

  ! user output of runtime
  call print_elapsed_time()

  ! Transfer fields from GPU card to host for further analysis
  if (GPU_MODE) call it_transfer_from_GPU()

!----  close energy file
  if (OUTPUT_ENERGY .and. myrank == 0) close(IOUT_ENERGY)

  ! safety check of last time loop increment
  if (it /= it_end) then
    print *,'Error time increments: it_end = ',it_end,' and last it = ',it,' do not match!'
    call exit_MPI(myrank,'Error invalid time increment ending')
  endif

  end subroutine iterate_time_undoatt

!
!-------------------------------------------------------------------------------------------------
!
! compute the strain in the whole crust/mantle and inner core domains
!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_compute_strain_att()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none
  ! local parameters
  integer :: ispec

  ! computes strain based on forward wavefield displ
  if (.not. GPU_MODE) then

    ! checks
    if (USE_DEVILLE_PRODUCTS_VAL) then

      ! inner core
      do ispec = 1, NSPEC_INNER_CORE
        call compute_element_strain_att_Dev(ispec,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
                                            displ_inner_core,veloc_inner_core,0._CUSTOM_REAL, &
                                            ibool_inner_core, &
                                            hprime_xx,hprime_xxT, &
                                            xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                            etax_inner_core,etay_inner_core,etaz_inner_core, &
                                            gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                            epsilondev_xx_inner_core, &
                                            epsilondev_yy_inner_core, &
                                            epsilondev_xy_inner_core, &
                                            epsilondev_xz_inner_core, &
                                            epsilondev_yz_inner_core, &
                                            NSPEC_INNER_CORE_STRAIN_ONLY,eps_trace_over_3_inner_core)
      enddo
      ! crust mantle
      do ispec = 1, NSPEC_crust_mantle
        call compute_element_strain_att_Dev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                                            displ_crust_mantle,veloc_crust_mantle,0._CUSTOM_REAL, &
                                            ibool_crust_mantle, &
                                            hprime_xx,hprime_xxT, &
                                            xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                            etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                            gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                            epsilondev_xx_crust_mantle, &
                                            epsilondev_yy_crust_mantle, &
                                            epsilondev_xy_crust_mantle, &
                                            epsilondev_xz_crust_mantle, &
                                            epsilondev_yz_crust_mantle, &
                                            NSPEC_CRUST_MANTLE_STRAIN_ONLY,eps_trace_over_3_crust_mantle)
      enddo

    else

      ! inner core
      do ispec = 1, NSPEC_INNER_CORE
        call compute_element_strain_att_noDev(ispec,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
                                              displ_inner_core,veloc_inner_core,0._CUSTOM_REAL, &
                                              ibool_inner_core, &
                                              hprime_xx,hprime_yy,hprime_zz, &
                                              xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                              etax_inner_core,etay_inner_core,etaz_inner_core, &
                                              gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                              epsilondev_xx_inner_core, &
                                              epsilondev_yy_inner_core, &
                                              epsilondev_xy_inner_core, &
                                              epsilondev_xz_inner_core, &
                                              epsilondev_yz_inner_core, &
                                              NSPEC_INNER_CORE_STRAIN_ONLY,eps_trace_over_3_inner_core)
      enddo
      ! crust mantle
      do ispec = 1, NSPEC_CRUST_MANTLE
        call compute_element_strain_att_noDev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                                              displ_crust_mantle,veloc_crust_mantle,0._CUSTOM_REAL, &
                                              ibool_crust_mantle, &
                                              hprime_xx,hprime_yy,hprime_zz, &
                                              xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                              etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                              gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                              epsilondev_xx_crust_mantle, &
                                              epsilondev_yy_crust_mantle, &
                                              epsilondev_xy_crust_mantle, &
                                              epsilondev_xz_crust_mantle, &
                                              epsilondev_yz_crust_mantle, &
                                              NSPEC_CRUST_MANTLE_STRAIN_ONLY,eps_trace_over_3_crust_mantle)
      enddo
    endif

  else

    ! calculates strains on GPU
    ! note: deltat is zero, thus strain is computed based on < displ(:,:) > rather than < displ(:,:) + deltat * veloc(:,:) >
    !       nevertheless, we implement < displ(:,:) + deltat * veloc(:,:) > in order to have a more general calculation
    !       as done in the CPU routine as well
    call compute_strain_gpu(Mesh_pointer,0._CUSTOM_REAL,1)

  endif ! GPU_MODE

  end subroutine it_compute_strain_att

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_compute_strain_att_backward()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none

  ! local parameters
  integer :: ispec

  ! computes strain based on backward/reconstructed wavefield b_displ
  if (.not. GPU_MODE) then

    ! checks
    if (USE_DEVILLE_PRODUCTS_VAL) then

      ! inner core
      do ispec = 1, NSPEC_INNER_CORE
        call compute_element_strain_att_Dev(ispec,NGLOB_INNER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                                            b_displ_inner_core,b_veloc_inner_core,0._CUSTOM_REAL, &
                                            ibool_inner_core, &
                                            hprime_xx,hprime_xxT, &
                                            xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                            etax_inner_core,etay_inner_core,etaz_inner_core, &
                                            gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                            b_epsilondev_xx_inner_core, &
                                            b_epsilondev_yy_inner_core, &
                                            b_epsilondev_xy_inner_core, &
                                            b_epsilondev_xz_inner_core, &
                                            b_epsilondev_yz_inner_core, &
                                            NSPEC_INNER_CORE_STRAIN_ONLY,b_eps_trace_over_3_inner_core)
      enddo

      ! crust mantle
      do ispec = 1, NSPEC_CRUST_MANTLE
        call compute_element_strain_att_Dev(ispec,NGLOB_CRUST_MANTLE_ADJOINT,NSPEC_CRUST_MANTLE_ADJOINT, &
                                            b_displ_crust_mantle,b_veloc_crust_mantle,0._CUSTOM_REAL, &
                                            ibool_crust_mantle, &
                                            hprime_xx,hprime_xxT, &
                                            xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                            etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                            gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                            b_epsilondev_xx_crust_mantle, &
                                            b_epsilondev_yy_crust_mantle, &
                                            b_epsilondev_xy_crust_mantle, &
                                            b_epsilondev_xz_crust_mantle, &
                                            b_epsilondev_yz_crust_mantle, &
                                            NSPEC_CRUST_MANTLE_STRAIN_ONLY,b_eps_trace_over_3_crust_mantle)
      enddo

    else

      ! inner core
      do ispec = 1, NSPEC_INNER_CORE
        call compute_element_strain_att_noDev(ispec,NGLOB_INNER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                                              b_displ_inner_core,b_veloc_inner_core,0._CUSTOM_REAL, &
                                              ibool_inner_core, &
                                              hprime_xx,hprime_yy,hprime_zz, &
                                              xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                              etax_inner_core,etay_inner_core,etaz_inner_core, &
                                              gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                              b_epsilondev_xx_inner_core, &
                                              b_epsilondev_yy_inner_core, &
                                              b_epsilondev_xy_inner_core, &
                                              b_epsilondev_xz_inner_core, &
                                              b_epsilondev_yz_inner_core, &
                                              NSPEC_INNER_CORE_STRAIN_ONLY,b_eps_trace_over_3_inner_core)
      enddo
      ! crust mantle
      do ispec = 1, NSPEC_crust_mantle
        call compute_element_strain_att_noDev(ispec,NGLOB_CRUST_MANTLE_ADJOINT,NSPEC_CRUST_MANTLE_ADJOINT, &
                                              b_displ_crust_mantle,b_veloc_crust_mantle,0._CUSTOM_REAL, &
                                              ibool_crust_mantle, &
                                              hprime_xx,hprime_yy,hprime_zz, &
                                              xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                              etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                              gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                              b_epsilondev_xx_crust_mantle, &
                                              b_epsilondev_yy_crust_mantle, &
                                              b_epsilondev_xy_crust_mantle, &
                                              b_epsilondev_xz_crust_mantle, &
                                              b_epsilondev_yz_crust_mantle, &
                                              NSPEC_CRUST_MANTLE_STRAIN_ONLY,b_eps_trace_over_3_crust_mantle)
      enddo
    endif

  else

    ! calculates strains on GPU
    call compute_strain_gpu(Mesh_pointer,0._CUSTOM_REAL,3)

  endif ! GPU_MODE

  end subroutine it_compute_strain_att_backward

