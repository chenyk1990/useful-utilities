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

  program xspecfem3D

  use specfem_par

  implicit none

!=======================================================================!
!                                                                       !
!   specfem3D is a 3-D spectral-element solver for the Earth.           !
!   It uses a mesh generated by program meshfem3D                       !
!                                                                       !
!=======================================================================!
!
! If you use this code for your own research, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {Communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! @ARTICLE{PeKoLuMaLeCaLeMaLiBlNiBaTr11,
! author = {Daniel Peter and Dimitri Komatitsch and Yang Luo and Roland Martin
!     and Nicolas {Le Goff} and Emanuele Casarotti and Pieyre {Le Loher}
!     and Federica Magnoni and Qinya Liu and C\'eline Blitz and Tarje Nissen-Meyer
!     and Piero Basini and Jeroen Tromp},
! title = {Forward and adjoint simulations of seismic wave propagation on fully
!     unstructured hexahedral meshes},
! journal={Geophys. J. Int.},
! year = {2011},
! volume = {186},
! pages = {721-739},
! number = {2},
! doi = {10.1111/j.1365-246X.2011.05044.x}}
!
! or
!
! @ARTICLE{VaCaSaKoVi99,
! author = {R. Vai and J. M. Castillo-Covarrubias and F. J. S\'anchez-Sesma and
! D. Komatitsch and J. P. Vilotte},
! title = {Elastic wave propagation in an irregularly layered medium},
! journal = {Soil Dynamics and Earthquake Engineering},
! year = {1999},
! volume = {18},
! pages = {11-18},
! number = {1},
! doi = {10.1016/S0267-7261(98)00027-X}}
!
! @ARTICLE{LeChKoHuTr09,
! author = {Shiann Jong Lee and Yu Chang Chan and Dimitri Komatitsch and Bor
! Shouh Huang and Jeroen Tromp},
! title = {Effects of realistic surface topography on seismic ground motion
! in the {Y}angminshan region of {T}aiwan based upon the spectral-element
! method and {LiDAR DTM}},
! journal = {Bull. Seismol. Soc. Am.},
! year = {2009},
! volume = {99},
! pages = {681-693},
! number = {2A},
! doi = {10.1785/0120080264}}
!
! @ARTICLE{LeChLiKoHuTr08,
! author = {Shiann Jong Lee and How Wei Chen and Qinya Liu and Dimitri Komatitsch
! and Bor Shouh Huang and Jeroen Tromp},
! title = {Three-Dimensional Simulations of Seismic Wave Propagation in the
! {T}aipei Basin with Realistic Topography Based upon the Spectral-Element Method},
! journal = {Bull. Seismol. Soc. Am.},
! year = {2008},
! volume = {98},
! pages = {253-264},
! number = {1},
! doi = {10.1785/0120070033}}
!
! @ARTICLE{LeKoHuTr09,
! author = {S. J. Lee and Dimitri Komatitsch and B. S. Huang and J. Tromp},
! title = {Effects of topography on seismic wave propagation: An example from
! northern {T}aiwan},
! journal = {Bull. Seismol. Soc. Am.},
! year = {2009},
! volume = {99},
! pages = {314-325},
! number = {1},
! doi = {10.1785/0120080020}}
!
! @ARTICLE{KoErGoMi10,
! author = {Dimitri Komatitsch and Gordon Erlebacher and Dominik G\"oddeke and
! David Mich\'ea},
! title = {High-order finite-element seismic wave propagation modeling with
! {MPI} on a large {GPU} cluster},
! journal = {J. Comput. Phys.},
! year = {2010},
! volume = {229},
! pages = {7692-7714},
! number = {20},
! doi = {10.1016/j.jcp.2010.06.024}}
!
! @ARTICLE{KoGoErMi10,
! author = {Dimitri Komatitsch and Dominik G\"oddeke and Gordon Erlebacher and
! David Mich\'ea},
! title = {Modeling the propagation of elastic waves using spectral elements
! on a cluster of 192 {GPU}s},
! journal = {Computer Science Research and Development},
! year = {2010},
! volume = {25},
! pages = {75-82},
! number = {1-2},
! doi = {10.1007/s00450-010-0109-1}}
!
! @ARTICLE{KoMiEr09,
! author = {Dimitri Komatitsch and David Mich\'ea and Gordon Erlebacher},
! title = {Porting a high-order finite-element earthquake modeling application
! to {NVIDIA} graphics cards using {CUDA}},
! journal = {Journal of Parallel and Distributed Computing},
! year = {2009},
! volume = {69},
! pages = {451-460},
! number = {5},
! doi = {10.1016/j.jpdc.2009.01.006}}
!
! @INCOLLECTION{ChKoViCaVaFe07,
! author = {Emmanuel Chaljub and Dimitri Komatitsch and Jean-Pierre Vilotte and
! Yann Capdeville and Bernard Valette and Gaetano Festa},
! title = {Spectral Element Analysis in Seismology},
! booktitle = {Advances in Wave Propagation in Heterogeneous Media},
! publisher = {Elsevier - Academic Press},
! year = {2007},
! editor = {Ru-Shan Wu and Val\'erie Maupin},
! volume = {48},
! series = {Advances in Geophysics},
! pages = {365-419}}
!
! @ARTICLE{KoVi98,
! author={D. Komatitsch and J. P. Vilotte},
! title={The spectral-element method: an efficient tool to simulate the seismic response of 2{D} and 3{D} geological structures},
! journal={Bull. Seismol. Soc. Am.},
! year=1998,
! volume=88,
! number=2,
! pages={368-392}}
!
! @ARTICLE{KoTr99,
! author={D. Komatitsch and J. Tromp},
! year=1999,
! title={Introduction to the spectral-element method for 3-{D} seismic wave propagation},
! journal={Geophys. J. Int.},
! volume=139,
! number=3,
! pages={806-822},
! doi={10.1046/j.1365-246x.1999.00967.x}}
!
! @ARTICLE{KoRiTr02,
! author={D. Komatitsch and J. Ritsema and J. Tromp},
! year=2002,
! title={The Spectral-Element Method, {B}eowulf Computing, and Global Seismology},
! journal={Science},
! volume=298,
! number=5599,
! pages={1737-1742},
! doi={10.1126/science.1076024}}
!
! @ARTICLE{KoTr02a,
! author={D. Komatitsch and J. Tromp},
! year=2002,
! title={Spectral-Element Simulations of Global Seismic Wave Propagation{-I. V}alidation},
! journal={Geophys. J. Int.},
! volume=149,
! number=2,
! pages={390-412},
! doi={10.1046/j.1365-246X.2002.01653.x}}
!
! @ARTICLE{KoTr02b,
! author={D. Komatitsch and J. Tromp},
! year=2002,
! title={Spectral-Element Simulations of Global Seismic Wave Propagation{-II. 3-D} Models, Oceans, Rotation, and Self-Gravitation},
! journal={Geophys. J. Int.},
! volume=150,
! pages={303-318},
! number = 1,
! doi={10.1046/j.1365-246X.2002.01716.x}}
!
! and/or another article from http://web.univ-pau.fr/~dkomati1/publications.html
!
!
! If you use the kernel capabilities of the code, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {Communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! @ARTICLE{PeKoLuMaLeCaLeMaLiBlNiBaTr11,
! author = {Daniel Peter and Dimitri Komatitsch and Yang Luo and Roland Martin
!     and Nicolas {Le Goff} and Emanuele Casarotti and Pieyre {Le Loher}
!     and Federica Magnoni and Qinya Liu and C\'eline Blitz and Tarje Nissen-Meyer
!     and Piero Basini and Jeroen Tromp},
! title = {Forward and adjoint simulations of seismic wave propagation on fully
!     unstructured hexahedral meshes},
! journal={Geophys. J. Int.},
! year = {2011},
! volume = {186},
! pages = {721-739},
! number = {2},
! doi = {10.1111/j.1365-246X.2011.05044.x}}
!
! @ARTICLE{LiTr06,
! author={Qinya Liu and Jeroen Tromp},
! title={Finite-frequency kernels based on adjoint methods},
! journal={Bull. Seismol. Soc. Am.},
! year=2006,
! volume=96,
! number=6,
! pages={2383-2397},
! doi={10.1785/0120060041}}
!
! If you use 3-D model S20RTS, please cite:
!
! @ARTICLE{RiVa00,
! author={J. Ritsema and H. J. {Van Heijst}},
! year=2000,
! title={Seismic imaging of structural heterogeneity in {E}arth's mantle: Evidence for large-scale mantle flow},
! journal={Science Progress},
! volume=83,
! pages={243-259}}
!
! Reference frame - convention:
! ----------------------------
!
! The code uses the following convention for the reference frame:
!
!  - X axis is East
!  - Y axis is North
!  - Z axis is up
!
! Note that this convention is different from both the Aki-Richards convention
! and the Harvard CMT convention.
!
! Let us recall that the Aki-Richards convention is:
!
!  - X axis is North
!  - Y axis is East
!  - Z axis is down
!
! and that the Harvard CMT convention is:
!
!  - X axis is South
!  - Y axis is East
!  - Z axis is up
!
! To report bugs or suggest improvements to the code, please send an email
! to Jeroen Tromp <jtromp AT princeton.edu> and/or use our online
! bug tracking system at http://www.geodynamics.org/roundup .
!
! Evolution of the code:
! ---------------------
!
! v. 7.0, many developers, January 2015:
!     simultaneous MPI runs, ADIOS file I/O support, ASDF seismograms, new seismogram names, tomography tools,
!     CUDA and OpenCL GPU support, CEM model support, updates AK135 model, binary topography files,
!     fixes geocentric/geographic conversions, updates ellipticity and gravity factors, git versioning system.
!
! v. 6.0, Daniel Peter (ETH Z\"urich, Switzerland), Dimitri Komatitsch and Zhinan Xie (CNRS / University of Marseille, France),
!     Elliott Sales de Andrade (University of Toronto, Canada), and many others, in particular from Princeton University, USA,
!     April 2014:
!     more flexible MPI implementation, GPU support, exact undoing of attenuation, LDDRK4-6 higher-order time scheme, etc...
!
! v. 5.1, Dimitri Komatitsch, University of Toulouse, France and Ebru Bozdag, Princeton University, USA, February 2011:
!     non blocking MPI for much better scaling on large clusters;
!     new convention for the name of seismograms, to conform to the IRIS standard;
!     new directory structure
!
! v. 5.0, many developers, February 2010:
!     new moho mesh stretching honoring crust2.0 moho depths,
!     new attenuation assignment, new SAC headers, new general crustal models,
!     faster performance due to Deville routines and enhanced loop unrolling,
!     slight changes in code structure (see also trivia at program start)
!
! v. 4.0 David Michea and Dimitri Komatitsch, University of Pau, France, February 2008:
!      first port to GPUs using CUDA, new doubling brick in the mesh, new perfectly load-balanced mesh,
!      more flexible routines for mesh design, new inflated central cube
!      with optimized shape, far fewer mesh files saved by the mesher,
!      global arrays sorted to speed up the simulation, seismos can be
!      written by the master, one more doubling level at the bottom
!      of the outer core if needed (off by default)
!
! v. 3.6 Many people, many affiliations, September 2006:
!      adjoint and kernel calculations, fixed IASP91 model,
!      added AK135 and 1066a, fixed topography/bathymetry routine,
!      new attenuation routines, faster and better I/Os on very large
!      systems, many small improvements and bug fixes, new "configure"
!      script, new user's manual etc.
!
! v. 3.5 Dimitri Komatitsch, Brian Savage and Jeroen Tromp, Caltech, July 2004:
!      any size of chunk, 3D attenuation, case of two chunks,
!      more precise topography/bathymetry model, new Par_file structure
!
! v. 3.4 Dimitri Komatitsch and Jeroen Tromp, Caltech, August 2003:
!      merged global and regional codes, no iterations in fluid, better movies
!
! v. 3.3 Dimitri Komatitsch, Caltech, September 2002:
!      flexible mesh doubling in outer core, inlined code, OpenDX support
!
! v. 3.2 Jeroen Tromp, Caltech, July 2002:
!      multiple sources and flexible PREM reading
!
! v. 3.1 Dimitri Komatitsch, Caltech, June 2002:
!      vectorized loops in solver and merged central cube
!
! v. 3.0 Dimitri Komatitsch and Jeroen Tromp, Caltech, May 2002:
!   ported to SGI and Compaq, double precision solver, more general anisotropy
!
! v. 2.3 Dimitri Komatitsch and Jeroen Tromp, Caltech, August 2001:
!                       gravity, rotation, oceans and 3-D models
!
! v. 2.2 Dimitri Komatitsch and Jeroen Tromp, Caltech, USA, March 2001:
!                       final MPI package
!
! v. 2.0 Dimitri Komatitsch, Harvard, USA, January 2000: MPI code for the globe
!
! v. 1.0 Dimitri Komatitsch, UNAM, Mexico, June 1999: first MPI code for a chunk
!
! Jeroen Tromp and Dimitri Komatitsch, Harvard, USA, July 1998: first chunk solver using OpenMP on a Sun machine
!
! Dimitri Komatitsch, IPG Paris, France, December 1996: first 3-D solver for the CM-5 Connection Machine,
!    parallelized on 128 processors using Connection Machine Fortran
!
! From Dahlen and Tromp (1998):
! ----------------------------
!
! Gravity is approximated by solving eq (3.259) without the Phi_E' term
! The ellipsoidal reference model is that of section 14.1
! The transversely isotropic expression for PREM is that of eq (8.190)
!
! Formulation in the fluid (acoustic) outer core:
! -----------------------------------------------
!
! In case of an acoustic medium, a displacement potential Chi is used
! as in Chaljub and Valette, Geophysical Journal International, vol. 158,
! p. 131-141 (2004) and *NOT* a velocity potential as in Komatitsch and Tromp,
! Geophysical Journal International, vol. 150, p. 303-318 (2002).
! This permits acoustic-elastic coupling based on a non-iterative time scheme.
! Displacement if we ignore gravity is then: u = grad(Chi)
! (In the context of the Cowling approximation displacement is
! u = grad(rho * Chi) / rho, *not* u = grad(Chi).)
! Velocity is then: v = grad(Chi_dot)       (Chi_dot being the time derivative of Chi)
! and pressure is: p = - rho * Chi_dot_dot  (Chi_dot_dot being the time second derivative of Chi).
! The source in an acoustic element is a pressure source.
! The potential in the outer core is called displ_outer_core for simplicity.
! Its first time derivative is called veloc_outer_core.
! Its second time derivative is called accel_outer_core.

! *************************************************
! ************** PROGRAM STARTS HERE **************
! *************************************************
!
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!
! trivia about the programming style adopted here:
!
! note 1: for performance reasons, we try to use as much from the stack memory as possible.
!             This is done to avoid memory fragmentation and also to optimize performance.
!             Stack memory is a place in computer memory where all the variables that are declared
!             and initialized **before** runtime are stored. Our static array allocation will use that one.
!             All variables declared within our main routine also will be stored on the stack.
!
!             the heap is the section of computer memory where all the variables created or initialized
!             **at** runtime are stored. it is used for dynamic memory allocation.
!
!             stack is much faster than the heap.
!
!             when calling a function, additional storage will be allocated for the variables in that function.
!             that storage will be allocated in the heap memory segment.
!
!             most routine calls here will have rather long argument lists, probably because of this performance criteria.
!             using modules/common data blocks together with dynamic allocation will put data into heap memory,
!             thus it has longer latency to access variables than stack memory variables.
!
!             however, declaring the static arrays needed in compute_forces_crust_mantle_Dev()
!             like e.g. sum_terms, tempx1,... in this main routine and
!             passing them along as arguments to the routine makes the code slower.
!             it seems that this stack/heap criterion is more complicated,
!             and inlining functions is a performance criteria as well.
!
!             for vectorization, we assume that arrays have contiguous memory allocated. this holds true for most compilers and
!             static memory allocation. however, note that dynamically allocated memory could in principle be non-contiguous.
!
!             another reason why the use of modules is restricted is to make the code thread safe.
!             having different threads access the same data structure and modifying it at the same time
!             would lead to problems. passing arguments is a way to avoid such complications.
!
! note 2: Most of the computation time is spent
!             inside the time loop (mainly in the compute_forces_crust_mantle_Dev() routine).
!             Any code performance tuning will be most effective in there.
!
! note 3: Fortran is a code language that uses column-first ordering for arrays,
!             e.g., it stores a(i,j) in this order: a(1,1),a(2,1),a(3,1),...,a(1,2),a(2,2),a(3,2),..
!             it is therefore more efficient to have the inner-loop over i, and the outer loop over j
!
! note 4: Deville et al. (2002) routines significantly reduce the total number of memory accesses
!             required to perform matrix-matrix products at the spectral element level.
!             For most compilers and hardware, will result in a significant speedup (> 30% or more, sometimes twice faster).
!
! note 5: whenever adding some new code, please make sure to use
!             spaces rather than tabs. Tabulators are in principle not allowed in Fortran95.
!
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!
  ! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call init_mpi()

  ! force Flush-To-Zero if available to avoid very slow Gradual Underflow trapping
#ifndef __MIC__
  call force_ftz()
#endif

  ! initializes simulation parameters
  call initialize_simulation()

  ! sets up reference element GLL points/weights/derivatives
  call setup_GLL_points()

  ! starts reading the databases
  call read_mesh_databases()

  ! reads topography & bathymetry & ellipticity
  call read_topography_bathymetry()

  ! prepares sources and receivers
  call setup_sources_receivers()

  ! sets up and precomputes simulation arrays
  call prepare_timerun()

  ! steps through time iterations
  if (UNDO_ATTENUATION) then
    call iterate_time_undoatt()
  else
    call iterate_time()
  endif

  ! saves last time frame and finishes kernel calculations
  call finalize_simulation()

  ! stop all the MPI processes, and exit
  call finalize_mpi()

  end program xspecfem3D

