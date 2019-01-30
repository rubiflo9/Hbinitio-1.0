module global
  implicit none
  ! ------------------------------------------
  type t_potential
     double precision,allocatable :: ext(:) ! external potential
     double precision,allocatable :: hartree(:) ! hartreel potential
     double precision,allocatable :: Vx(:) ! exchange potential
     double precision,allocatable :: perturb(:) ! perturbation potential
     double precision,allocatable :: tot(:) ! perturbation potential
     double precision::EX,Ehartree
  end type t_potential
  ! ------------------------------------------
  type t_box
     character(len=32)::shape
     double precision::width
     double precision::center(3)
     double precision::radius
  end type t_box
  ! --------------------------------------------------------
  !
  !     PERTURBATION data type
  !
  ! -------------------------------------------------------
  type t_perturbation
     character(len=32)::shape
     double precision::Intensity
     double precision::sigma
     double precision::location(3)
  end type t_perturbation
  ! --------------------------------------------------------
  !
  !     PARAM data type
  !
  ! -------------------------------------------------------
  type t_param
     logical::restart
     character(len=32)::scheme    ! numerov || davidson
     character(len=32)::prefix
     character (len=1024)::filenameeigen
     character (len=1024)::filenrj
     character (len=1024)::inputfile
     logical::init_wf
     logical::extrapol
     integer::extrap_add
     integer::ieof
     integer::loopmax
     integer::nvecmin
     integer::nvecmax
     integer::Nx
     ! eigenvectors to converge
     integer::nvec_to_cvg
     integer,allocatable::list_idx_to_cvg(:)
     integer::noccstate
     double precision,allocatable::occupation(:)
     ! accurency
     double precision :: ETA
     type(t_box)::box
     type(t_perturbation)::perturb
     integer:: dim !dimension of the mesh 1(1D), 2(2D) or 3(3D)
     logical:: hartree
     logical::exchange
     double precision::Z
  end type t_param
  !------------------------------------------
  type t_nrj
     double precision::last
     double precision::previous
     double precision::dnrj
  end type t_nrj
  type tt_cvg
     double precision::nrj
     double precision::nrjold
     double precision::dnrj
     double precision:: resi
     logical::cvg
  end type tt_cvg
  type t_cvg
     integer::nwfc
     type(t_nrj)::total_nrj
     type(tt_cvg),allocatable::wfc(:)
     integer,allocatable::list_idx(:)
     double precision :: ETA
     integer :: nvec_to_cvg
     integer :: ncvg
     logical::cvg
  end type t_cvg
    !------------------------------------------
  type t_point
     double precision::q(3) ! coordinate
     double precision::d     ! distance from the center of the cell
     double precision::val   ! value of the Hartree potential
  end type t_point
  !------------------------------------------
  type t_ijk_to_idx
     integer::n
     double precision::q(3)
     logical::active
  end type t_ijk_to_idx
  !------------------------------------------
  type t_node
     integer::n_neighbors
     integer,allocatable::list_neighbors(:)
     integer::n_bound
     integer,allocatable::list_bound(:)
     logical::usefull_unactive
     logical::active
     integer::i,j,k
     double precision::q(3)
  end type t_node
  ! -----------------------------------------------------------
  !
  !   MESH data type
  !
  ! ----------------------------------------------------------
  type t_mesh
     integer :: Nx,Ny,Nz,Ntot,nactive,nunactive
     integer,allocatable :: list_bound(:,:),n_bound(:)  ! list_bound is linked with bound(:) 
     double precision :: dx,dy,dz,dv
     type(t_box)::box
     type(t_perturbation)::perturb
     integer :: dim
     integer::nbound
     integer::n_usefull_unactive
     type(t_point),allocatable::bound(:)
     type(t_ijk_to_idx),allocatable::ijk_to_idx(:,:,:)  ! from (i,j,k) -> n
     type(t_node),allocatable::node(:)
  end type t_mesh
  !------------------------------------------
  type t_wavefunction
     integer :: nwfc
     double precision,allocatable::wfc(:,:)
     double precision,allocatable::eps(:)
     double precision,allocatable ::epsprev(:),deps(:) ! eigenvalues
  end type t_wavefunction
  !------------------------------------------
  type t_molecule
     type(t_wavefunction)::wf       ! wavefunctions of the molecule
     type(t_mesh)::mesh               !
     type(t_potential)::pot
     double precision,allocatable::rho(:)
  end type t_molecule

contains
end module global
