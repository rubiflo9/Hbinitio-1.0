  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program Hbinitio
  !$ use OMP_LIB
  use time_tracking
  use global
  use poten
  use IO
  use param_mod
  use mesh_mod
  use davidson_mod
  use numerov_mod_dev
  implicit none
  !  include 'mpif.h'
  type(t_cvg) :: cvg
  type(t_time) :: time_spent
  type(t_param)::param,param2
  !------------------------------------------
  type t_perturb
     double precision,allocatable::coeff(:,:)
  end type t_perturb
  type(t_perturb)::perturb
  type(t_molecule):: molecule,molecule2
  !------------------------------------------
  integer :: i,j,k,n,l
  double precision::x,y,z
  !  integer::ierr,my_id,num_procs
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CHARACTER(len=32) :: arg
  
  i = 0
  DO
     CALL get_command_argument(i, arg)
     IF (LEN_TRIM(arg) == 0) EXIT
     if(i.eq.1) param%inputfile=TRIM(arg)
     i = i+1
  END DO
  
  call time_tracking_init(time_spent)
  !  call mpi_init(ierr )
  !  call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
  !  call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

  call read_param(param)
  
  call new_molecule(molecule,param)
    
  ! call read_pp(pp)
  !       filename='pot_ext.cube'
  !       call save_cube_3D(pot_ext,filename,m)

  cvg%nwfc=param%nvecmin
  allocate(cvg%wfc(cvg%nwfc))
  do i=1,cvg%nwfc
     cvg%wfc%cvg=.FALSE.
  end do
  cvg%nvec_to_cvg=param%nvec_to_cvg
  cvg%ETA=param%ETA
  allocate(cvg%list_idx(param%nvec_to_cvg))
  do i=1,param%nvec_to_cvg
     cvg%list_idx(i)=param%list_idx_to_cvg(i)
  end do
  
  allocate(perturb%coeff(cvg%nvec_to_cvg,cvg%nvec_to_cvg))          
    
  if(param%scheme.eq.'numerov') then
     print *,'Starting NUMEROV scheme'
     call numerov(molecule,cvg,param)
  end if
  
  if(param%scheme.eq.'davidson') then
     print *,'Starting DAVIDSON scheme'
     call davidson(param,molecule%mesh,cvg,molecule,molecule%pot,time_spent)
     
     ! if(param%extrapol) then
     !    call read_param(param2)
     !    param2%dim=param%dim
     !    param2%box_width=param%box_width
     !    param2%Nx=param%Nx+param%extrap_add
     !    param2%init_wf=.FALSE.
     !    call new_mesh(mesh2,param2)  
     !    call init_pot(mesh2,param2,pot2)
     !    call new_molecule(molecule2,param2,mesh2)
        
     !    do k=1,mesh2%Nz
     !       do i=1,mesh2%Nx
     !          do j=1,mesh2%Ny
     !             n=j+(i-1)*mesh2%Ny+(k-1)*mesh2%Ny*mesh2%Nx;             
     !             x=i*mesh2%dx
     !             y=j*mesh2%dY
     !             z=k*mesh2%dZ
     !             do l=1,param2%nvecmin
     !                molecule2%wf%wfc(n,l)=interpolate(x,y,z,mesh,molecule,l)
     !             end do
     !          end do
     !       end do
     !    end do
     !    !    call save_cube_3D(molecule2%wf%wfc(:,1),'essai.cube',mesh2)
     !    !   stop
     !    call davidson(param2,mesh2,cvg,molecule2,pot2,time_spent)
     ! endif
  end if

  !--------------------------------------------------------------------------
  !
  ! perturbation theory
  !
  !--------------------------------------------------------------------------
  ! call calc_coeff(param,pot,mesh,molecule,perturb)
  ! do i=1,cvg%nvec_to_cvg
  !    print *,perturb%coeff(i,i),perturb%coeff(i,i)+molecule%wf%eps(i)
  ! end do

  ! do i=1,cvg%nvec_to_cvg
  !    print *,(perturb%coeff(i,j),j=1,cvg%nvec_to_cvg)
  ! end do
  !--------------------------------------------------------------------------
  !
  ! end of Hbinitio
  !
  !--------------------------------------------------------------------------
  deallocate(molecule%wf%epsprev)
  deallocate(molecule%wf%deps)
  deallocate(molecule%wf%eps)
  deallocate(molecule%pot%ext)
  deallocate(molecule%pot%perturb)
  deallocate(molecule%pot%tot)
  deallocate(perturb%coeff)
!  call free_mesh(molecule%mesh)
  call cpu_time(time_spent%end)
  if (cvg%ncvg.ge.cvg%nvec_to_cvg) print *,'Main > Job DONE !'
  print '("Main > Total Time = ",e16.6," seconds.")',time_spent%end-time_spent%start
!  call mpi_finalize(ierr)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
  ! --------------------------------------------------------------------------------------
  !
  !             new_molecule()
  !
  ! --------------------------------------------------------------------------------------
  subroutine new_molecule(molecule,param)
    implicit none
    type(t_molecule)::molecule
    type(t_param)::param

    call new_mesh(molecule%mesh,param)
    call init_pot(molecule%mesh,molecule%pot)
    molecule%wf%nwfc=param%nvecmin
    allocate(molecule%wf%eps(molecule%wf%nwfc))
    allocate(molecule%wf%epsprev(molecule%wf%nwfc))
    allocate(molecule%wf%deps(molecule%wf%nwfc))
    allocate(molecule%wf%wfc(molecule%mesh%nactive,molecule%wf%nwfc))
    allocate(molecule%rho(molecule%mesh%nactive))
  end subroutine new_molecule
  ! --------------------------------------------------------------------------------------
  !
  !             V_from_wfc()
  !
  ! --------------------------------------------------------------------------------------
       subroutine V_from_wfc(param,molecule,mesh,V)
         double precision, external :: ddot
         integer :: lorb
         type(t_param)::param
         type(t_molecule)::molecule
         type(t_mesh)::mesh
         double precision::V(:,:),normloc
          do lorb=1,param%nvecmin
             call dcopy(mesh%nactive,molecule%wf%wfc(:,lorb),1,V(:,lorb),1) ! g->d
             normloc=ddot(mesh%nactive,V(:,lorb),1,V(:,lorb),1)
             normloc=1.0/sqrt(normloc)
             call dscal(mesh%nactive,normloc,V(:,lorb),1)
          end do
!    GS%nindep=1
!    call GramSchmidt(V,Vdump,nvec,m,GS)

        end subroutine V_from_wfc
  ! --------------------------------------------------------------------------------------
  !
  !              calc_coeff()
  !
  ! --------------------------------------------------------------------------------------
  subroutine calc_coeff(param,pot,mesh,molecule,perturb)
    implicit none
    type(t_perturb)::perturb
    type(t_molecule)::molecule
    type(t_mesh)::mesh
    type(t_potential)::pot
    type(t_param)::param
    integer::i,j
    do i=1,param%nvec_to_cvg
       do j=1,param%nvec_to_cvg
          perturb%coeff(i,j)=mesh%dv*&
               sum(molecule%wf%wfc(:,i)*pot%perturb*molecule%wf%wfc(:,j))
       end do
    end do
  end subroutine calc_coeff
  

  
  
end program Hbinitio
