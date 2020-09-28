subroutine UCIS(spin_conserved,spin_flip,nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,dipole_int,eHF)

! Perform configuration interaction single calculation`

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_abab(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart,nspin)

! Local variables

  logical                       :: dump_matrix = .false.
  logical                       :: dump_trans = .false.
  integer                       :: ispin
  double precision              :: lambda

  integer                       :: nS_aa,nS_bb,nS_sc
  double precision,allocatable  :: A_sc(:,:)
  double precision,allocatable  :: Omega_sc(:)

  integer                       :: nS_ab,nS_ba,nS_sf
  double precision,allocatable  :: A_sf(:,:)
  double precision,allocatable  :: Omega_sf(:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|      Configuration Interaction Singles       |'
  write(*,*)'************************************************'
  write(*,*)

! Adiabatic connection scaling

  lambda = 1d0

!----------------------------!
! Spin-conserved transitions !
!----------------------------!

  if(spin_conserved) then

    ispin = 1

    ! Memory allocation

    nS_aa = nS(1)
    nS_bb = nS(2)
    nS_sc = nS_aa + nS_bb

    allocate(A_sc(nS_sc,nS_sc),Omega_sc(nS_sc))

    call unrestricted_linear_response_A_matrix(ispin,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,lambda,eHF, & 
                                               ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,A_sc)
 
    if(dump_matrix) then
      print*,'CIS matrix (spin-conserved transitions)'
      call matout(nS_sc,nS_sc,A_sc)
      write(*,*)
    endif

    call diagonalize_matrix(nS_sc,A_sc,Omega_sc)
    call print_excitation('UCIS        ',5,nS_sc,Omega_sc)
 
    if(dump_trans) then
      print*,'Spin-conserved CIS transition vectors'
      call matout(nS_sc,nS_sc,A_sc)
      write(*,*)
    endif

    deallocate(A_sc,Omega_sc)

  endif

!-----------------------!
! Spin-flip transitions !
!-----------------------!

  if(spin_flip) then

    ispin = 2

    ! Memory allocation

    nS_ab = (nO(1) - nC(1))*(nV(2) - nR(2))
    nS_ba = (nO(2) - nC(2))*(nV(1) - nR(1))
    nS_sf = nS_ab + nS_ba

    allocate(A_sf(nS_sf,nS_sf),Omega_sf(nS_sf))

    call unrestricted_linear_response_A_matrix(ispin,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sf,lambda,eHF, & 
                                               ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,A_sf)
 
    if(dump_matrix) then
      print*,'CIS matrix (spin-conserved transitions)'
      call matout(nS_sf,nS_sf,A_sf)
      write(*,*)
    endif

    call diagonalize_matrix(nS_sf,A_sf,Omega_sf)
    call print_excitation('UCIS        ',6,nS_sf,Omega_sf)
 
    if(dump_trans) then
      print*,'Spin-flip CIS transition vectors'
      call matout(nS_sf,nS_sf,A_sf)
      write(*,*)
    endif

    deallocate(A_sf,Omega_sf)

  endif

end subroutine UCIS
