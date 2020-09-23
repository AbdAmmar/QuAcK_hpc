subroutine URPAx(doACFDT,exchange_kernel,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS,ENuc,EUHF, & 
                 ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,e)

! Perform random phase approximation calculation with exchange (aka TDHF) in the unrestricted formalism

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  double precision,intent(in)   :: eta
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_abab(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin

  integer                       :: nS_aa,nS_bb,nS_sc
  double precision,allocatable  :: Omega_sc(:)
  double precision,allocatable  :: XpY_sc(:,:)
  double precision,allocatable  :: XmY_sc(:,:)

  integer                       :: nS_ab,nS_ba,nS_sf
  double precision,allocatable  :: Omega_sf(:)
  double precision,allocatable  :: XpY_sf(:,:)
  double precision,allocatable  :: XmY_sf(:,:)

  double precision              :: rho_sc,rho_sf
  double precision              :: EcRPAx(nspin)
  double precision              :: EcAC(nspin)

! Hello world

  write(*,*)
  write(*,*)'*********************************************************************'
  write(*,*)'| Unrestricted random phase approximation calculation with exchange |'
  write(*,*)'*********************************************************************'
  write(*,*)

! Initialization

  EcRPAx(:) = 0d0
  EcAC(:)   = 0d0

! Spin-conserved transitions

  if(spin_conserved) then 

    ispin = 1

    ! Memory allocation

    nS_aa = nS(1)
    nS_bb = nS(2)
    nS_sc = nS_aa + nS_bb

    allocate(Omega_sc(nS_sc),XpY_sc(nS_sc,nS_sc),XmY_sc(nS_sc,nS_sc))

    call unrestricted_linear_response(ispin,.false.,.false.,.false.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,1d0,e, & 
                                      ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,rho_sc,EcRPAx(ispin),Omega_sc,XpY_sc,XmY_sc)
    call print_excitation('URPAx  ',5,nS_sc,Omega_sc)
!   call print_transition_vectors(nBas,nC,nO,nV,nR,nS,Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))


  endif

! Spin-flip transitions

  if(spin_flip) then

    ispin = 2

    ! Memory allocation

    nS_ab = (nO(1) - nC(1))*(nV(2) - nR(2))
    nS_ba = (nO(2) - nC(2))*(nV(1) - nR(1))
    nS_sf = nS_ab + nS_ba

    allocate(Omega_sf(nS_sf),XpY_sf(nS_sf,nS_sf),XmY_sf(nS_sf,nS_sf))

    call unrestricted_linear_response(ispin,.false.,.false.,.false.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sf,1d0,e, &
                                      ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,rho_sf,EcRPAx(ispin),Omega_sf,XpY_sf,XmY_sf)
    call print_excitation('URPAx  ',6,nS_sf,Omega_sf)
!   call print_transition_vectors(nBas,nC,nO,nV,nR,nS,Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))


  endif

! if(exchange_kernel) then

!   EcRPAx(1) = 0.5d0*EcRPAx(1)
!   EcRPAx(2) = 1.5d0*EcRPAx(2)

! end if

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@URPAx correlation energy (spin-conserved) =',EcRPAx(1)
  write(*,'(2X,A50,F20.10)') 'Tr@URPAx correlation energy (spin-flip)      =',EcRPAx(2)
  write(*,'(2X,A50,F20.10)') 'Tr@URPAx correlation energy                  =',EcRPAx(1) + EcRPAx(2)
  write(*,'(2X,A50,F20.10)') 'Tr@URPAx total energy                        =',ENuc + EUHF + EcRPAx(1) + EcRPAx(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

! if(doACFDT) then

!   write(*,*) '-------------------------------------------------------'
!   write(*,*) 'Adiabatic connection version of RPAx correlation energy'
!   write(*,*) '-------------------------------------------------------'
!   write(*,*)

!   call ACFDT(exchange_kernel,.false.,.false.,.false.,.false.,.false.,singlet,triplet,eta, &
!              nBas,nC,nO,nV,nR,nS,ERI,e,e,Omega,XpY,XmY,rho,EcAC)

!   if(exchange_kernel) then

!     EcAC(1) = 0.5d0*EcAC(1)
!     EcAC(2) = 1.5d0*EcAC(2)

!   end if

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10)') 'AC@RPAx correlation energy (singlet) =',EcAC(1)
!   write(*,'(2X,A50,F20.10)') 'AC@RPAx correlation energy (triplet) =',EcAC(2)
!   write(*,'(2X,A50,F20.10)') 'AC@RPAx correlation energy           =',EcAC(1) + EcAC(2)
!   write(*,'(2X,A50,F20.10)') 'AC@RPAx total energy                 =',ENuc + EUHF + EcAC(1) + EcAC(2)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

! end if

end subroutine URPAx
