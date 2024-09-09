subroutine phACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,ERI,XpY,XmY,EcAC)

! Compute the correlation energy via the adiabatic connection formula

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: exchange_kernel
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: ia,jb,kc
  integer                       :: ia0,jb0,ia_tmp
  double precision              :: delta_spin
  double precision              :: delta_Kx
  double precision              :: delta_tmp
  double precision,allocatable  :: Ap(:,:)
  double precision,allocatable  :: Bp(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: Y(:,:)
  double precision,allocatable  :: tmp1(:,:), tmp2(:,:)
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: EcAC

! Singlet or triplet manifold?

  delta_spin = 0d0
  if(ispin == 1) delta_spin = +1d0
  if(ispin == 2) delta_spin = -1d0

! Exchange kernel

  delta_Kx = 0d0
  if(exchange_kernel) delta_Kx = 1d0
  
! Memory allocation

  allocate(Ap(nS,nS),Bp(nS,nS))

  ! Compute Aiajb = (ia|bj) and Biajb = (ia|jb)
!  ia = 0
!  do i=nC+1,nO
!    do a=nO+1,nBas-nR
!      ia = ia + 1
!      jb = 0
!      do j=nC+1,nO
!        do b=nO+1,nBas-nR
!          jb = jb + 1
!          Ap(ia,jb) = (1d0 + delta_spin)*ERI(i,b,a,j) - delta_Kx*ERI(i,b,j,a)
!          Bp(ia,jb) = (1d0 + delta_spin)*ERI(i,j,a,b) - delta_Kx*ERI(i,j,b,a)
!        enddo
!      enddo
!    enddo
!  enddo

  ia_tmp = nBas - nR - nO
  delta_tmp = 1d0 + delta_spin

  !$OMP PARALLEL                   &
  !$OMP DEFAULT(NONE)              &
  !$OMP PRIVATE(i,a,ia,j,b,jb,jb0) &
  !$OMP SHARED(nC,nO,nV,nR,nBas,ia_tmp,delta_tmp,delta_Kx,ERI,Ap,Bp)
  !$OMP DO COLLAPSE(2)
  do i = nC+1, nO
    do a = nO+1, nBas-nR
      ia = (i - nC - 1) * ia_tmp + a - nO
      do j = nC+1, nO
        jb0 = (j - nC - 1) * ia_tmp
        do b = nO+1, nBas-nR
          jb = jb0 + b - nO
          Ap(ia,jb) = delta_tmp * ERI(i,b,a,j) - delta_Kx * ERI(i,b,j,a)
          Bp(ia,jb) = delta_tmp * ERI(i,j,a,b) - delta_Kx * ERI(i,j,b,a)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

! Compute Tr(K x P_lambda)

  !allocate(X(nS,nS),Y(nS,nS))
  !X(:,:) = 0.5d0*(XpY(:,:) + XmY(:,:))
  !Y(:,:) = 0.5d0*(XpY(:,:) - XmY(:,:))
  !EcAC = trace_matrix(nS,matmul(X,matmul(Bp,transpose(Y))) + matmul(Y,matmul(Bp,transpose(X)))) &
  !     + trace_matrix(nS,matmul(X,matmul(Ap,transpose(X))) + matmul(Y,matmul(Ap,transpose(Y)))) &
  !     - trace_matrix(nS,Ap)
  !deallocate(X,Y)

  allocate(tmp1(nS,nS),tmp2(nS,nS))
  call dgemm("N","T",nS,nS,nS,1.d0,Ap+Bp,nS,XpY(1,1),nS,0.d0,tmp1(1,1),nS)
  call dgemm("N","N",nS,nS,nS,0.5d0,XpY(1,1),nS,tmp1(1,1),nS,0.d0,tmp2(1,1),nS)
  call dgemm("N","T",nS,nS,nS,1.d0,Ap-Bp,nS,XmY(1,1),nS,0.d0,tmp1(1,1),nS)
  call dgemm("N","N",nS,nS,nS,0.5d0,XmY(1,1),nS,tmp1(1,1),nS,1.d0,tmp2(1,1),nS)
  EcAC = trace_matrix(nS,tmp2(1,1)) - trace_matrix(nS,Ap(1,1))
  deallocate(tmp1,tmp2)


  deallocate(Ap,Bp)

end subroutine 
