subroutine unrestricted_linear_response_A_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,lambda, & 
                                                 e,ERI_aa,ERI_ab,ERI_bb,A_lr)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dRPA
  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: ERI_aa(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_ab(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_bb(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: delta_dRPA
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,a,b,ia,jb

! Output variables

  double precision,intent(out)  :: A_lr(nSt,nSt)

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

!-----------------------------------------------
! Build A matrix for spin-conserving transitions
!-----------------------------------------------

  if(ispin == 1) then 

    ! alpha-alpha block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1
 
            A_lr(ia,jb) = (e(a,1) - e(i,1))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                        + lambda*ERI_aa(i,b,a,j) - (1d0 - delta_dRPA)*lambda*ERI_aa(i,b,j,a)

          end  do
        end  do
      end  do
    end  do

    ! alpha-beta block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1
 
            A_lr(ia,nSa+jb) = lambda*ERI_ab(i,b,a,j) 

          end  do
        end  do
      end  do
    end  do

    ! beta-alpha block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1
 
            A_lr(nSa+ia,jb) = lambda*ERI_ab(b,i,j,a)

          end  do
        end  do
      end  do
    end  do

    ! beta-beta block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1
 
            A_lr(nSa+ia,nSa+jb) = (e(a,2) - e(i,2))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                                + lambda*ERI_bb(i,b,a,j) - (1d0 - delta_dRPA)*lambda*ERI_bb(i,b,j,a)

          end  do
        end  do
      end  do
    end  do

  end if


end subroutine unrestricted_linear_response_A_matrix
