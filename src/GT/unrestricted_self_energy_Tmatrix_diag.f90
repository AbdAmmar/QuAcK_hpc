subroutine unrestricted_self_energy_Tmatrix_diag(ispin,eta,nBas,nC,nO,nV,nR,nH,nP,e,Omega1,rho1,Omega2,rho2,EcGM,SigT)

! Compute diagonal of the correlation part of the T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nH,ispin
  integer,intent(in)            :: nP
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: Omega1(nP)
  double precision,intent(in)   :: rho1(nBas,nBas,nP,nspin)
  double precision,intent(in)   :: Omega2(nH)
  double precision,intent(in)   :: rho2(nBas,nBas,nH,nspin)

! Local variables

  integer                       :: i,j,a,b,p,cd,kl
  double precision              :: eps

! Output variables

  double precision,intent(inout)  :: EcGM(nspin)
  double precision,intent(inout)  :: SigT(nBas,nspin)

!----------------------------------------------
! Occupied part of the T-matrix self-energy 
!----------------------------------------------

  if(ispin==1) then 

    do p=nC(1)+1,nBas-nR(1)
      do i=nC(1)+1,nO(1)
        do cd=1,nP
          eps = e(p,1) + e(i,1) - Omega1(cd)
          SigT(p,1) = SigT(p,1) + rho1(p,i,cd,1)**2*eps/(eps**2 + eta**2)
        enddo
      enddo
    enddo

  end if  
 
!beta part
  
  if(ispin==2) then 

    do p=nC(2)+1,nBas-nR(2)
      do i=nC(2)+1,nO(2)
        do cd=1,nP
          eps = e(p,2) + e(i,2) - Omega1(cd)
          SigT(p,2) = SigT(p,2) + rho1(p,i,cd,2)**2*eps/(eps**2 + eta**2)
        enddo
      enddo
    enddo

  end if 
  
!----------------------------------------------
! Virtual part of the T-matrix self-energy
!----------------------------------------------
  
  !alpha part 
  
  if(ispin==1) then

    do p=nC(1)+1,nBas-nR(1)
      do a=nO(1)+1,nBas-nR(1)
        do kl=1,nH
          eps = e(p,1) + e(a,1) - Omega2(kl)
          SigT(p,1) = SigT(p,1) + rho2(p,a,kl,1)**2*eps/(eps**2 + eta**2)
        enddo
      enddo
    enddo
  
  end if

  !alpha part 

  if(ispin==2) then

    do p=nC(2)+1,nBas-nR(2)
      do a=nO(2)+1,nBas-nR(2)
        do kl=1,nH
          eps = e(p,2) + e(a,2) - Omega2(kl)
          SigT(p,2) = SigT(p,2) + rho2(p,a,kl,2)**2*eps/(eps**2 + eta**2)
        enddo
      enddo
    enddo

  end if 

!----------------------------------------------
! Galitskii-Migdal correlation energy
!----------------------------------------------
  if(ispin==1) then

    do i=nC(1)+1,nO(1)
      do j=nC(1)+1,nO(1)
        do cd=1,nP
          eps = e(i,1) + e(j,1) - Omega1(cd)
          EcGM(1) = EcGM(1) + rho1(i,j,cd,1)*rho1(i,j,cd,1)*eps/(eps**2 + eta**2)
        enddo
      enddo
    enddo 

    do a=nO(1)+1,nBas-nR(1)
      do b=nO(1)+1,nBas-nR(1)
        do kl=1,nH
          eps = e(a,1) + e(b,1) - Omega2(kl)
          EcGM(1) = EcGM(1) - rho2(a,b,kl,1)*rho2(a,b,kl,1)*eps/(eps**2 + eta**2)
        enddo
      enddo
    enddo

  end if

  if(ispin==2) then

    do i=nC(2)+1,nO(2)
      do j=nC(2)+1,nO(2)
        do cd=1,nP
          eps = e(i,2) + e(j,2) - Omega1(cd)
          EcGM(2) = EcGM(2) + rho1(i,j,cd,2)*rho1(i,j,cd,2)*eps/(eps**2 + eta**2)
        enddo
      enddo
    enddo

    do a=nO(2)+1,nBas-nR(2)
      do b=nO(2)+1,nBas-nR(2)
        do kl=1,nH
          eps = e(a,2) + e(b,2) - Omega2(kl)
          EcGM(2) = EcGM(2) - rho2(a,b,kl,2)*rho2(a,b,kl,2)*eps/(eps**2 + eta**2)
        enddo
      enddo
    enddo

  end if

end subroutine unrestricted_self_energy_Tmatrix_diag
