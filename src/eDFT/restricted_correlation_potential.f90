subroutine restricted_correlation_potential(rung,DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

! Compute the correlation potential 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Local variables

  double precision,allocatable  :: FcLDA(:,:)
  double precision,allocatable  :: FcGGA(:,:)
  double precision              :: aC

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas)

! Memory allocation

  select case (rung)

!   Hartree calculation

    case(0) 

      Fc(:,:) = 0d0

!   LDA functionals

    case(1) 

      call lda_correlation_potential(DFA,nEns,wEns(:),nGrid,weight(:),nBas,AO(:,:),rho(:),Fc(:,:))

!   GGA functionals

    case(2) 

      call gga_correlation_potential(DFA,nEns,wEns(:),nGrid,weight(:),nBas,AO(:,:),dAO(:,:,:),rho(:),drho(:,:),Fc(:,:))

!   Hybrid functionals

    case(4) 

      allocate(FcLDA(nBas,nBas),FcGGA(nBas,nBas))

      aC = 0.81d0

      call lda_correlation_potential(DFA,nEns,wEns(:),nGrid,weight(:),nBas,AO(:,:),rho(:),FcLDA(:,:))
      call gga_correlation_potential(DFA,nEns,wEns(:),nGrid,weight(:),nBas,AO(:,:),dAO(:,:,:),rho(:),drho(:,:),FcGGA(:,:))

      Fc(:,:) = FcLDA(:,:) + aC*(FcGGA(:,:) - FcLDA(:,:))  

!   Hartree-Fock calculation

    case(666) 

      Fc(:,:) = 0d0

  end select

end subroutine restricted_correlation_potential
