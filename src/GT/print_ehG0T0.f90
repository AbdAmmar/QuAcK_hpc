subroutine print_ehG0T0(nBas,nO,eHF,ENuc,ERHF,SigC,Z,eGT,EcRPA,EcGM)

! Print one-electron energies and other stuff for G0W0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: EcGM
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: SigC(nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: eGT(nBas)

  integer                            :: p,HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGT(LUMO)-eGT(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)'  One-shot eh G0T0 calculation'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_c (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',SigC(p)*HaToeV,'|',Z(p),'|',eGT(p)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A30,F15.6,A3)') 'ehG0T0 HOMO      energy:',eGT(HOMO)*HaToeV,' eV'
  write(*,'(2X,A30,F15.6,A3)') 'ehG0T0 LUMO      energy:',eGT(LUMO)*HaToeV,' eV'
  write(*,'(2X,A30,F15.6,A3)') 'ehG0T0 HOMO-LUMO gap   :',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A30,F15.6,A3)') 'RPA@ehG0T0 total energy      :',ENuc + ERHF + EcRPA,' au'
  write(*,'(2X,A30,F15.6,A3)') 'RPA@ehG0T0 correlation energy:',EcRPA,' au'
  write(*,'(2X,A30,F15.6,A3)') 'GM@ehG0T0  total energy      :',ENuc + ERHF + EcGM,' au'
  write(*,'(2X,A30,F15.6,A3)') 'GM@ehG0T0  correlation energy:',EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine
