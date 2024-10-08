
! ---

subroutine orthogonalization_matrix(nBas, S, X)

! Compute the orthogonalization matrix X

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: S(nBas,nBas)

! Local variables

  logical                       :: debug
  double precision,allocatable  :: UVec(:,:),Uval(:)
  double precision,parameter    :: thresh = 1d-6
  integer,parameter             :: ortho_type = 1

  integer                       :: i

! Output variables

  double precision,intent(out)  :: X(nBas,nBas)

  debug = .false.

! Type of orthogonalization ortho_type
!
!  1 = Lowdin
!  2 = Canonical
!  3 = SVD
!

  allocate(Uvec(nBas,nBas),Uval(nBas))

  if(ortho_type == 1) then

    !
    ! S V = V s   where 
    !
    !     V.T V = 1 and s > 0 (S is positive def)
    !
    ! S = V s V.T
    !   = V s^0.5 s^0.5 V.T
    !   = V s^0.5 V.T V s^0.5 V.T
    !   = S^0.5 S^0.5 
    !
    ! where
    !
    !     S^0.5 = V s^0.5 V.T
    !
    ! X = S^(-0.5)
    !   = V s^(-0.5) V.T
    !

!   write(*,*)
!   write(*,*) ' Lowdin orthogonalization'
!   write(*,*)

    Uvec = S
    call diagonalize_matrix(nBas, Uvec, Uval)

    do i = 1, nBas

      if(Uval(i) < thresh) then 

        write(*,*) 'Eigenvalue',i,' is very small in Lowdin orthogonalization = ',Uval(i)

      end if

      Uval(i) = 1d0 / dsqrt(Uval(i))

    end do
    
    call ADAt(nBas, Uvec(1,1), Uval(1), X(1,1))

  elseif(ortho_type == 2) then

!   write(*,*)
!   write(*,*) 'Canonical orthogonalization'
!   write(*,*)

    Uvec = S
    call diagonalize_matrix(nBas, Uvec, Uval)

    do i = 1, nBas

      if(Uval(i) > thresh) then 

        Uval(i) = 1d0 / dsqrt(Uval(i))

      else

        write(*,*) ' Eigenvalue',i,'too small for canonical orthogonalization'

      end if

    end do
    
    call AD(nBas, Uvec, Uval)
    X = Uvec
 
  elseif(ortho_type == 3) then

!   write(*,*)
!   write(*,*) ' SVD-based orthogonalization NYI'
!   write(*,*)
   
!   Uvec = S
!   call diagonalize_matrix(nBas,Uvec,Uval)

!   do i=1,nBas
!     if(Uval(i) > thresh) then 
!       Uval(i) = 1d0/sqrt(Uval(i))
!     else
!       write(*,*) 'Eigenvalue',i,'too small for canonical orthogonalization'
!     end if
!   end do
!   
!   call AD(nBas,Uvec,Uval)
!   X = Uvec
 
  end if

! Print results

  if(debug) then

    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Orthogonalization matrix'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,X)
    write(*,*)

  end if

end subroutine 

! ---

